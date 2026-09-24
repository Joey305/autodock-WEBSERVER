#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
import shlex
import sys
from datetime import datetime
from pathlib import Path
from typing import List

from hpc_profiles import packaged_profile_or_default, render_lsf_header, render_setup_block, replace_profile
from lsf_templates import COMPACTED_SDF_HTML_BODY, CSV_RESOLVER_BODY, PYMOL_BODY

HERE = Path(__file__).resolve().parent
DEFAULT_PROFILE = packaged_profile_or_default(HERE)
DEFAULT_ENV = render_setup_block(DEFAULT_PROFILE).strip()


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", Path(name).stem)


def parse_index_list(text: str, n: int) -> list[int]:
    text = text.strip()
    if not text:
        return []
    picks: set[int] = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            if a.strip().isdigit() and b.strip().isdigit():
                lo, hi = int(a), int(b)
                if lo > hi:
                    lo, hi = hi, lo
                for value in range(lo, hi + 1):
                    if 1 <= value <= n:
                        picks.add(value - 1)
        elif part.isdigit():
            value = int(part)
            if 1 <= value <= n:
                picks.add(value - 1)
    return sorted(picks)


def input_default(prompt: str, default: str) -> str:
    value = input(f"{prompt} [{default}]: ").strip()
    return value if value else str(default)


def discover_score_csvs(base: Path = Path(".")) -> list[Path]:
    patterns = [
        "ALL_Docking_Results_with_provenance_*.csv",
        "*vina_docking_scores_sorted.csv",
    ]
    seen: set[Path] = set()
    out: list[Path] = []
    for pattern in patterns:
        for path in sorted(base.glob(pattern), key=lambda p: p.stat().st_mtime, reverse=True):
            if not path.is_file() or path.name.endswith("_UNSORTED.tmp.csv") or ".tmp." in path.name:
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            out.append(path)
    return out


def show_indexed(items: list[Path], title: str) -> None:
    print(f"\n{title}")
    for i, item in enumerate(items, start=1):
        print(f" [{i}] {item.name}")


def _display_path(path: Path, base_dir: Path) -> str:
    try:
        return str(path.resolve().relative_to(base_dir.resolve()))
    except ValueError:
        return str(path)


def discover_receptor_root_candidates(csv_path: Path) -> list[Path]:
    search_bases = [Path.cwd(), csv_path.resolve().parent]
    candidates: list[Path] = []
    seen: set[Path] = set()

    for base in search_bases:
        if not base.exists():
            continue
        for path in sorted(base.iterdir(), key=lambda item: item.name.lower()):
            if not path.is_dir() or not path.name.startswith("Receptor"):
                continue
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            candidates.append(resolved)

    return candidates


def prompt_receptor_roots(csv_path: Path) -> list[str]:
    base_dir = csv_path.resolve().parent
    candidates = discover_receptor_root_candidates(csv_path)
    if not candidates:
        roots_default = "Receptors Receptors_PDBQT ."
        receptor_roots_raw = input_default("Receptor search roots (space-separated)", roots_default)
        return receptor_roots_raw.split()

    print("\nSelect receptor search folder(s):")
    for index, path in enumerate(candidates, start=1):
        print(f" [{index}] {_display_path(path, base_dir)}")
    print("Enter one or more indexes, e.g. 1,2, or type folder paths directly.")
    print("Press Enter to use all detected Receptor* folders.")

    while True:
        raw = input("Receptor search roots: ").strip()
        if not raw:
            return [_display_path(path, base_dir) for path in candidates]

        selected: list[str] = []
        bad_tokens: list[str] = []
        for token in raw.replace(",", " ").split():
            if token.isdigit():
                index = int(token)
                if 1 <= index <= len(candidates):
                    selected.append(_display_path(candidates[index - 1], base_dir))
                else:
                    bad_tokens.append(token)
                continue

            typed_path = Path(token).expanduser()
            if not typed_path.is_absolute():
                typed_path = (base_dir / typed_path).resolve()
            if typed_path.exists() and typed_path.is_dir():
                selected.append(_display_path(typed_path, base_dir))
            else:
                bad_tokens.append(token)

        if selected and not bad_tokens:
            return list(dict.fromkeys(selected))
        if bad_tokens:
            print(f"Could not resolve: {', '.join(bad_tokens)}")
        print("Please choose valid indexes from the list or type existing folder paths.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build LSF jobs for 5C_BuildPymolSesh.py and 5_COMPACTED_SDF_HTML.py."
    )
    parser.add_argument("--csvs", default="", help="Comma-separated score CSVs. Omit for interactive selection.")
    parser.add_argument("--auto", action="store_true", help="Use all discovered score CSVs.")
    parser.add_argument(
        "--tasks",
        default="both",
        choices=["both", "pymol", "compacted_sdf_html"],
        help="Which postprocessing jobs to write.",
    )
    parser.add_argument("--receptor-roots", default="", help="Space-separated receptor search roots.")
    parser.add_argument("--outdir", default=".", help="Output root for both builders.")
    parser.add_argument("--obabel-bin", default="", help="Path to obabel executable.")
    parser.add_argument("--hydrogen-mode", choices=["none", "nonpolar", "all"], default="none")
    parser.add_argument("--allow-reference-fallback", action="store_true")
    parser.add_argument("--include-previous-outputs-as-reference", action="store_true")
    parser.add_argument("--min-mcs-fraction", type=float, default=0.80)

    parser.add_argument("--pymol-mode", choices=["per_ligand", "per_receptor"], default="per_receptor")
    parser.add_argument("--pymol-top", type=int, default=5)
    parser.add_argument("--pymol-dry-run-selection", action="store_true")

    parser.add_argument("--top-ligands", type=int, default=25)
    parser.add_argument("--top-poses", type=int, default=25)
    parser.add_argument("--project-name", default="Docking_HTML_Viz_Project_Compacted_SDF")
    parser.add_argument("--page-title", default="Compacted SDF Docking Visualization Project")

    parser.add_argument("--workers", type=int, default=DEFAULT_PROFILE.workers, help=f"CPU cores per job (default: {DEFAULT_PROFILE.workers})")
    parser.add_argument("--queue", default=DEFAULT_PROFILE.queue, help=f"LSF queue (default: {DEFAULT_PROFILE.queue})")
    parser.add_argument("--project", default=DEFAULT_PROFILE.project, help=f"LSF project (default: {DEFAULT_PROFILE.project})")
    parser.add_argument("--walltime", default=DEFAULT_PROFILE.vina_walltime, help=f"Walltime (default: {DEFAULT_PROFILE.vina_walltime})")
    parser.add_argument("--mem-per-core", default=str(DEFAULT_PROFILE.mem_per_core_mb), help=f"MB per core for rusage[mem=...] (default: {DEFAULT_PROFILE.mem_per_core_mb})")
    parser.add_argument("--email", default=DEFAULT_PROFILE.email, help="Email for LSF notifications.")
    parser.add_argument("--env-activate", default=DEFAULT_ENV, help="Shell line(s) to activate the software environment.")
    return parser.parse_args()


def selected_csvs(args) -> list[Path]:
    if args.csvs.strip():
        csvs = [Path(item.strip()) for item in args.csvs.split(",") if item.strip()]
        missing = [str(path) for path in csvs if not path.is_file()]
        if missing:
            print("Missing CSV file(s): " + ", ".join(missing))
            sys.exit(2)
        return csvs

    candidates = discover_score_csvs()
    if args.auto:
        return candidates

    if not candidates:
        print("No docking score CSV files were found.")
        sys.exit(2)

    show_indexed(candidates, "Choose score CSV(s):")
    raw = input("Enter indices (comma-separated / ranges): ").strip()
    idxs = parse_index_list(raw, len(candidates))
    if not idxs:
        print("No CSVs selected.")
        sys.exit(2)
    return [candidates[i] for i in idxs]


def extra_flags(args, *, compacted: bool) -> str:
    flags: list[str] = []
    if args.obabel_bin.strip():
        flags.extend(["--obabel-bin", args.obabel_bin.strip()])
    if args.allow_reference_fallback:
        flags.append("--allow-reference-fallback")
    if args.include_previous_outputs_as_reference:
        flags.append("--include-previous-outputs-as-reference")
    flags.extend(["--min-mcs-fraction", str(args.min_mcs_fraction)])
    if not compacted and args.pymol_dry_run_selection:
        flags.append("--dry-run-selection")
    return shlex.join(flags)


def base_lsf_text(args, *, jobname: str, log_prefix: str) -> str:
    profile = replace_profile(
        DEFAULT_PROFILE,
        queue=args.queue,
        project=args.project,
        workers=args.workers,
        vina_cpus=args.workers,
        mem_per_core_mb=args.mem_per_core,
        vina_walltime=args.walltime,
        email=args.email,
        setup_commands=args.env_activate,
    )
    return (
        f"#!/bin/bash\n# Auto-generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        + render_lsf_header(
            profile=profile,
            jobname=jobname,
            log_prefix=log_prefix,
            walltime=args.walltime,
            workers=args.workers,
            mem_per_core_mb=args.mem_per_core,
        ).split("\n", 1)[1]
        + render_setup_block(profile)
        + f'PYBIN="{profile.python_command}"\nif [ -z "$PYBIN" ]; then\n  echo "No Python command configured"; exit 127\nfi\necho "Using Python: $PYBIN"\n\n'
    )


def write_pymol_lsf(csv_path: Path, args) -> Path:
    tag = sanitize_name(csv_path.name)
    text = (
        base_lsf_text(args, jobname=f"pymol_{tag}", log_prefix=f"pymol_{tag}")
        + f'export CSV="{csv_path}"\n'
        + f'export RECEPTOR_ROOTS="{args.receptor_roots}"\n'
        + f'export PYMOL_OUTDIR="{args.outdir}"\n'
        + f'export PYMOL_MODE="{args.pymol_mode}"\n'
        + f'export PYMOL_TOP="{args.pymol_top}"\n'
        + f'export HYDROGEN_MODE="{args.hydrogen_mode}"\n'
        + f'export PYMOL_EXTRA_ARGS="{extra_flags(args, compacted=False)}"\n'
        + CSV_RESOLVER_BODY
        + PYMOL_BODY
    )
    out = HERE / f"run_pymol_{tag}.lsf"
    out.write_text(text)
    return out


def write_compacted_sdf_html_lsf(csv_path: Path, args) -> Path:
    tag = sanitize_name(csv_path.name)
    text = (
        base_lsf_text(args, jobname=f"compact_sdf_{tag}", log_prefix=f"compact_sdf_{tag}")
        + f'export CSV="{csv_path}"\n'
        + f'export RECEPTOR_ROOTS="{args.receptor_roots}"\n'
        + f'export COMPACTED_OUTDIR="{args.outdir}"\n'
        + f'export COMPACTED_TOP_LIGANDS="{args.top_ligands}"\n'
        + f'export COMPACTED_TOP_POSES="{args.top_poses}"\n'
        + f'export COMPACTED_PROJECT_NAME="{args.project_name}"\n'
        + f'export COMPACTED_PAGE_TITLE="{args.page_title}"\n'
        + f'export HYDROGEN_MODE="{args.hydrogen_mode}"\n'
        + f'export COMPACTED_EXTRA_ARGS="{extra_flags(args, compacted=True)}"\n'
        + CSV_RESOLVER_BODY
        + COMPACTED_SDF_HTML_BODY
    )
    out = HERE / f"run_compacted_sdf_html_{tag}.lsf"
    out.write_text(text)
    return out


def write_submitter(paths: list[Path]) -> Path:
    sh = HERE / "submit_all_outputs.sh"
    lines = [
        "#!/bin/bash",
        "set -uo pipefail",
        'DIR="$(cd "$(dirname "$0")" && pwd)"',
        'cd "$DIR"',
        f'echo "Submitting {len(paths)} output builder job(s)..."',
        "fails=0",
    ]
    for path in paths:
        lines.append(f'echo "Submitting {path.name} ..."')
        lines.append(f'if ! bsub < "{path.name}"; then')
        lines.append(f'  echo "Failed: {path.name}"')
        lines.append("  fails=$((fails+1))")
        lines.append("fi")
    lines.append('echo "Done. Failed submissions: $fails"')
    sh.write_text("\n".join(lines) + "\n")
    os.chmod(sh, 0o755)
    return sh


def main():
    print("\n=== 5B_LSFbatch Output Builder LSF Generator ===")
    args = parse_args()
    csvs = selected_csvs(args)
    if not csvs:
        print("No CSVs selected.")
        sys.exit(2)

    if not args.csvs.strip() and not args.auto:
        args.tasks = input_default("Tasks (both/pymol/compacted_sdf_html)", args.tasks)
        if not args.receptor_roots.strip():
            args.receptor_roots = " ".join(prompt_receptor_roots(csvs[0]))
        args.outdir = input_default("Output root", args.outdir)
        args.pymol_top = int(input_default("PyMOL top N", str(args.pymol_top)))
        args.top_ligands = int(input_default("Compacted SDF top ligands", str(args.top_ligands)))
        args.top_poses = int(input_default("Compacted SDF top poses", str(args.top_poses)))
        args.queue = input_default("Queue", args.queue)
        args.project = input_default("Project", args.project)
        args.walltime = input_default("Walltime", args.walltime)
        args.workers = int(input_default("Workers", str(args.workers)))
        args.mem_per_core = input_default("Mem per core MB", str(args.mem_per_core))
        args.email = input_default("Email", args.email)

    if not args.receptor_roots.strip():
        args.receptor_roots = "Receptors Receptors_PDBQT ."

    paths: List[Path] = []
    for csv_path in csvs:
        if args.tasks in {"both", "pymol"}:
            path = write_pymol_lsf(csv_path, args)
            paths.append(path)
            print(f"Wrote {path.name}")
        if args.tasks in {"both", "compacted_sdf_html"}:
            path = write_compacted_sdf_html_lsf(csv_path, args)
            paths.append(path)
            print(f"Wrote {path.name}")

    submitter = write_submitter(paths)
    print(f"\nMaster submitter: {submitter.name}")
    print("Submit all with:")
    print("  ./submit_all_outputs.sh\n")


if __name__ == "__main__":
    main()
