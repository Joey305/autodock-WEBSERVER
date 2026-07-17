#!/usr/bin/env python3

import argparse, os, sys
import json
import csv
from pathlib import Path
from datetime import datetime
import re  # at top
from typing import Dict, List, Optional, Tuple

from hpc_profiles import packaged_profile_or_default, render_lsf_header, render_setup_block, replace_profile

HERE = Path(__file__).resolve().parent
DEFAULT_PROFILE = packaged_profile_or_default(HERE)
DEFAULT_ENV = render_setup_block(DEFAULT_PROFILE).strip()

SUBMITTER_HDR = """#!/bin/bash
set -euo pipefail
echo "Submitting {N} jobs..."
"""

def sanitize_name(p: str) -> str:
    # keep letters, digits, dot, underscore, dash; collapse everything else to "_"
    return re.sub(r'[^A-Za-z0-9._-]+', '_', Path(p).name)


def parse_args():
    p = argparse.ArgumentParser(
        description="Build per-target LSF jobs + a master submitter for 1_ConformerGeneration.py (CSV / folder of SDF|SMILES / single SDF / folder containing CSV)."
    )
    # Common knobs
    p.add_argument("--poses", type=int, default=64, help="Poses per molecule (default: 64)")
    p.add_argument("--workers", type=int, default=DEFAULT_PROFILE.workers, help=f"CPU cores per job (default: {DEFAULT_PROFILE.workers})")
    p.add_argument("--queue", default=DEFAULT_PROFILE.queue, help=f"LSF queue (default: {DEFAULT_PROFILE.queue})")
    p.add_argument("--project", default=DEFAULT_PROFILE.project, help=f"LSF project (default: {DEFAULT_PROFILE.project})")
    p.add_argument("--walltime", default=DEFAULT_PROFILE.confgen_walltime, help=f"Walltime (default: {DEFAULT_PROFILE.confgen_walltime})")
    p.add_argument("--mem-per-core", default=str(DEFAULT_PROFILE.mem_per_core_mb), help=f"MB per core for rusage[mem=...] (default: {DEFAULT_PROFILE.mem_per_core_mb})")
    p.add_argument("--email", default=DEFAULT_PROFILE.email, help="Email for notifications")
    p.add_argument("--env-activate", default=DEFAULT_ENV, help="Shell line to activate conda/env")
    p.add_argument("--obabel-bin", default="", help="Path to obabel (optional). Empty → PATH/OBABEL_BIN.")
    p.add_argument("--enumerate-protomers", action="store_true",
                   help="Pass through to 1_ConformerGeneration.py to enumerate protonation states with Dimorphite-DL.")
    p.add_argument("--ph-min", type=float, default=6.8,
                   help="Minimum pH for protomer enumeration (default: 6.8)")
    p.add_argument("--ph-max", type=float, default=7.4,
                   help="Maximum pH for protomer enumeration (default: 7.4)")
    p.add_argument("--ph-precision", type=float, default=0.5,
                   help="Dimorphite precision parameter (default: 0.5)")
    p.add_argument("--max-protomers", type=int, default=4,
                   help="Maximum protomer states to keep per ligand (default: 4)")
    p.add_argument("--max-tautomers", type=int, default=4,
                   help="Maximum tautomer states to keep per protomer (default: 4)")
    p.add_argument("--max-transforms", type=int, default=200,
                   help="RDKit tautomer transform cap (default: 200)")

    # Non-interactive fast-path (still supported)
    p.add_argument("--mode", choices=["1","2","3","4"], help="1=CSV, 2=Folder, 3=Single SDF, 4=Folder containing CSV")
    p.add_argument("--targets", default="", help="Comma-separated paths (CSV files for mode1; folders for mode2; SDF files for mode3; folders for mode4).")
    p.add_argument("--filetype", choices=["sdf","smiles"], help="Required for mode=2 to indicate folder contents.")
    p.add_argument("--csv-smiles-col", help="CSV SMILES column (mode=1)")
    p.add_argument("--csv-id-col", help="CSV ID column (mode=1)")
    p.add_argument("--csv-column-map", help="Optional JSON file mapping CSV folders/files to smiles/id column names for mode=4.")
    p.add_argument("--auto", action="store_true",
                   help="Auto-discover targets (mode=2 only): folders in the current directory containing ligand files for the selected --filetype.")
    return p.parse_args()

def discover_folders(filetype: str) -> List[str]:
    wanted_suffixes = {".sdf"} if filetype == "sdf" else {".smiles", ".smi"}
    out = []
    for d in sorted(Path(".").iterdir()):
        if not d.is_dir():
            continue
        if any(p.is_file() and p.suffix.lower() in wanted_suffixes for p in d.rglob("*")):
            out.append(d.name)
    return out

def discover_csv_folders() -> List[str]:
    out = []
    for d in sorted(Path(".").iterdir()):
        if not d.is_dir():
            continue
        if any(d.glob("*.csv")):
            out.append(d.name)
    return out

def read_csv_headers(csv_path: Path) -> List[str]:
    with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        return next(reader, []) or []

def resolve_single_csv_in_folder(folder: Path) -> Path:
    csv_files = sorted(folder.glob("*.csv"), key=lambda p: p.name.lower())
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {folder}")
    if len(csv_files) > 1:
        raise ValueError(f"Multiple CSV files found in {folder}; expected one")
    return csv_files[0]

def prompt_csv_columns(headers: List[str], target_label: str,
                       default_smiles: str = "", default_id: str = "") -> Tuple[str, str]:
    print(f"\nCSV columns for {target_label}:")
    for i, header in enumerate(headers):
        print(f" [{i}] {header}")

    header_lookup = {header.lower(): header for header in headers}
    guessed_smiles = default_smiles if default_smiles in headers else header_lookup.get("canonical_smiles") or header_lookup.get("smiles") or (headers[0] if headers else "")
    guessed_id = default_id if default_id in headers else header_lookup.get("molecule_chembl_id") or header_lookup.get("compound_id") or header_lookup.get("id") or (headers[1] if len(headers) > 1 else (headers[0] if headers else ""))

    def _resolve_column(raw_value: str, default_value: str) -> str:
        value = raw_value.strip()
        if not value:
            return default_value
        if value.isdigit():
            idx = int(value)
            if 0 <= idx < len(headers):
                return headers[idx]
        return value

    smiles_col = _resolve_column(input(f"Index or name of SMILES column? [{guessed_smiles}]: "), guessed_smiles)
    id_col = _resolve_column(input(f"Index or name of ID/name column? [{guessed_id}]: "), guessed_id)
    if smiles_col not in headers:
        print(f"❌ Column not found in {target_label}: {smiles_col}")
        sys.exit(2)
    if id_col not in headers:
        print(f"❌ Column not found in {target_label}: {id_col}")
        sys.exit(2)
    return smiles_col, id_col

def load_csv_column_map(path: str) -> Dict[str, Dict[str, str]]:
    payload = json.loads(Path(path).read_text())
    if not isinstance(payload, dict):
        raise ValueError("CSV column map JSON must be an object.")
    normalized: Dict[str, Dict[str, str]] = {}
    for key, value in payload.items():
        if not isinstance(value, dict):
            raise ValueError(f"CSV column map entry for {key} must be an object.")
        smiles_col = str(value.get("smiles_col", "")).strip()
        id_col = str(value.get("id_col", "")).strip()
        if not smiles_col or not id_col:
            raise ValueError(f"CSV column map entry for {key} must include smiles_col and id_col.")
        normalized[str(key)] = {"smiles_col": smiles_col, "id_col": id_col}
    return normalized

def csv_columns_for_target(target: str,
                           global_smiles_col: Optional[str],
                           global_id_col: Optional[str],
                           column_map: Optional[Dict[str, Dict[str, str]]] = None) -> Tuple[Optional[str], Optional[str]]:
    if column_map:
        target_path = Path(target)
        keys = [target, target_path.name, str(target_path.resolve())]
        for key in keys:
            entry = column_map.get(key)
            if entry:
                return entry["smiles_col"], entry["id_col"]
    return global_smiles_col, global_id_col

def list_with_indices(items: List[str], title: str):
    print(f"\n{title}")
    for i, name in enumerate(items, start=1):  # 1-based indices
        print(f" [{i}] {name}")

def parse_index_list(s: str, n: int) -> List[int]:
    """
    Parse '1,3,5-7' into [1,3,5,6,7]; 1-based -> return zero-based indices.
    """
    s = s.strip()
    if not s:
        return []
    picks = set()
    for part in s.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a, b = a.strip(), b.strip()
            if not a.isdigit() or not b.isdigit():
                continue
            lo, hi = int(a), int(b)
            if lo > hi: lo, hi = hi, lo
            for k in range(lo, hi + 1):
                if 1 <= k <= n:
                    picks.add(k - 1)
        else:
            if part.isdigit():
                k = int(part)
                if 1 <= k <= n:
                    picks.add(k - 1)
    return sorted(picks)

def build_state_flags(args) -> str:
    flags = [
        f" --ph-min {args.ph_min}",
        f" --ph-max {args.ph_max}",
        f" --ph-precision {args.ph_precision}",
        f" --max-protomers {args.max_protomers}",
        f" --max-tautomers {args.max_tautomers}",
        f" --max-transforms {args.max_transforms}",
    ]
    if args.enumerate_protomers:
        flags.append(" --enumerate-protomers")
    return "".join(flags)

def build_run_cmd(mode: str, target: str, poses: int, workers: int,
                  filetype: Optional[str], csv_smiles_col: Optional[str], csv_id_col: Optional[str],
                  args) -> str:
    state_flags = build_state_flags(args)
    if mode == "1":
        csv_flags = ""
        if csv_smiles_col:
            csv_flags += f" --smiles-col {csv_smiles_col}"
        if csv_id_col:
            csv_flags += f" --id-col {csv_id_col}"
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 1 --csv \"{target}\"{csv_flags} --num-confs {poses} --workers {workers}{state_flags}\n"
        )
    elif mode == "2":
        ft = filetype or "sdf"
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 2 --folder \"{target}\" --filetype {ft} --num-confs {poses} --workers {workers}{state_flags}\n"
        )
    elif mode == "3":
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 3 --sdf \"{target}\" --num-confs {poses} --workers {workers}{state_flags}\n"
        )
    else:
        csv_flags = ""
        if csv_smiles_col:
            csv_flags += f" --smiles-col {csv_smiles_col}"
        if csv_id_col:
            csv_flags += f" --id-col {csv_id_col}"
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 4 --folder \"{target}\"{csv_flags} --num-confs {poses} --workers {workers}{state_flags}\n"
        )

def write_lsf(name: str, run_cmd: str, args) -> Path:
    profile = replace_profile(
        DEFAULT_PROFILE,
        queue=args.queue,
        project=args.project,
        workers=args.workers,
        mem_per_core_mb=args.mem_per_core,
        confgen_walltime=args.walltime,
        email=args.email,
        setup_commands=args.env_activate,
    )
    obabel_line = f'export OBABEL_BIN="{args.obabel_bin}"\n' if args.obabel_bin.strip() else ""
    text = (
        f"#!/bin/bash\n# Auto-generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        + render_lsf_header(
            profile=profile,
            jobname=f"confgen_{name}",
            log_prefix=f"confgen_{name}",
            walltime=args.walltime,
            workers=args.workers,
            mem_per_core_mb=args.mem_per_core,
        ).split("\n", 1)[1]
        + render_setup_block(profile)
        + f'PYBIN="{profile.python_command}"\nif [ -z "$PYBIN" ]; then\n  echo "❌ No Python interpreter found after activating env"; exit 127\nfi\necho "Using Python: $PYBIN"\n\n'
        + obabel_line
        + run_cmd
    )
    out = HERE / f"run_confgen_{name}.lsf"
    out.write_text(text)
    return out

def write_submitter(paths: List[Path]) -> Path:
    sh = HERE / "submit_all_confgen.sh"
    lines = [SUBMITTER_HDR.format(N=len(paths))]
    for p in paths:
        lines.append(f'bsub < "{p.name}"')
    sh.write_text("\n".join(lines) + "\n")
    os.chmod(sh, 0o755)
    return sh

def prompt_mode_block():
    print("\nInput type?")
    print(" [1] CSV file(s)  (mode=1)")
    print(" [2] Directory(ies) of Ligands  (mode=2: SDF or SMILES)")
    print(" [3] Single SDF file(s) (mode=3)")
    print(" [4] Directory containing CSV file(s) (mode=4)")
    mode = input("Choose 1/2/3/4: ").strip()
    if mode not in ("1","2","3","4"):
        print("❌ Invalid choice."); sys.exit(2)

    targets = []
    filetype = None
    csv_smiles_col = None
    csv_id_col = None
    csv_column_map = None

    if mode == "1":
        print("\nRun scope?")
        print(" [1] One CSV")
        print(" [2] Multiple CSVs (comma-separated)")
        scope = input("Choose 1/2: ").strip()
        if scope == "1":
            t = input("CSV path: ").strip()
            targets = [t]
        else:
            t = input("CSV paths (comma-separated): ").strip()
            targets = [x.strip() for x in t.split(",") if x.strip()]
        csv_smiles_col = input("CSV SMILES column name: ").strip()
        csv_id_col = input("CSV ID column name: ").strip()

    elif mode == "2":
        # FOLDER MODE with index picker
        ft = input("Directory contents [sdf] or [smiles]? ").strip().lower()
        if ft not in ("sdf","smiles"):
            print("❌ Invalid filetype."); sys.exit(2)
        filetype = ft

        # Discover candidates and present index list
        candidates = discover_folders(filetype)
        if not candidates:
            print(f"❌ No directories found containing .{filetype} ligand files."); sys.exit(2)

        print("\nRun scope?")
        print(" [1] One Directory of Ligands")
        print(" [2] Multiple Directories of Ligands (choose by indices)")
        scope = input("Choose 1/2: ").strip()

        list_with_indices(candidates, title="Available Directories:")
        if scope == "1":
            s = input("Enter index of the directory to use: ").strip()
            idxs = parse_index_list(s, len(candidates))
            if len(idxs) != 1:
                print("❌ Please select exactly one index."); sys.exit(2)
            targets = [candidates[idxs[0]]]
        else:
            s = input("Enter comma-separated indices (supports ranges, e.g., 1,3,5-7): ").strip()
            idxs = parse_index_list(s, len(candidates))
            if not idxs:
                print("❌ No indices selected."); sys.exit(2)
            targets = [candidates[i] for i in idxs]

    elif mode == "3":
        print("\nRun scope?")
        print(" [1] One SDF file")
        print(" [2] Multiple SDF files (comma-separated)")
        scope = input("Choose 1/2: ").strip()
        if scope == "1":
            t = input("SDF file path: ").strip()
            targets = [t]
        else:
            t = input("SDF file paths (comma-separated): ").strip()
            targets = [x.strip() for x in t.split(",") if x.strip()]

    else:
        print("\nRun scope?")
        print(" [1] One folder containing a CSV")
        print(" [2] Multiple folders containing CSVs (choose by indices)")
        scope = input("Choose 1/2: ").strip()
        candidates = discover_csv_folders()
        if candidates:
            list_with_indices(candidates, title="Available CSV Folders:")
        if scope == "1":
            if candidates:
                s = input("Enter index of the folder to use, or type a path: ").strip()
                idxs = parse_index_list(s, len(candidates))
                targets = [candidates[idxs[0]]] if len(idxs) == 1 else [s]
            else:
                targets = [input("Folder path containing CSV: ").strip()]
        else:
            if candidates:
                s = input("Enter comma-separated indices (supports ranges, e.g., 1,3,5-7), or type paths separated by commas: ").strip()
                idxs = parse_index_list(s, len(candidates))
                targets = [candidates[i] for i in idxs] if idxs else [x.strip() for x in s.split(",") if x.strip()]
            else:
                s = input("Folder paths containing CSVs (comma-separated): ").strip()
                targets = [x.strip() for x in s.split(",") if x.strip()]
        if not targets:
            print("❌ No CSV folders selected."); sys.exit(2)

        target_paths = [Path(t) for t in targets]
        if len(target_paths) == 1:
            csv_path = resolve_single_csv_in_folder(target_paths[0])
            headers = read_csv_headers(csv_path)
            csv_smiles_col, csv_id_col = prompt_csv_columns(headers, target_paths[0].name)
        else:
            same_for_all = input("Use the same SMILES/ID columns for all selected CSV folders? [Y/n]: ").strip().lower()
            if same_for_all in ("", "y", "yes"):
                csv_path = resolve_single_csv_in_folder(target_paths[0])
                headers = read_csv_headers(csv_path)
                csv_smiles_col, csv_id_col = prompt_csv_columns(headers, target_paths[0].name)
            else:
                csv_column_map = {}
                for target_path in target_paths:
                    csv_path = resolve_single_csv_in_folder(target_path)
                    headers = read_csv_headers(csv_path)
                    smiles_col, id_col = prompt_csv_columns(headers, target_path.name)
                    csv_column_map[str(target_path)] = {"smiles_col": smiles_col, "id_col": id_col}

    return mode, targets, filetype, csv_smiles_col, csv_id_col, csv_column_map



def main():
    args = parse_args()

    # Non-interactive mode stays as before
    if args.mode and (args.targets or args.auto):
        mode = args.mode
        csv_column_map = load_csv_column_map(args.csv_column_map) if args.csv_column_map else None
        if mode == "2" and args.auto:
            if not args.filetype:
                print("❌ --auto requires --filetype for mode=2 (sdf|smiles)."); sys.exit(2)
            targets = discover_folders(args.filetype)
            filetype = args.filetype
        else:
            targets = [x.strip() for x in args.targets.split(",") if x.strip()]
            filetype = args.filetype
        csv_smiles_col = args.csv_smiles_col
        csv_id_col = args.csv_id_col
    else:
        # Interactive path with index picker for folder mode
        mode, targets, filetype, csv_smiles_col, csv_id_col, csv_column_map = prompt_mode_block()

    # Validate & filter
    valid = []
    for t in targets:
        p = Path(t)
        if mode == "1":
            if not p.is_file() or p.suffix.lower() != ".csv":
                print(f"⚠️ Skipping (not a CSV): {t}")
                continue
        elif mode == "2":
            if not p.is_dir():
                print(f"⚠️ Skipping (not a directory): {t}")
                continue
            if filetype == "sdf" and not any(x.is_file() and x.suffix.lower() == ".sdf" for x in p.rglob("*")):
                print(f"⚠️ Skipping (no .sdf in dir): {t}")
                continue
            if filetype == "smiles" and not any(
                x.is_file() and x.suffix.lower() in {".smiles", ".smi"} for x in p.rglob("*")
            ):
                print(f"⚠️ Skipping (no .smiles or .smi in dir): {t}")
                continue
        elif mode == "4":
            if not p.is_dir():
                print(f"⚠️ Skipping (not a directory): {t}")
                continue
            if not any(p.glob("*.csv")):
                print(f"⚠️ Skipping (no .csv in dir): {t}")
                continue
        else:
            if not p.is_file() or p.suffix.lower() != ".sdf":
                print(f"⚠️ Skipping (not an .sdf file): {t}")
                continue
        valid.append(t)

    if not valid:
        print("❌ No valid targets to build."); sys.exit(2)

    # Build per-target LSF and master submitter
    lsf_paths = []
    for t in valid:
        name = sanitize_name(t)
        target_smiles_col, target_id_col = csv_columns_for_target(t, csv_smiles_col, csv_id_col, csv_column_map)
        if mode in ("1", "4") and (not target_smiles_col or not target_id_col):
            print(f"❌ Missing CSV column settings for target: {t}")
            sys.exit(2)
        run_cmd = build_run_cmd(mode, t, args.poses, args.workers, filetype, target_smiles_col, target_id_col, args)
        lsf = write_lsf(name, run_cmd, args)
        lsf_paths.append(lsf)
        print(f"✅ Wrote {lsf.name}")

    submitter = write_submitter(lsf_paths)
    print(f"\n✅ Master submitter: {submitter.name}")
    print("Submit all with:")
    print("  ./submit_all_confgen.sh\n")

if __name__ == "__main__":
    main()
