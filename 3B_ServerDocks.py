#!/usr/bin/env python3
from __future__ import annotations

import sys
import os
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List

from hpc_profiles import packaged_profile_or_default, render_lsf_header, render_setup_block, replace_profile

HERE = Path(__file__).resolve().parent
DEFAULT_PROFILE = packaged_profile_or_default(HERE)

def _ignore_candidate_dir(d: Path) -> bool:
    name = d.name
    return (
        not name
        or name.startswith(".")
        or name == "__pycache__"
        or name.startswith("logs")
        or name.startswith("Ligands_TMP_SDF_")
    )

# ---------- utilities ----------
def all_dirs_sorted() -> List[Path]:
    return sorted([p for p in Path(".").iterdir() if p.is_dir()], key=lambda x: x.name.lower())

def _ligand_dir_rank(d: Path):
    name = d.name.lower()
    return (
        0 if "_ligands_pdbqt_" in name else 1,
        0 if "ligand" in name else 1,
        0 if "batch" in name or "output" in name else 1,
        name,
    )

def _receptor_dir_rank(d: Path):
    name = d.name.lower()
    return (
        0 if name.startswith("receptors") else 1,
        0 if "receptor" in name else 1,
        name,
    )

def ligand_dirs_only() -> List[Path]:
    """
    Show directory names only, ranked so likely ligand output folders appear first.
    """
    return sorted(
        [
            p for p in Path(".").iterdir()
            if p.is_dir()
            and not _ignore_candidate_dir(p)
        ],
        key=_ligand_dir_rank,
    )

def receptor_dirs_only() -> List[Path]:
    """
    Show directory names only, ranked so likely receptor folders appear first.
    """
    return sorted(
        [
            p for p in Path(".").iterdir()
            if p.is_dir()
            and not _ignore_candidate_dir(p)
        ],
        key=_receptor_dir_rank,
    )

def list_files(pattern: str) -> List[Path]:
    return sorted([p for p in Path(".").glob(pattern) if p.is_file()], key=lambda x: x.name.lower())

def show_indexed(items: List[Path], title: str):
    print(f"\n{title}")
    for i, p in enumerate(items, start=1):
        print(f" [{i}] {p.name}")

def parse_index_list(s: str, n: int) -> List[int]:
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
            if a.strip().isdigit() and b.strip().isdigit():
                lo, hi = int(a), int(b)
                if lo > hi:
                    lo, hi = hi, lo
                for k in range(lo, hi + 1):
                    if 1 <= k <= n:
                        picks.add(k - 1)
        else:
            if part.isdigit():
                k = int(part)
                if 1 <= k <= n:
                    picks.add(k - 1)

    return sorted(picks)

def input_default(prompt: str, default):
    s = input(f"{prompt} [{default}]: ").strip()
    return s if s else str(default)

def ligand_tag(name: str) -> str:
    m = re.match(r"^Ligands_CPD(\d+)_Ligands", name)
    return f"CPD{m.group(1)}" if m else name

def receptor_tag(name: str) -> str:
    return name[len("Receptors_"):] if name.startswith("Receptors_") else name

def write_lsf(
    jobtag: str,
    receptors: str,
    ligands: str,
    centers_csv: str,
    poses: int,
    queue: str,
    project: str,
    walltime: str,
    workers: int,
    mem_per_core: int,
    email: str
) -> Path:
    profile = replace_profile(
        DEFAULT_PROFILE,
        queue=queue,
        project=project,
        workers=workers,
        mem_per_core_mb=mem_per_core,
        vina_walltime=walltime,
        email=email,
    )
    vina_pin = f'export VINA_EXE="{profile.vina_executable}"\n' if profile.vina_executable else ""
    txt = (
        f"#!/bin/bash\n# Auto-generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
        + render_lsf_header(
            profile=profile,
            jobname=f"vina_{jobtag}",
            log_prefix=f"vina_{jobtag}",
            walltime=walltime,
            workers=workers,
            mem_per_core_mb=mem_per_core,
        ).split("\n", 1)[1]
        + render_setup_block(profile)
        + vina_pin
        + f'PYBIN="{profile.python_command}"\nif [ -z "$PYBIN" ]; then\n  echo "No Python command configured"; exit 127\nfi\n'
        + '"$PYBIN" 3_Complete_batch_docking.py \\\n'
        + f'  --receptors "{receptors}" \\\n'
        + f'  --ligands   "{ligands}" \\\n'
        + f'  --centers_csv "{centers_csv}" \\\n'
        + f"  --poses {poses}\n"
    )
    out = HERE / f"run_vina_{jobtag}.lsf"
    out.write_text(txt)
    return out

def write_submitter(paths: list[Path]) -> Path:
    sh = HERE / "submit_all_vina.sh"
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        f'echo "Submitting {len(paths)} docking jobs..."'
    ]
    for p in paths:
        lines.append(f'bsub < "{p.name}"')
    sh.write_text("\n".join(lines) + "\n")
    os.chmod(sh, 0o755)
    return sh

# ---------- main ----------
def main():
    print("\n=== Vina Job Builder (Fast Mode) ===")

    lig_candidates = ligand_dirs_only()
    rec_candidates = receptor_dirs_only()
    csv_candidates = list_files("*.csv")

    if not lig_candidates:
        print("❌ No candidate ligand directories were found.")
        sys.exit(2)

    if not rec_candidates:
        print("❌ No candidate receptor directories were found.")
        sys.exit(2)

    if not csv_candidates:
        print("❌ No CSV files found in current working directory.")
        sys.exit(2)

    # Ask ligand scope
    lig_scope = input(
        "\nLigands scope?\n"
        " [1] One Directory of Ligands\n"
        " [2] Multiple Directories of Ligands\n"
        "Choose 1/2: "
    ).strip()

    if lig_scope == "1":
        show_indexed(lig_candidates, "Choose ONE Ligand Directory:")
        s = input("Enter index: ").strip()
        idxs = parse_index_list(s, len(lig_candidates))
        if len(idxs) != 1:
            print("❌ Please select exactly one ligand directory.")
            sys.exit(2)
        lig_dirs = [lig_candidates[idxs[0]]]
    elif lig_scope == "2":
        show_indexed(lig_candidates, "Choose MULTIPLE Ligand Directories:")
        s = input("Enter comma-separated indices (ranges ok, e.g. 1,3,5-7): ").strip()
        idxs = parse_index_list(s, len(lig_candidates))
        if not idxs:
            print("❌ No ligand directories selected.")
            sys.exit(2)
        lig_dirs = [lig_candidates[i] for i in idxs]
    else:
        print("❌ Invalid choice for ligand scope.")
        sys.exit(2)

    # Ask receptor mode
    mode = input(
        "\nReceptors input?\n"
        " [1] One directory with MULTIPLE receptor .pdbqt files\n"
        " [2] MULTIPLE directories each with ONE receptor .pdbqt\n"
        "Choose 1/2: "
    ).strip()

    poses = int(input_default("\nDocking poses per ligand", 9))

    # Defaults for LSF
    queue = DEFAULT_PROFILE.queue
    project = DEFAULT_PROFILE.project
    walltime = DEFAULT_PROFILE.vina_walltime
    workers = DEFAULT_PROFILE.workers
    mem_per_core = DEFAULT_PROFILE.mem_per_core_mb
    email = DEFAULT_PROFILE.email

    out_paths = []

    if mode == "1":
        show_indexed(rec_candidates, "Select the Receptors Directory:")
        s = input("Enter index: ").strip()
        idxs = parse_index_list(s, len(rec_candidates))
        if len(idxs) != 1:
            print("❌ Select exactly one receptor directory.")
            sys.exit(2)

        receptors_dir = rec_candidates[idxs[0]]

        print("\nCenters CSV options:")
        for i, p in enumerate(csv_candidates, start=1):
            print(f" [{i}] {p.name}")
        s = input("Enter index of centers CSV for this receptors directory: ").strip()
        idxs = parse_index_list(s, len(csv_candidates))
        if len(idxs) != 1:
            print("❌ Select exactly one centers CSV.")
            sys.exit(2)

        centers_csv = csv_candidates[idxs[0]].name
        rec_tag = receptor_tag(receptors_dir.name)

        for lig_dir in lig_dirs:
            lt = ligand_tag(lig_dir.name)
            jobtag = f"{rec_tag}_{lt}"
            p = write_lsf(
                jobtag=jobtag,
                receptors=receptors_dir.name,
                ligands=lig_dir.name,
                centers_csv=centers_csv,
                poses=poses,
                queue=queue,
                project=project,
                walltime=walltime,
                workers=workers,
                mem_per_core=mem_per_core,
                email=email
            )
            out_paths.append(p)
            print(f"✅ Wrote {p.name}")

    elif mode == "2":
        show_indexed(rec_candidates, "Select Receptor Directories:")
        s = input("Enter indices (comma-separated / ranges): ").strip()
        idxs = parse_index_list(s, len(rec_candidates))
        if not idxs:
            print("❌ No receptor directories selected.")
            sys.exit(2)

        chosen_recs = [rec_candidates[i] for i in idxs]

        rec_to_cent = {}
        for rd in chosen_recs:
            print(f"\nCenters CSV for receptor directory: {rd.name}")
            for i, p in enumerate(csv_candidates, start=1):
                print(f" [{i}] {p.name}")
            s = input("Enter index: ").strip()
            sel = parse_index_list(s, len(csv_candidates))
            if len(sel) != 1:
                print("❌ Select exactly one centers CSV.")
                sys.exit(2)
            rec_to_cent[rd] = csv_candidates[sel[0]].name

        for rd in chosen_recs:
            rec_tag = receptor_tag(rd.name)
            centers_csv = rec_to_cent[rd]

            for lig_dir in lig_dirs:
                lt = ligand_tag(lig_dir.name)
                jobtag = f"{rec_tag}_{lt}"
                p = write_lsf(
                    jobtag=jobtag,
                    receptors=rd.name,
                    ligands=lig_dir.name,
                    centers_csv=centers_csv,
                    poses=poses,
                    queue=queue,
                    project=project,
                    walltime=walltime,
                    workers=workers,
                    mem_per_core=mem_per_core,
                    email=email
                )
                out_paths.append(p)
                print(f"✅ Wrote {p.name}")
    else:
        print("❌ Invalid choice for receptor mode.")
        sys.exit(2)

    sub = write_submitter(out_paths)
    print(f"\n✅ Master submitter: {sub.name}")
    print("Submit all with:\n  ./submit_all_vina.sh\n")

if __name__ == "__main__":
    main()
