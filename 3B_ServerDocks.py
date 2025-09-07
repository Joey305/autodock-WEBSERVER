#!/usr/bin/env python3
import os, sys
from pathlib import Path
from datetime import datetime
import re

HERE = Path(__file__).resolve().parent

# ---------- utilities ----------
def count_pdbqt(d: Path) -> int:
    return sum(1 for _ in d.glob("*.pdbqt"))

def list_dirs_all() -> list[Path]:
    return sorted([p for p in Path(".").iterdir() if p.is_dir()], key=lambda x: x.name.lower())

def prioritize_dirs_for_ligands(dirs: list[Path]) -> list[Path]:
    def key(p: Path):
        has_pdbqt = count_pdbqt(p) > 0
        lig_prefix = p.name.startswith("Ligands_")
        # lower key sorts earlier
        return (
            0 if (lig_prefix and has_pdbqt) else
            1 if lig_prefix else
            2 if has_pdbqt else
            3,
            p.name.lower()
        )
    return sorted(dirs, key=key)

def prioritize_dirs_for_receptors(dirs: list[Path]) -> list[Path]:
    def key(p: Path):
        has_pdbqt = count_pdbqt(p) > 0
        rec_prefix = p.name.startswith("Receptors")
        return (
            0 if (rec_prefix and has_pdbqt) else
            1 if rec_prefix else
            2 if has_pdbqt else
            3,
            p.name.lower()
        )
    return sorted(dirs, key=key)

def list_files(pattern: str) -> list[Path]:
    return sorted([p for p in Path(".").glob(pattern) if p.is_file()], key=lambda x: x.name.lower())

def list_with_indices(items: list[Path], title: str, show_pdbqt_counts: bool = True):
    print(f"\n{title}")
    for i, p in enumerate(items, start=1):
        suffix = ""
        if show_pdbqt_counts:
            n = count_pdbqt(p)
            suffix = f" ({n} pdbqt)" if n else " (0 pdbqt)"
        print(f" [{i}] {p.name}{suffix}")

def parse_index_list(s: str, n: int) -> list[int]:
    s = s.strip()
    if not s: return []
    picks = set()
    for part in s.split(","):
        part = part.strip()
        if not part: continue
        if "-" in part:
            a, b = part.split("-", 1)
#!/usr/bin/env python3
import sys, os, re
from pathlib import Path
from datetime import datetime

HERE = Path(__file__).resolve().parent

# ---------- utilities ----------
def count_pdbqt(d: Path) -> int:
    return sum(1 for _ in d.glob("*.pdbqt"))

def all_dirs_sorted() -> list[Path]:
    return sorted([p for p in Path(".").iterdir() if p.is_dir()], key=lambda x: x.name.lower())

def prioritize_for_ligands(dirs: list[Path]) -> list[Path]:
    def key(p: Path):
        has_pdbqt = count_pdbqt(p) > 0
        lig_prefix = p.name.startswith("Ligands_")
        return (
            0 if (lig_prefix and has_pdbqt) else
            1 if lig_prefix else
            2 if has_pdbqt else
            3,
            p.name.lower()
        )
    return sorted(dirs, key=key)

def prioritize_for_receptors(dirs: list[Path]) -> list[Path]:
    def key(p: Path):
        has_pdbqt = count_pdbqt(p) > 0
        rec_prefix = p.name.startswith("Receptors")
        return (
            0 if (rec_prefix and has_pdbqt) else
            1 if rec_prefix else
            2 if has_pdbqt else
            3,
            p.name.lower()
        )
    return sorted(dirs, key=key)

def list_files(pattern: str) -> list[Path]:
    return sorted([p for p in Path(".").glob(pattern) if p.is_file()], key=lambda x: x.name.lower())

def show_indexed(items: list[Path], title: str, show_pdbqt_counts: bool = True):
    print(f"\n{title}")
    for i, p in enumerate(items, start=1):
        suffix = ""
        if show_pdbqt_counts:
            n = count_pdbqt(p)
            suffix = f" ({n} pdbqt)" if n else " (0 pdbqt)"
        print(f" [{i}] {p.name}{suffix}")

def parse_index_list(s: str, n: int) -> list[int]:
    s = s.strip()
    if not s: return []
    picks = set()
    for part in s.split(","):
        part = part.strip()
        if not part: continue
        if "-" in part:
            a, b = part.split("-", 1)
            if a.strip().isdigit() and b.strip().isdigit():
                lo, hi = int(a), int(b)
                if lo > hi: lo, hi = hi, lo
                for k in range(lo, hi+1):
                    if 1 <= k <= n: picks.add(k-1)
        else:
            if part.isdigit():
                k = int(part)
                if 1 <= k <= n: picks.add(k-1)
    return sorted(picks)

def input_default(prompt, default):
    s = input(f"{prompt} [{default}]: ").strip()
    return s if s else str(default)

def ligand_tag(name: str) -> str:
    m = re.match(r"^Ligands_CPD(\d+)_Ligands", name)
    return f"CPD{m.group(1)}" if m else name

def receptor_tag(name: str) -> str:
    return name[len("Receptors_"):] if name.startswith("Receptors_") else name

LSF_TEMPLATE = """#!/bin/bash
# Auto-generated: {timestamp}
#BSUB -J vina_{jobtag}
#BSUB -P {project}
#BSUB -o logs/vina_{jobtag}_%J.out
#BSUB -e logs/vina_{jobtag}_%J.err
#BSUB -W {walltime}
#BSUB -q {queue}
#BSUB -n {workers}
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem={mem_per_core}]"
#BSUB -B
#BSUB -N
#BSUB -u {email}

set -euo pipefail

cd "$LS_SUBCWD"
echo "PWD: $(pwd)"
mkdir -p logs

# Conda env with vina
source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh
conda activate vina_env

# Optional pin
export VINA_EXE="$HOME/miniconda3/envs/vina_env/bin/vina"

python 3_Complete_batch_docking.py \\
  --receptors "{receptors}" \\
  --ligands   "{ligands}" \\
  --centers_csv "{centers_csv}" \\
  --poses {poses}
"""

def write_lsf(jobtag: str, receptors: str, ligands: str, centers_csv: str,
              poses: int, queue: str, project: str, walltime: str,
              workers: int, mem_per_core: int, email: str) -> Path:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    txt = LSF_TEMPLATE.format(
        timestamp=ts, jobtag=jobtag, project=project, walltime=walltime,
        queue=queue, workers=workers, mem_per_core=mem_per_core, email=email,
        receptors=receptors, ligands=ligands, centers_csv=centers_csv, poses=poses
    )
    out = HERE / f"run_vina_{jobtag}.lsf"
    out.write_text(txt)
    return out

def write_submitter(paths: list[Path]) -> Path:
    sh = HERE / "submit_all_vina.sh"
    lines = ['#!/bin/bash', 'set -euo pipefail', f'echo "Submitting {len(paths)} docking jobs..."']
    for p in paths:
        lines.append(f'bsub < "{p.name}"')
    sh.write_text("\n".join(lines) + "\n")
    os.chmod(sh, 0o755)
    return sh

# ---------- main ----------
def main():
    print("\n=== Vina Job Builder ===")

    # Ask FIRST, then list
    lig_scope = input("\nLigands scope?\n [1] One Directory of Ligands\n [2] Multiple Directories of Ligands\nChoose 1/2: ").strip()
    all_d = all_dirs_sorted()
    lig_candidates = prioritize_for_ligands(all_d)

    if not lig_candidates:
        print("❌ No directories found."); sys.exit(2)

    if lig_scope == "1":
        show_indexed(lig_candidates, "Choose ONE Ligand Directory (all dirs; prioritized):")
        s = input("Enter index: ").strip()
        idxs = parse_index_list(s, len(lig_candidates))
        if len(idxs) != 1:
            print("❌ Please select exactly one index."); sys.exit(2)
        lig_dirs = [lig_candidates[idxs[0]]]
    else:
        show_indexed(lig_candidates, "Choose MULTIPLE Ligand Directories (all dirs; prioritized):")
        s = input("Enter comma-separated indices (ranges ok, e.g., 1,3,5-7): ").strip()
        idxs = parse_index_list(s, len(lig_candidates))
        if not idxs:
            print("❌ No indices selected."); sys.exit(2)
        lig_dirs = [lig_candidates[i] for i in idxs]

    # Receptors: ask mode FIRST, then list
    mode = input("\nReceptors input?\n [1] One directory with MULTIPLE receptor .pdbqt files\n [2] MULTIPLE directories each with ONE receptor .pdbqt\nChoose 1/2: ").strip()

    rec_candidates = prioritize_for_receptors(all_d)
    csv_candidates = list_files("*.csv")
    if not csv_candidates:
        print("❌ No CSV files found in CWD."); sys.exit(2)

    poses = int(input_default("\nDocking poses per ligand", 9))

    # Defaults for LSF
    queue = "hihg"
    project = "brd"
    walltime = "120:00"
    workers = 16
    mem_per_core = 2000
    email = "jxs794@miami.edu"

    out_paths = []

    if mode == "1":
        show_indexed(rec_candidates, "Select the Receptors Directory (all dirs; prioritized):")
        s = input("Enter index: ").strip()
        idxs = parse_index_list(s, len(rec_candidates))
        if len(idxs) != 1:
            print("❌ Select exactly one directory."); sys.exit(2)
        receptors_dir = rec_candidates[idxs[0]]

        # centers for this directory
        print("\nCenters CSV options:")
        for i, p in enumerate(csv_candidates, start=1):
            print(f" [{i}] {p.name}")
        s = input("Enter index of centers CSV for this receptors directory: ").strip()
        idxs = parse_index_list(s, len(csv_candidates))
        if len(idxs) != 1:
            print("❌ Select exactly one centers CSV."); sys.exit(2)
        centers_csv = csv_candidates[idxs[0]].name

        rec_tag = receptor_tag(receptors_dir.name)
        for lig_dir in lig_dirs:
            lt = ligand_tag(lig_dir.name)
            jobtag = f"{rec_tag}_{lt}"
            p = write_lsf(
                jobtag=jobtag, receptors=receptors_dir.name, ligands=lig_dir.name,
                centers_csv=centers_csv, poses=poses,
                queue=queue, project=project, walltime=walltime,
                workers=workers, mem_per_core=mem_per_core, email=email
            )
            out_paths.append(p)
            print(f"✅ Wrote {p.name}")

    else:
        show_indexed(rec_candidates, "Select Receptor Directories (all dirs; prioritized):")
        s = input("Enter indices (comma-separated / ranges): ").strip()
        idxs = parse_index_list(s, len(rec_candidates))
        if not idxs:
            print("❌ No receptor directories selected."); sys.exit(2)
        chosen_recs = [rec_candidates[i] for i in idxs]

        # centers per receptor directory
        rec_to_cent = {}
        for rd in chosen_recs:
            print(f"\nCenters CSV for receptor directory: {rd.name}")
            for i, p in enumerate(csv_candidates, start=1):
                print(f" [{i}] {p.name}")
            s = input("Enter index: ").strip()
            sel = parse_index_list(s, len(csv_candidates))
            if len(sel) != 1:
                print("❌ Select exactly one centers CSV."); sys.exit(2)
            rec_to_cent[rd] = csv_candidates[sel[0]].name

        for rd in chosen_recs:
            rec_tag = receptor_tag(rd.name)
            centers_csv = rec_to_cent[rd]
            for lig_dir in lig_dirs:
                lt = ligand_tag(lig_dir.name)
                jobtag = f"{rec_tag}_{lt}"
                p = write_lsf(
                    jobtag=jobtag, receptors=rd.name, ligands=lig_dir.name,
                    centers_csv=centers_csv, poses=poses,
                    queue=queue, project=project, walltime=walltime,
                    workers=workers, mem_per_core=mem_per_core, email=email
                )
                out_paths.append(p)
                print(f"✅ Wrote {p.name}")

    sub = write_submitter(out_paths)
    print(f"\n✅ Master submitter: {sub.name}")
    print("Submit all with:\n  ./submit_all_vina.sh\n")

if __name__ == "__main__":
    main()
