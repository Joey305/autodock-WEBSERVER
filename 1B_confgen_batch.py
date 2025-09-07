#!/usr/bin/env python3
import argparse, os, sys
from pathlib import Path
from datetime import datetime
import re  # at top


HERE = Path(__file__).resolve().parent
DEFAULT_ENV = "source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh && conda activate vina_env"

LSF_TPL = """#!/bin/bash
# Auto-generated: {timestamp}
#BSUB -J confgen_{name}
#BSUB -P {project}
#BSUB -o logs/confgen_{name}_%J.out
#BSUB -e logs/confgen_{name}_%J.err
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
mkdir -p logs

# Activate environment
{env_activate}

# Pick a Python interpreter from the (now-activated) env
PYBIN="$(command -v python3 || command -v python)"
if [ -z "$PYBIN" ]; then
  echo "❌ No Python interpreter found after activating env"; exit 127
fi
echo "Using Python: $PYBIN"


# Optional obabel override (leave empty to rely on PATH/OBABEL_BIN)
{obabel_line}
{run_cmd}
"""

SUBMITTER_HDR = """#!/bin/bash
set -euo pipefail
echo "Submitting {N} jobs..."
"""

def sanitize_name(p: str) -> str:
    # keep letters, digits, dot, underscore, dash; collapse everything else to "_"
    return re.sub(r'[^A-Za-z0-9._-]+', '_', Path(p).name)


def parse_args():
    p = argparse.ArgumentParser(
        description="Build per-target LSF jobs + a master submitter for 1_ConformerGeneration.py (CSV / folder of SDF|SMILES / single SDF)."
    )
    # Common knobs
    p.add_argument("--poses", type=int, default=64, help="Poses per molecule (default: 64)")
    p.add_argument("--workers", type=int, default=16, help="CPU cores per job (default: 16)")
    p.add_argument("--queue", default="hihg", help="LSF queue (default: hihg)")
    p.add_argument("--project", default="brd", help="LSF project (default: brd)")
    p.add_argument("--walltime", default="48:00", help="Walltime (default: 48:00)")
    p.add_argument("--mem-per-core", default="2000", help="MB per core for rusage[mem=...] (default: 2000)")
    p.add_argument("--email", default="jxs794@miami.edu", help="Email for notifications")
    p.add_argument("--env-activate", default=DEFAULT_ENV, help="Shell line to activate conda/env")
    p.add_argument("--obabel-bin", default="", help="Path to obabel (optional). Empty → PATH/OBABEL_BIN.")

    # Non-interactive fast-path (still supported)
    p.add_argument("--mode", choices=["1","2","3"], help="1=CSV, 2=Folder, 3=Single SDF")
    p.add_argument("--targets", default="", help="Comma-separated paths (CSV files for mode1; folders for mode2; SDF files for mode3).")
    p.add_argument("--filetype", choices=["sdf","smiles"], help="Required for mode=2 to indicate folder contents.")
    p.add_argument("--csv-smiles-col", help="CSV SMILES column (mode=1)")
    p.add_argument("--csv-id-col", help="CSV ID column (mode=1)")
    p.add_argument("--auto", action="store_true",
                   help="Auto-discover targets (mode=2 only): folders matching Ligands_CPD* containing *.sdf or *.smiles per --filetype.")
    return p.parse_args()

def discover_folders(filetype: str) -> list[str]:
    out = []
    for d in sorted(Path(".").glob("Ligands_CPD*")):
        if not d.is_dir():
            continue
        if filetype == "sdf":
            if any(d.glob("*.sdf")):
                out.append(d.name)
        else:
            if any(d.glob("*.smiles")):
                out.append(d.name)
    return out

def list_with_indices(items: list[str], title: str):
    print(f"\n{title}")
    for i, name in enumerate(items, start=1):  # 1-based indices
        print(f" [{i}] {name}")

def parse_index_list(s: str, n: int) -> list[int]:
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

def build_run_cmd(mode: str, target: str, poses: int, workers: int,
                  filetype: str | None, csv_smiles_col: str | None, csv_id_col: str | None) -> str:
    if mode == "1":
        csv_flags = ""
        if csv_smiles_col:
            csv_flags += f" --smiles-col {csv_smiles_col}"
        if csv_id_col:
            csv_flags += f" --id-col {csv_id_col}"
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 1 --csv \"{target}\"{csv_flags} --num-confs {poses} --workers {workers}\n"
        )
    elif mode == "2":
        ft = filetype or "sdf"
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 2 --folder \"{target}\" --filetype {ft} --num-confs {poses} --workers {workers}\n"
        )
    else:
        return (
            '"$PYBIN" 1_ConformerGeneration.py'
            f" --mode 3 --sdf \"{target}\" --num-confs {poses} --workers {workers}\n"
        )

def write_lsf(name: str, run_cmd: str, args) -> Path:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    obabel_line = f'export OBABEL_BIN="{args.obabel_bin}"\n' if args.obabel_bin.strip() else ""
    text = LSF_TPL.replace("{run_cmd}", run_cmd).format(
        timestamp=ts,
        name=name,
        project=args.project,
        walltime=args.walltime,
        queue=args.queue,
        workers=args.workers,
        mem_per_core=args.mem_per_core,
        email=args.email,
        env_activate=args.env_activate,
        obabel_line=obabel_line
    )
    out = HERE / f"run_confgen_{name}.lsf"
    out.write_text(text)
    return out

def write_submitter(paths: list[Path]) -> Path:
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
    mode = input("Choose 1/2/3: ").strip()
    if mode not in ("1","2","3"):
        print("❌ Invalid choice."); sys.exit(2)

    targets = []
    filetype = None
    csv_smiles_col = None
    csv_id_col = None

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
            print(f"❌ No directories found matching Ligands_CPD* with .{filetype} files."); sys.exit(2)

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

    else:
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

    return mode, targets, filetype, csv_smiles_col, csv_id_col



def main():
    args = parse_args()

    # Non-interactive mode stays as before
    if args.mode and (args.targets or args.auto):
        mode = args.mode
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
        mode, targets, filetype, csv_smiles_col, csv_id_col = prompt_mode_block()

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
            if filetype == "sdf" and not any(p.glob("*.sdf")):
                print(f"⚠️ Skipping (no .sdf in dir): {t}")
                continue
            if filetype == "smiles" and not any(p.glob("*.smiles")):
                print(f"⚠️ Skipping (no .smiles in dir): {t}")
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
        run_cmd = build_run_cmd(mode, t, args.poses, args.workers, filetype, args.csv_smiles_col, args.csv_id_col)
        lsf = write_lsf(name, run_cmd, args)
        lsf_paths.append(lsf)
        print(f"✅ Wrote {lsf.name}")

    submitter = write_submitter(lsf_paths)
    print(f"\n✅ Master submitter: {submitter.name}")
    print("Submit all with:")
    print("  ./submit_all_confgen.sh\n")

if __name__ == "__main__":
    main()
