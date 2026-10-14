#!/usr/bin/env python3
"""Create a resumable, memory-safe LSF array from interrupted Vina batches.

This utility is deliberately read-only with respect to the existing docking
results.  It trusts a docking only when its prior run log recorded success and
the corresponding out.pdbqt is present and non-empty; everything else is
placed in a new recovery shard.
"""
from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path


SUCCESS_RE = re.compile(r"^✅\s+(.+?)\s+→", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd(), help="IE1 docking directory")
    parser.add_argument("--recovery-root", type=Path, required=True, help="New directory for manifests, scripts, and results")
    parser.add_argument("--chunk-size", type=int, default=1000, help="Ligands per array task")
    parser.add_argument("--workers", type=int, default=7, help="Concurrent Vina processes per task")
    parser.add_argument("--cpus", type=int, default=8, help="LSF CPU slots per task")
    parser.add_argument("--memory-mb-per-slot", type=int, default=8000, help="LSF rusage memory per slot")
    parser.add_argument("--array-limit", type=int, default=100, help="Maximum simultaneous array tasks")
    parser.add_argument("--queue", default="normal")
    parser.add_argument("--project", default="brd")
    parser.add_argument("--walltime", default="48:00")
    parser.add_argument("--email", default="jxs794@miami.edu")
    parser.add_argument("--seed", type=int, default=676245557)
    parser.add_argument("--compact-results", action="store_true", help="Write compact shard CSVs instead of per-ligand pose trees.")
    parser.add_argument("--compact-scratch-root", help="Optional scratch root; defaults to VINA_SCRATCH_ROOT or each shard results directory.")
    return parser.parse_args()


def parse_batch_name(log_path: Path, ligand_dirs: list[Path]) -> tuple[str, Path, Path]:
    """Return receptor name, ligand directory, and matching old result directory."""
    stem = log_path.stem.removeprefix("run_log_")
    receptor_prefix, _, _ = stem.partition("_chembl_phase4_approved_smiles_part_")
    if not receptor_prefix.startswith("Receptors_"):
        raise ValueError(f"Cannot determine receptor from {log_path.name}")
    receptor = receptor_prefix.removeprefix("Receptors_")
    ligand_dir = next((path for path in ligand_dirs if path.name in stem), None)
    if ligand_dir is None:
        raise ValueError(f"Cannot determine ligand directory from {log_path.name}")
    old_results = log_path.parent / f"Docking_Results_{stem}"
    if not old_results.is_dir():
        raise FileNotFoundError(f"Expected result directory is missing: {old_results}")
    return receptor, ligand_dir, old_results


def valid_prior_output(results_dir: Path, receptor: str, ligand_name: str) -> bool:
    output = results_dir / receptor / ligand_name / "out.pdbqt"
    return output.is_file() and output.stat().st_size > 0


def write_array_script(args: argparse.Namespace, task_count: int) -> None:
    recovery_root = args.recovery_root.resolve()
    script = recovery_root / "run_recovery_array.lsf"
    compact_arg = "  --compact-results \\\n" if args.compact_results else ""
    scratch_export = (
        f'export VINA_SCRATCH_ROOT="${{VINA_SCRATCH_ROOT:-{args.compact_scratch_root}}}"\n'
        if args.compact_results and args.compact_scratch_root else ""
    )
    script.write_text(
        f"""#!/bin/bash
# Auto-generated recovery array. Existing successful docking outputs are preserved.
#BSUB -J IE1_recovery[1-{task_count}]%{args.array_limit}
#BSUB -P {args.project}
#BSUB -o logs/IE1_recovery_%J.out
#BSUB -e logs/IE1_recovery_%J.err
#BSUB -W {args.walltime}
#BSUB -q {args.queue}
#BSUB -n {args.cpus}
#BSUB -R \"span[hosts=1]\"
#BSUB -R \"rusage[mem={args.memory_mb_per_slot}]\"
#BSUB -B
#BSUB -N
#BSUB -u {args.email}

set -euo pipefail
cd {recovery_root}
mkdir -p logs results
source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh
conda activate vina_env
export VINA_EXE=\"$HOME/miniconda3/envs/vina_env/bin/vina\"
{scratch_export}

line=$(sed -n \"${{LSB_JOBINDEX}}p\" tasks.tsv)
IFS=$'\\t' read -r receptor manifest result_dir <<< \"$line\"
if [ -z \"${{receptor:-}}\" ] || [ ! -f \"$manifest\" ]; then
  echo \"Invalid task index ${{LSB_JOBINDEX}}\" >&2
  exit 2
fi

python3 ../3_Complete_batch_docking.py \\
  --receptors ../Receptors \\
  --receptor-file \"../Receptors/${{receptor}}.pdbqt\" \\
  --ligands .. \\
  --ligand-list \"$manifest\" \\
  --results-dir \"$result_dir\" \\
  --centers_csv ../vina_centers.csv \\
  --poses 10 \\
  --exhaustiveness 8 \\
  --min-rmsd 1 \\
  --energy-range 100 \\
  --seed {args.seed} \\
{compact_arg}
  --workers {args.workers}
""",
        encoding="utf-8",
    )
    script.chmod(0o755)
    submit = recovery_root / "submit_recovery_array.sh"
    submit.write_text("#!/bin/bash\nset -euo pipefail\nbsub < run_recovery_array.lsf\n", encoding="utf-8")
    submit.chmod(0o755)


def main() -> None:
    args = parse_args()
    if args.chunk_size < 1 or args.workers < 1 or args.cpus < args.workers:
        raise SystemExit("chunk size/workers must be positive, and cpus must be at least workers")

    root = args.root.resolve()
    args.recovery_root = args.recovery_root.resolve()
    if args.recovery_root.exists() and any(args.recovery_root.iterdir()):
        raise SystemExit(f"Recovery directory must be new or empty: {args.recovery_root}")
    ligand_dirs = sorted(root.glob("chembl_phase4_approved_smiles_part_*_Ligands_PDBQT_*"))
    logs = sorted(root.glob("run_log_Receptors_*.txt"))
    if len(ligand_dirs) != 15 or len(logs) != 60:
        raise SystemExit(f"Expected 15 ligand directories and 60 run logs; found {len(ligand_dirs)} and {len(logs)}")

    missing_by_receptor: dict[str, list[Path]] = defaultdict(list)
    retained = 0
    for log_path in logs:
        receptor, ligand_dir, old_results = parse_batch_name(log_path, ligand_dirs)
        prior_successes = set(SUCCESS_RE.findall(log_path.read_text(encoding="utf-8", errors="replace")))
        for ligand_path in sorted(ligand_dir.glob("*.pdbqt")):
            ligand_name = ligand_path.stem
            if ligand_name in prior_successes and valid_prior_output(old_results, receptor, ligand_name):
                retained += 1
            else:
                missing_by_receptor[receptor].append(ligand_path.resolve())

    args.recovery_root.mkdir(parents=True, exist_ok=True)
    manifests = args.recovery_root / "manifests"
    manifests.mkdir()
    tasks: list[tuple[str, Path, Path]] = []
    for receptor in sorted(missing_by_receptor):
        paths = missing_by_receptor[receptor]
        for index in range(0, len(paths), args.chunk_size):
            shard_number = index // args.chunk_size + 1
            manifest = manifests / f"{receptor}_shard_{shard_number:04d}.txt"
            manifest.write_text("\n".join(str(path) for path in paths[index:index + args.chunk_size]) + "\n", encoding="utf-8")
            results_dir = args.recovery_root / "results" / receptor / f"shard_{shard_number:04d}"
            tasks.append((receptor, manifest, results_dir))

    tasks_path = args.recovery_root / "tasks.tsv"
    tasks_path.write_text(
        "".join(f"{receptor}\t{manifest}\t{results_dir}\n" for receptor, manifest, results_dir in tasks),
        encoding="utf-8",
    )
    write_array_script(args, len(tasks))
    expected_missing = sum(len(paths) for paths in missing_by_receptor.values())
    print(f"Retained successful dockings: {retained}")
    print(f"Dockings to recover: {expected_missing}")
    print(f"Array tasks: {len(tasks)} (up to {args.array_limit} simultaneous, {args.workers} Vina workers/task)")
    for receptor in sorted(missing_by_receptor):
        print(f"  {receptor}: {len(missing_by_receptor[receptor])} dockings, {math.ceil(len(missing_by_receptor[receptor]) / args.chunk_size)} shards")


if __name__ == "__main__":
    main()
