#!/usr/bin/env python3
from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from datetime import datetime

HERE = Path(__file__).resolve().parent

LSF_TEMPLATE = """#!/bin/bash
# Auto-generated: {timestamp}
#BSUB -J parse_{jobtag}
#BSUB -P {project}
#BSUB -o logs/parse_{jobtag}_%J.out
#BSUB -e logs/parse_{jobtag}_%J.err
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

source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh
conda activate vina_env

python 4_ParseScores.py \\
  --dir "{results_dir}" \\
  --workers {workers} \\
  --heartbeat {heartbeat} {fallback_flag}
"""

def sanitize_name(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s)

def parse_index_list(s: str, n: int) -> list[int]:
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

def result_dirs() -> list[Path]:
    return sorted(
        [p for p in Path(".").glob("Docking_Results_*") if p.is_dir()],
        key=lambda x: x.name.lower()
    )

def suffix_from_results_dir(results_dir: Path) -> str:
    prefix = "Docking_Results_"
    if not results_dir.name.startswith(prefix):
        raise ValueError(f"Unexpected results dir name: {results_dir.name}")
    return results_dir.name[len(prefix):]

def matching_run_log(results_dir: Path) -> Path | None:
    suffix = suffix_from_results_dir(results_dir)
    exact = results_dir.parent / f"run_log_{suffix}.txt"
    if exact.is_file():
        return exact
    cands = sorted(results_dir.parent.glob(f"run_log_{suffix}*.txt"))
    return cands[0] if cands else None

def has_job_duration(results_dir: Path) -> bool:
    suffix = suffix_from_results_dir(results_dir)
    return (results_dir.parent / f"job_duration_{suffix}.txt").is_file()

def has_docking_summary(results_dir: Path) -> bool:
    suffix = suffix_from_results_dir(results_dir)
    return (results_dir.parent / f"docking_summary_{suffix}.txt").is_file()

def has_parsed_csv(results_dir: Path) -> bool:
    patt = f"{results_dir.name}_*_vina_docking_scores_sorted.csv"
    return any(results_dir.parent.glob(patt))

def build_status_line(i: int, d: Path) -> str:
    flags = []
    if matching_run_log(d):
        flags.append("RUNLOG")
    if has_job_duration(d):
        flags.append("DURATION")
    if has_docking_summary(d):
        flags.append("SUMMARY")
    if has_parsed_csv(d):
        flags.append("PARSED")
    if not flags:
        flags.append("NO_MARKERS")
    return f" [{i}] {d.name}   [{' | '.join(flags)}]"

def write_lsf(results_dir: Path, queue: str, project: str, walltime: str,
              workers: int, mem_per_core: int, email: str,
              heartbeat: int, fallback_crawl: bool) -> Path:
    suffix = suffix_from_results_dir(results_dir)
    jobtag = sanitize_name(suffix)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    fallback_flag = "--fallback-crawl" if fallback_crawl else ""

    txt = LSF_TEMPLATE.format(
        timestamp=ts,
        jobtag=jobtag,
        project=project,
        walltime=walltime,
        queue=queue,
        workers=workers,
        mem_per_core=mem_per_core,
        email=email,
        results_dir=results_dir.name,
        heartbeat=heartbeat,
        fallback_flag=fallback_flag
    )

    out = HERE / f"run_parse_{jobtag}.lsf"
    out.write_text(txt)
    return out

def write_submitter(paths: list[Path]) -> Path:
    sh = HERE / "submit_all_PARSEscores.sh"
    lines = [
        "#!/bin/bash",
        "set -uo pipefail",
        'DIR="$(cd "$(dirname "$0")" && pwd)"',
        'cd "$DIR"',
        f'echo "Submitting {len(paths)} parse jobs..."',
        "fails=0"
    ]
    for p in paths:
        lines.append(f'echo "Submitting {p.name} ..."')
        lines.append(f'if ! bsub < "{p.name}"; then')
        lines.append(f'  echo "❌ Failed: {p.name}"')
        lines.append("  fails=$((fails+1))")
        lines.append("fi")
    lines.append('echo "Done. Failed submissions: $fails"')
    sh.write_text("\n".join(lines) + "\n")
    os.chmod(sh, 0o755)
    return sh

def input_default(prompt: str, default: str) -> str:
    s = input(f"{prompt} [{default}]: ").strip()
    return s if s else str(default)

def main():
    print("\n=== 4B_LSFbatch Parse LSF Builder (Interactive) ===")

    dirs = result_dirs()
    if not dirs:
        print("❌ No Docking_Results_* directories found.")
        sys.exit(2)

    print("\nDetected result directories:")
    print("  RUNLOG   = matching run_log_* exists")
    print("  DURATION = matching job_duration_* exists")
    print("  SUMMARY  = matching docking_summary_* exists")
    print("  PARSED   = at least one parsed CSV already exists\n")

    for i, d in enumerate(dirs, start=1):
        print(build_status_line(i, d))

    print("\nSelection mode?")
    print(" [0] ALL directories")
    print(" [1] Only directories with RUNLOG")
    print(" [2] Only directories without PARSED CSV")
    print(" [3] Only directories with RUNLOG and without PARSED CSV")
    print(" [4] Custom indices")
    mode = input("Choose 0/1/2/3/4 [0]: ").strip() or "0"

    if mode == "0":
        selected = dirs
    elif mode == "1":
        selected = [d for d in dirs if matching_run_log(d)]
    elif mode == "2":
        selected = [d for d in dirs if not has_parsed_csv(d)]
    elif mode == "3":
        selected = [d for d in dirs if matching_run_log(d) and not has_parsed_csv(d)]
    elif mode == "4":
        raw = input("Enter indices (e.g. 1,2,5-8): ").strip()
        idxs = parse_index_list(raw, len(dirs))
        selected = [dirs[i] for i in idxs]
    else:
        print("❌ Invalid choice.")
        sys.exit(2)

    if not selected:
        print("❌ No directories selected.")
        sys.exit(2)

    print("\nChosen directories:")
    for d in selected:
        print("  -", d.name)

    use_fallback = input("\nAdd --fallback-crawl to parse jobs? [y/N]: ").strip().lower() == "y"
    skip_existing = input("Skip directories that already have parsed CSVs? [y/N]: ").strip().lower() == "y"

    if skip_existing:
        selected = [d for d in selected if not has_parsed_csv(d)]

    if not selected:
        print("❌ Nothing left after filtering.")
        sys.exit(2)

    queue = input_default("Queue", "general")
    project = input_default("Project", "brd")
    walltime = input_default("Walltime", "200:00")
    workers = int(input_default("Workers", "16"))
    mem_per_core = int(input_default("Mem per core MB", "2000"))
    heartbeat = int(input_default("Heartbeat", "1000"))
    email = input_default("Email", "jxs794@miami.edu")

    out_paths = []
    for d in selected:
        p = write_lsf(
            results_dir=d,
            queue=queue,
            project=project,
            walltime=walltime,
            workers=workers,
            mem_per_core=mem_per_core,
            email=email,
            heartbeat=heartbeat,
            fallback_crawl=use_fallback
        )
        out_paths.append(p)
        print(f"✅ Wrote {p.name}")

    sub = write_submitter(out_paths)
    print(f"\n✅ Master submitter: {sub.name}")
    print("Submit all with:")
    print(f"  cd {HERE}")
    print("  ./submit_all_PARSEscores.sh\n")

if __name__ == "__main__":
    main()
