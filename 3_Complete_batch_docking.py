#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import csv
import glob
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import platform
from datetime import datetime
import argparse
import re
from pathlib import Path

from ligand_manifest import find_ligand_state_manifests, load_ligand_state_manifest, merge_ligand_metadata

RECEPTOR_SUFFIXES = (".pdbqt", ".pdb", ".mol2")

# ---------- UI helpers (only used if flags not provided) ----------
def choose(prompt, items):
    print(f"\n{prompt}", flush=True)
    for i, item in enumerate(items):
        print(f" [{i}] {item}", flush=True)
    idx = int(input("Enter index: "))
    return items[idx]

def progress_bar(done, total, width=40, successes=0, failures=0):
    """Render a single-line progress bar with percentage and counts."""
    if total == 0:
        total = 1
    ratio = done / total
    filled = int(ratio * width)
    bar = "#" * filled + "-" * (width - filled)
    pct = int(ratio * 100)
    msg = f"\r[{bar}] {pct:3d}%  ({done}/{total})  ✅ {successes}  ❌ {failures}"
    sys.stdout.write(msg)
    sys.stdout.flush()
    if done == total:
        sys.stdout.write("\n")
        sys.stdout.flush()

# ---------- naming helpers ----------
def _sanitize_tag(s: str) -> str:
    # Keep alnum, dash, underscore, dot; replace others with _
    return "".join(c if (c.isalnum() or c in "-_.") else "_" for c in s).strip("_") or "X"

def _receptor_tag_from_dir(receptor_dir: str) -> str:
    base = Path(receptor_dir).name
    if base.startswith("Receptors_"):
        base = base[len("Receptors_"):]
    return _sanitize_tag(base or "Receptors")

def _ligand_tag_from_dir(ligand_dir: str) -> str:
    base = Path(ligand_dir).name
    m = re.match(r"^Ligands_CPD(\d+)_Ligands", base)
    if m:
        return f"CPD{m.group(1)}"
    return _sanitize_tag(base or "Ligands")


def _strip_known_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in RECEPTOR_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _normalize_receptor_key(name: str) -> str:
    stem = _strip_known_suffix(os.path.basename(name))
    return re.sub(r"[^a-z0-9]+", "", stem.lower())


def _candidate_receptor_keys(name: str) -> set[str]:
    raw = _strip_known_suffix(os.path.basename(name))
    variants = {raw, raw.lower(), re.sub(r"[^A-Za-z0-9]+", "", raw)}
    base = re.split(r"[_-]+", raw, maxsplit=1)[0].strip()
    if base:
        variants.update({base, base.lower(), re.sub(r"[^A-Za-z0-9]+", "", base)})
    return {v for v in variants if v}


def _resolve_receptor_file(pdbid: str, receptor_files: list[str]) -> tuple[str | None, str]:
    requested = (pdbid or "").strip()
    if not requested:
        return None, "empty PDB_ID"

    exact_targets = [requested, *[requested + ext for ext in RECEPTOR_SUFFIXES]]
    by_basename = {os.path.basename(path): path for path in receptor_files}
    for target in exact_targets:
        match = by_basename.get(os.path.basename(target))
        if match:
            return match, "exact"

    requested_stem = _strip_known_suffix(requested).lower()
    stem_matches = [
        path for path in receptor_files
        if _strip_known_suffix(os.path.basename(path)).lower() == requested_stem
    ]
    if len(stem_matches) == 1:
        return stem_matches[0], "stem"

    requested_norm = _normalize_receptor_key(requested)
    norm_matches = [
        path for path in receptor_files
        if _normalize_receptor_key(os.path.basename(path)) == requested_norm
    ]
    if len(norm_matches) == 1:
        return norm_matches[0], "normalized"

    requested_keys = _candidate_receptor_keys(requested)
    fuzzy_matches = []
    for path in receptor_files:
        file_keys = _candidate_receptor_keys(os.path.basename(path))
        if any(
            req == got or req.startswith(got) or got.startswith(req)
            for req in requested_keys
            for got in file_keys
        ):
            fuzzy_matches.append(path)
    if len(fuzzy_matches) == 1:
        return fuzzy_matches[0], "fuzzy"

    return None, "missing" if not fuzzy_matches else "ambiguous"

# ---------- Docking worker ----------
def run_docking(job):
    vina_exe = job["vina_exe"]

    # Write Vina config
    with open(job["config_path"], "w") as cfg:
        cfg.write(f"receptor = {job['receptor_file']}\n")
        cfg.write(f"ligand = {job['ligand_file']}\n")
        cfg.write(f"center_x = {job['cx']}\n")
        cfg.write(f"center_y = {job['cy']}\n")
        cfg.write(f"center_z = {job['cz']}\n")
        cfg.write(f"size_x = 20\n")
        cfg.write(f"size_y = 20\n")
        cfg.write(f"size_z = 20\n")
        cfg.write(f"num_modes = {job['num_modes']}\n")
        cfg.write(f"out = {job['output_pdbqt']}\n")

    # Keep each Vina process to 1 thread; pool controls total concurrency
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"

    result = subprocess.run(
        [vina_exe, "--config", job["config_path"], "--cpu", "1"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env
    )

    # Save whichever stream has content
    with open(job["output_log"], "wb") as log_file:
        if result.stdout:
            log_file.write(result.stdout)
        if result.stderr:
            if result.stdout:
                log_file.write(b"\n--- STDERR ---\n")
            log_file.write(result.stderr)

    if result.returncode == 0:
        return True, f"✅ {job['ligand']} → {job['pdbid']}"
    else:
        err = result.stderr.decode(errors='ignore').strip()
        return False, f"❌ {job['ligand']} vs {job['pdbid']}: {err}"


def resolve_vina_executable(cli_vina_exe: str | None = None) -> str:
    candidates = []
    if cli_vina_exe:
        candidates.append(cli_vina_exe)
    env_vina = os.environ.get("VINA_EXE", "").strip()
    if env_vina:
        candidates.append(env_vina)
    path_vina = shutil.which("vina")
    if path_vina:
        candidates.append(path_vina)
    conda_prefix = os.environ.get("CONDA_PREFIX", "").strip()
    if conda_prefix:
        candidates.append(os.path.join(conda_prefix, "bin", "vina"))

    seen = set()
    for candidate in candidates:
        if not candidate:
            continue
        candidate = os.path.abspath(candidate) if os.path.sep in candidate else candidate
        if candidate in seen:
            continue
        seen.add(candidate)
        probe = candidate
        if probe == "vina":
            resolved = shutil.which("vina")
            probe = resolved if resolved else probe
        if probe != "vina" and not os.path.exists(probe):
            continue
        try:
            result = subprocess.run(
                [probe, "--version"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
                text=True,
            )
        except OSError:
            continue
        if result.returncode == 0:
            return os.path.abspath(probe) if probe != "vina" else probe

    raise SystemExit(
        "Cannot find AutoDock Vina executable.\n"
        "Install/download Vina so `vina --version` works, or set VINA_EXE=/path/to/vina, "
        "or pass --vina-exe /path/to/vina.\n"
        "Do not use `pip install vina` as the primary install path on systems where it "
        "tries to build from source and fails on Boost."
    )


def build_jobs(results_dir, ligand_dir, receptor_dir, grid_rows, num_modes, vina_exe):
    jobs = []
    ligand_files = glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
    receptor_files = []
    for suffix in RECEPTOR_SUFFIXES:
        receptor_files.extend(glob.glob(os.path.join(receptor_dir, f"*{suffix}")))
    ligand_names = [os.path.splitext(os.path.basename(f))[0] for f in ligand_files]
    manifest_map = {}
    for manifest_path in find_ligand_state_manifests([Path(ligand_dir), Path(ligand_dir).parent]):
        manifest_map.update(load_ligand_state_manifest(manifest_path))

    if not ligand_names:
        raise FileNotFoundError(f"No ligands (*.pdbqt) found in {ligand_dir}")
    if not receptor_files:
        raise FileNotFoundError(f"No receptors (*{', *'.join(RECEPTOR_SUFFIXES)}) found in {receptor_dir}")

    for row in grid_rows:
        pdbid = row["PDB_ID"]
        cx, cy, cz = row["X"], row["Y"], row["Z"]
        receptor_file, _ = _resolve_receptor_file(pdbid, receptor_files)
        if not receptor_file:
            continue

        for ligand in ligand_names:
            ligand_file = os.path.join(ligand_dir, f"{ligand}.pdbqt")
            if not os.path.exists(ligand_file):
                continue

            safe_pdbid = os.path.basename(receptor_file)
            safe_pdbid = safe_pdbid.replace(".pdbqt", "").replace(".pdb", "").replace(".mol2", "")
            output_subdir = os.path.join(results_dir, safe_pdbid, ligand)
            os.makedirs(output_subdir, exist_ok=True)

            jobs.append(
                {
                    "ligand": ligand,
                    "ligand_file": ligand_file,
                    "receptor_file": receptor_file,
                    "output_pdbqt": os.path.join(output_subdir, "out.pdbqt"),
                    "output_log": os.path.join(output_subdir, "log.txt"),
                    "config_path": os.path.join(output_subdir, "config.txt"),
                    "cx": cx,
                    "cy": cy,
                    "cz": cz,
                    "pdbid": pdbid,
                    "num_modes": num_modes,
                    "vina_exe": vina_exe,
                    "ligand_metadata": merge_ligand_metadata(ligand, manifest_row=manifest_map.get(ligand)),
                }
            )
    return jobs

def cpu_count_cgroup_aware():
    """Respect cgroup/LSF limits; fall back to os.cpu_count()."""
    try:
        n = len(os.sched_getaffinity(0))
    except Exception:
        n = os.cpu_count() or 1
    # LSF often sets this; honor it if present
    lsf_n = os.environ.get("LSB_DJOB_NUMPROC")
    if lsf_n:
        try:
            n = min(n, int(lsf_n))
        except ValueError:
            pass
    return max(1, n)

def parse_args():
    ap = argparse.ArgumentParser(
        description="Batch AutoDock Vina runner (scheduler-friendly). "
                    "Provide flags for non-interactive use; falls back to interactive otherwise."
    )
    ap.add_argument("--receptors", help="Path to receptor folder")
    ap.add_argument("--ligands", help="Path to ligand folder with .pdbqt ligands")
    ap.add_argument("--centers_csv", help="Path to vina centers CSV with headers PDB_ID,X,Y,Z")
    ap.add_argument("--poses", type=int, help="num_modes per ligand (e.g., 9, 20, 64)")
    ap.add_argument("--vina-exe", help="Path to AutoDock Vina executable")
    ap.add_argument("--reserve_cores", type=int, default=1,
                    help="How many cores to reserve for system/IO (default: 1)")
    return ap.parse_args()

# ---------- Main ----------
def main():
    args = parse_args()
    start_time = time.time()
    start_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    timestamp_tag = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    cwd = os.getcwd()
    vina_exe = resolve_vina_executable(args.vina_exe)

    # Resolve inputs (flags preferred; else interactive)
    if args.receptors and args.ligands and args.centers_csv and args.poses:
        receptor_dir = os.path.abspath(args.receptors)
        ligand_dir   = os.path.abspath(args.ligands)
        vina_csv     = os.path.abspath(args.centers_csv)
        num_modes    = int(args.poses)
    else:
        dirs = [d for d in os.listdir(cwd) if os.path.isdir(d)]
        receptor_dir = choose("📁 Select receptor folder:", dirs)
        ligand_dir   = choose("📁 Select ligand folder:", dirs)
        csv_files    = [f for f in os.listdir(cwd) if f.endswith(".csv")]
        vina_csv     = choose("📄 Select vina_centers.csv file:", csv_files)
        num_modes    = int(input("\n🔢 How many poses per ligand? (e.g., 9, 20, 64): "))

        receptor_dir = os.path.join(cwd, receptor_dir)
        ligand_dir   = os.path.join(cwd, ligand_dir)
        vina_csv     = os.path.join(cwd, vina_csv)

    # --- Read grid centers FIRST (so we can validate before building names/dirs) ---
    with open(vina_csv, "r", newline="") as f:
        reader = csv.DictReader(f)
        grid_rows = list(reader)

    required_cols = {"PDB_ID", "X", "Y", "Z"}
    if not required_cols.issubset(set(reader.fieldnames or [])):
        raise ValueError(f"CSV {vina_csv} must have headers: {sorted(required_cols)}")

    # --- Build naming tags from folder names only (no centers tag) ---
    rec_tag = _receptor_tag_from_dir(receptor_dir)
    lig_tag = _ligand_tag_from_dir(ligand_dir)

    # --- Outputs with receptor/ligand/poses/timestamp ---
    results_dir = os.path.join(
        cwd, f"Docking_Results_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}"
    )
    config_dir  = os.path.join(
        cwd, f"Docking_Configs_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}"
    )
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(config_dir, exist_ok=True)

    # Run log (append all results; keeps terminal quiet) — tagged
    run_log_path = os.path.join(
        cwd, f"run_log_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt"
    )
    run_log = open(run_log_path, "a", encoding="utf-8")

    jobs = build_jobs(results_dir, ligand_dir, receptor_dir, grid_rows, num_modes, vina_exe)
    receptor_files = []
    for suffix in RECEPTOR_SUFFIXES:
        receptor_files.extend(glob.glob(os.path.join(receptor_dir, f"*{suffix}")))
    missing_receptor_msgs = []
    for row in grid_rows:
        receptor_file, match_kind = _resolve_receptor_file(row["PDB_ID"], receptor_files)
        if not receptor_file:
            missing_receptor_msgs.append(
                f"❌ Unmatched receptor center: {row['PDB_ID']} (match={match_kind})"
            )
    for msg in missing_receptor_msgs:
        run_log.write(msg + "\n")
    run_log.flush()

    ncpus = cpu_count_cgroup_aware()
    reserve = max(0, int(args.reserve_cores))
    max_workers = max(1, ncpus - reserve)

    print(f"\n🧠 Detected {ncpus} schedulable cores. "
          f"Running up to {max_workers} concurrent Vina jobs (reserve {reserve}).", flush=True)
    print(f"🗂  Results dir: {results_dir}", flush=True)
    print(f"⚙️  Configs dir:  {config_dir}\n", flush=True)
    print(f"🧪 Vina executable: {vina_exe}", flush=True)

    successes, failures = 0, 0
    results_log = []
    total_jobs = len(jobs)

    if total_jobs == 0:
        available = ", ".join(sorted(os.path.basename(path) for path in receptor_files)[:10]) or "(none found)"
        details = "; ".join(missing_receptor_msgs[:5]) if missing_receptor_msgs else "no receptor-center rows matched"
        raise RuntimeError(
            "No jobs prepared. Check receptors/ligands/CSV inputs. "
            f"Available receptors: {available}. Details: {details}"
        )

    # Initial progress bar line
    progress_bar(0, total_jobs, successes=0, failures=0)

    # Concurrency limited by max_workers; quiet terminal, log to file
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_docking, job) for job in jobs]
        done_count = 0
        for future in as_completed(futures):
            ok, line = future.result()
            results_log.append(line)
            run_log.write(line + "\n")
            run_log.flush()

            if ok:
                successes += 1
            else:
                failures += 1

            done_count += 1
            progress_bar(done_count, total_jobs, successes=successes, failures=failures)

    end_time = time.time()
    duration_min = round((end_time - start_time) / 60, 2)
    duration_sec = round(end_time - start_time, 2)
    end_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_log.close()

    # Duration file in working directory (tagged)
    duration_txt = os.path.join(
        cwd, f"job_duration_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt"
    )
    with open(duration_txt, "w", encoding="utf-8") as df:
        df.write(f"Start Time : {start_stamp}\n")
        df.write(f"End Time   : {end_stamp}\n")
        df.write(f"Duration   : {duration_min} minutes\n")
        df.write(f"Duration_s : {duration_sec} seconds\n")
        df.write(f"Jobs Total : {total_jobs}\n")
        df.write(f"Successes  : {successes}\n")
        df.write(f"Failures   : {failures}\n")

    # Summary file (tagged)
    summary_path = os.path.join(
        cwd, f"docking_summary_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt"
    )
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Docking Summary Report\n")
        f.write("=======================\n")
        f.write(f"Start Time : {start_stamp}\n")
        f.write(f"End Time   : {end_stamp}\n")
        f.write(f"Total Time : {duration_min} minutes\n")
        f.write(f"Jobs Run   : {total_jobs}\n")
        f.write(f"Successes  : {successes}\n")
        f.write(f"Failures   : {failures}\n\n")
        f.write("System Info:\n")
        f.write(f" - Platform: {platform.system()} {platform.release()}\n")
        f.write(f" - CPU Cores (schedulable): {ncpus}\n")
        f.write(f" - Max Workers (cores used): {max_workers}\n")
        f.write(f" - Threads per Vina job: 1 (via --cpu 1 / OMP_NUM_THREADS=1)\n\n")
        f.write("Inputs/Outputs:\n")
        f.write(f" - Receptor dir: {receptor_dir}\n")
        f.write(f" - Ligand dir  : {ligand_dir}\n")
        f.write(f" - Centers CSV : {vina_csv}\n")
        f.write(f" - Results dir : {results_dir}\n")
        f.write(f" - Configs dir : {config_dir}\n")
        f.write(f" - num_modes   : {num_modes}\n\n")
        f.write("Results:\n--------\n")
        for line in results_log:
            f.write(line + "\n")

    print(f"\n⏱️  Duration file : {duration_txt}", flush=True)
    print(f"🧾 Run log       : {run_log_path}", flush=True)
    print(f"📄 Summary saved : {summary_path}", flush=True)
    print(f"✅ Done. Results in: {results_dir}", flush=True)

if __name__ == "__main__":
    main()












# #!/usr/bin/env python3
# import os
# import sys
# import csv
# import glob
# import subprocess
# from concurrent.futures import ProcessPoolExecutor, as_completed
# import time
# import platform
# from datetime import datetime
# import argparse

# # ---------- UI helpers (only used if flags not provided) ----------
# def choose(prompt, items):
#     print(f"\n{prompt}", flush=True)
#     for i, item in enumerate(items):
#         print(f" [{i}] {item}", flush=True)
#     idx = int(input("Enter index: "))
#     return items[idx]

# def progress_bar(done, total, width=40, successes=0, failures=0):
#     """Render a single-line progress bar with percentage and counts."""
#     if total == 0:
#         total = 1
#     ratio = done / total
#     filled = int(ratio * width)
#     bar = "#" * filled + "-" * (width - filled)
#     pct = int(ratio * 100)
#     msg = f"\r[{bar}] {pct:3d}%  ({done}/{total})  ✅ {successes}  ❌ {failures}"
#     sys.stdout.write(msg)
#     sys.stdout.flush()
#     if done == total:
#         sys.stdout.write("\n")
#         sys.stdout.flush()

# # ---------- Docking worker ----------
# def run_docking(job):
#     # Prefer $VINA_EXE, else fall back to "vina" in PATH (your conda env provides it)
#     vina_exe = os.environ.get("VINA_EXE", "vina")

#     # Write Vina config
#     with open(job["config_path"], "w") as cfg:
#         cfg.write(f"receptor = {job['receptor_file']}\n")
#         cfg.write(f"ligand = {job['ligand_file']}\n")
#         cfg.write(f"center_x = {job['cx']}\n")
#         cfg.write(f"center_y = {job['cy']}\n")
#         cfg.write(f"center_z = {job['cz']}\n")
#         cfg.write(f"size_x = 20\n")
#         cfg.write(f"size_y = 20\n")
#         cfg.write(f"size_z = 20\n")
#         cfg.write(f"num_modes = {job['num_modes']}\n")
#         cfg.write(f"out = {job['output_pdbqt']}\n")

#     # Keep each Vina process to 1 thread; pool controls total concurrency
#     env = os.environ.copy()
#     env["OMP_NUM_THREADS"] = "1"

#     result = subprocess.run(
#         [vina_exe, "--config", job["config_path"], "--cpu", "1"],
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         env=env
#     )

#     # Save whichever stream has content
#     with open(job["output_log"], "wb") as log_file:
#         if result.stdout:
#             log_file.write(result.stdout)
#         if result.stderr:
#             if result.stdout:
#                 log_file.write(b"\n--- STDERR ---\n")
#             log_file.write(result.stderr)

#     if result.returncode == 0:
#         return True, f"✅ {job['ligand']} → {job['pdbid']}"
#     else:
#         err = result.stderr.decode(errors='ignore').strip()
#         return False, f"❌ {job['ligand']} vs {job['pdbid']}: {err}"

# def cpu_count_cgroup_aware():
#     """Respect cgroup/LSF limits; fall back to os.cpu_count()."""
#     try:
#         n = len(os.sched_getaffinity(0))
#     except Exception:
#         n = os.cpu_count() or 1
#     # LSF often sets this; honor it if present
#     lsf_n = os.environ.get("LSB_DJOB_NUMPROC")
#     if lsf_n:
#         try:
#             n = min(n, int(lsf_n))
#         except ValueError:
#             pass
#     return max(1, n)

# def parse_args():
#     ap = argparse.ArgumentParser(
#         description="Batch AutoDock Vina runner (scheduler-friendly). "
#                     "Provide flags for non-interactive use; falls back to interactive otherwise."
#     )
#     ap.add_argument("--receptors", help="Path to receptor folder")
#     ap.add_argument("--ligands", help="Path to ligand folder with .pdbqt ligands")
#     ap.add_argument("--centers_csv", help="Path to vina centers CSV with headers PDB_ID,X,Y,Z")
#     ap.add_argument("--poses", type=int, help="num_modes per ligand (e.g., 9, 20, 64)")
#     ap.add_argument("--reserve_cores", type=int, default=1,
#                     help="How many cores to reserve for system/IO (default: 1)")
#     return ap.parse_args()

# # ---------- Main ----------
# def main():
#     args = parse_args()
#     start_time = time.time()
#     start_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     timestamp_tag = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
#     cwd = os.getcwd()

#     # Resolve inputs (flags preferred; else interactive)
#     if args.receptors and args.ligands and args.centers_csv and args.poses:
#         receptor_dir = os.path.abspath(args.receptors)
#         ligand_dir   = os.path.abspath(args.ligands)
#         vina_csv     = os.path.abspath(args.centers_csv)
#         num_modes    = int(args.poses)
#     else:
#         dirs = [d for d in os.listdir(cwd) if os.path.isdir(d)]
#         receptor_dir = choose("📁 Select receptor folder:", dirs)
#         ligand_dir   = choose("📁 Select ligand folder:", dirs)
#         csv_files    = [f for f in os.listdir(cwd) if f.endswith(".csv")]
#         vina_csv     = choose("📄 Select vina_centers.csv file:", csv_files)
#         num_modes    = int(input("\n🔢 How many poses per ligand? (e.g., 9, 20, 64): "))

#         receptor_dir = os.path.join(cwd, receptor_dir)
#         ligand_dir   = os.path.join(cwd, ligand_dir)
#         vina_csv     = os.path.join(cwd, vina_csv)

#     # Outputs with poses + timestamp
#     results_dir = os.path.join(cwd, f"Docking_Results_{num_modes}Poses_{timestamp_tag}")
#     config_dir  = os.path.join(cwd, f"Docking_Configs_{num_modes}Poses_{timestamp_tag}")
#     os.makedirs(results_dir, exist_ok=True)
#     os.makedirs(config_dir, exist_ok=True)

#     # Run log (append all results; keeps terminal quiet)
#     run_log_path = os.path.join(cwd, f"run_log_{timestamp_tag}.txt")
#     run_log = open(run_log_path, "a", encoding="utf-8")

#     # Read grid centers
#     with open(vina_csv, "r", newline="") as f:
#         reader = csv.DictReader(f)
#         grid_rows = list(reader)

#     # Validate headers
#     required_cols = {"PDB_ID", "X", "Y", "Z"}
#     if not required_cols.issubset(set(reader.fieldnames or [])):
#         raise ValueError(f"CSV {vina_csv} must have headers: {sorted(required_cols)}")

#     # Prepare jobs
#     jobs = []
#     ligand_files = glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
#     ligand_names = [os.path.splitext(os.path.basename(f))[0] for f in ligand_files]

#     if not ligand_names:
#         raise FileNotFoundError(f"No ligands (*.pdbqt) found in {ligand_dir}")

#     for row in grid_rows:
#         pdbid = row["PDB_ID"]
#         cx, cy, cz = row["X"], row["Y"], row["Z"]
#         receptor_file = os.path.join(receptor_dir, pdbid)

#         if not os.path.exists(receptor_file):
#             # try common extensions
#             alt = None
#             for ext in (".pdbqt", ".pdb", ".mol2"):
#                 cand = receptor_file + ext
#                 if os.path.exists(cand):
#                     alt = cand
#                     break
#             if alt:
#                 receptor_file = alt
#             else:
#                 msg = f"❌ Missing receptor: {receptor_file}(.pdbqt/.pdb/.mol2)"
#                 run_log.write(msg + "\n"); run_log.flush()
#                 continue

#         for ligand in ligand_names:
#             ligand_file = os.path.join(ligand_dir, f"{ligand}.pdbqt")
#             if not os.path.exists(ligand_file):
#                 msg = f"❌ Missing ligand: {ligand_file}"
#                 run_log.write(msg + "\n"); run_log.flush()
#                 continue

#             safe_pdbid = os.path.basename(receptor_file)
#             safe_pdbid = safe_pdbid.replace(".pdbqt", "").replace(".pdb", "").replace(".mol2", "")

#             output_subdir = os.path.join(results_dir, safe_pdbid, ligand)
#             os.makedirs(output_subdir, exist_ok=True)

#             output_pdbqt = os.path.join(output_subdir, "out.pdbqt")
#             output_log   = os.path.join(output_subdir, "log.txt")
#             config_path  = os.path.join(output_subdir, "config.txt")

#             jobs.append({
#                 "ligand": ligand,
#                 "ligand_file": ligand_file,
#                 "receptor_file": receptor_file,
#                 "output_pdbqt": output_pdbqt,
#                 "output_log": output_log,
#                 "config_path": config_path,
#                 "cx": cx,
#                 "cy": cy,
#                 "cz": cz,
#                 "pdbid": pdbid,
#                 "num_modes": num_modes
#             })

#     ncpus = cpu_count_cgroup_aware()
#     reserve = max(0, int(args.reserve_cores))
#     max_workers = max(1, ncpus - reserve)

#     print(f"\n🧠 Detected {ncpus} schedulable cores. "
#           f"Running up to {max_workers} concurrent Vina jobs (reserve {reserve}).", flush=True)
#     print(f"🗂  Results dir: {results_dir}", flush=True)
#     print(f"⚙️  Configs dir:  {config_dir}\n", flush=True)

#     successes, failures = 0, 0
#     results_log = []
#     total_jobs = len(jobs)

#     if total_jobs == 0:
#         raise RuntimeError("No jobs prepared. Check receptors/ligands/CSV inputs.")

#     # Initial progress bar line
#     progress_bar(0, total_jobs, successes=0, failures=0)

#     # Concurrency limited by max_workers; quiet terminal, log to file
#     with ProcessPoolExecutor(max_workers=max_workers) as executor:
#         futures = [executor.submit(run_docking, job) for job in jobs]
#         done_count = 0
#         for future in as_completed(futures):
#             ok, line = future.result()
#             results_log.append(line)
#             run_log.write(line + "\n")
#             run_log.flush()

#             if ok:
#                 successes += 1
#             else:
#                 failures += 1

#             done_count += 1
#             progress_bar(done_count, total_jobs, successes=successes, failures=failures)

#     end_time = time.time()
#     duration_min = round((end_time - start_time) / 60, 2)
#     duration_sec = round(end_time - start_time, 2)
#     end_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
#     run_log.close()

#     # Duration file in working directory
#     duration_txt = os.path.join(cwd, f"job_duration_{timestamp_tag}.txt")
#     with open(duration_txt, "w", encoding="utf-8") as df:
#         df.write(f"Start Time : {start_stamp}\n")
#         df.write(f"End Time   : {end_stamp}\n")
#         df.write(f"Duration   : {duration_min} minutes\n")
#         df.write(f"Duration_s : {duration_sec} seconds\n")
#         df.write(f"Jobs Total : {total_jobs}\n")
#         df.write(f"Successes  : {successes}\n")
#         df.write(f"Failures   : {failures}\n")

#     # Summary file (also quiet)
#     summary_path = os.path.join(cwd, f"docking_summary_{timestamp_tag}.txt")
#     with open(summary_path, "w", encoding="utf-8") as f:
#         f.write("Docking Summary Report\n")
#         f.write("=======================\n")
#         f.write(f"Start Time : {start_stamp}\n")
#         f.write(f"End Time   : {end_stamp}\n")
#         f.write(f"Total Time : {duration_min} minutes\n")
#         f.write(f"Jobs Run   : {total_jobs}\n")
#         f.write(f"Successes  : {successes}\n")
#         f.write(f"Failures   : {failures}\n\n")
#         f.write("System Info:\n")
#         f.write(f" - Platform: {platform.system()} {platform.release()}\n")
#         f.write(f" - CPU Cores (schedulable): {ncpus}\n")
#         f.write(f" - Max Workers (cores used): {max_workers}\n")
#         f.write(f" - Threads per Vina job: 1 (via --cpu 1 / OMP_NUM_THREADS=1)\n\n")
#         f.write("Inputs/Outputs:\n")
#         f.write(f" - Receptor dir: {receptor_dir}\n")
#         f.write(f" - Ligand dir  : {ligand_dir}\n")
#         f.write(f" - Centers CSV : {vina_csv}\n")
#         f.write(f" - Results dir : {results_dir}\n")
#         f.write(f" - Configs dir : {config_dir}\n")
#         f.write(f" - num_modes   : {num_modes}\n\n")
#         f.write("Results:\n--------\n")
#         for line in results_log:
#             f.write(line + "\n")

#     print(f"\n⏱️  Duration file : {duration_txt}", flush=True)
#     print(f"🧾 Run log       : {run_log_path}", flush=True)
#     print(f"📄 Summary saved : {summary_path}", flush=True)
#     print(f"✅ Done. Results in: {results_dir}", flush=True)

# if __name__ == "__main__":
#     main()
