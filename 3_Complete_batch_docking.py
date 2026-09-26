#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import platform
import re
import secrets
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

from ligand_manifest import (
    find_ligand_state_manifests,
    load_ligand_state_manifest,
    merge_ligand_metadata,
)

RECEPTOR_SUFFIXES = (".pdbqt", ".pdb", ".mol2")

# AutoDock Vina defaults. We write them explicitly so downstream analyses know
# exactly how the docking ensemble was produced.
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_MIN_RMSD = 1.0
DEFAULT_ENERGY_RANGE = 3.0

# A deliberately wide output window for downstream ensemble analyses such as
# 8_Protacability.py. This does NOT change Vina's search/scoring; it prevents
# the normal 3 kcal/mol output window from discarding retained modes merely
# because they score >3 kcal/mol above mode 1.
FULL_ENSEMBLE_ENERGY_RANGE = 100.0

# Keep generated seeds positive and safely inside a signed 32-bit integer.
MAX_VINA_SEED = 2_147_483_647


# ---------- UI helpers (only used if flags not provided) ----------
def choose(prompt, items):
    print(f"\n{prompt}", flush=True)
    for i, item in enumerate(items):
        print(f" [{i}] {item}", flush=True)
    idx = int(input("Enter index: "))
    return items[idx]


def input_default(prompt: str, default):
    value = input(f"{prompt} [{default}]: ").strip()
    return value if value else str(default)


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


# ---------- Vina parameter / seed helpers ----------
def _validate_vina_parameters(
    num_modes: int,
    exhaustiveness: int,
    min_rmsd: float,
    energy_range: float,
) -> None:
    if num_modes <= 0:
        raise ValueError("--poses must be greater than zero")
    if exhaustiveness <= 0:
        raise ValueError("--exhaustiveness must be greater than zero")
    if min_rmsd <= 0:
        raise ValueError("--min-rmsd must be greater than zero")
    if energy_range <= 0:
        raise ValueError("--energy-range must be greater than zero")


def _normalize_batch_seed(seed: int | None) -> int:
    """
    Return a positive reproducible batch seed.

    None or 0 means AUTO. The returned batch seed is recorded in provenance
    and used to derive a deterministic seed for every receptor-ligand job.
    """
    if seed in (None, 0):
        return secrets.randbelow(MAX_VINA_SEED - 1) + 1
    if seed < 1 or seed > MAX_VINA_SEED:
        raise ValueError(f"--seed must be between 1 and {MAX_VINA_SEED}, or 0/AUTO")
    return int(seed)


def _make_job_seed(
    batch_seed: int,
    receptor_key: str,
    ligand_key: str,
    cx: str,
    cy: str,
    cz: str,
) -> int:
    """
    Derive a stable positive seed for one docking job.

    Including receptor, ligand, and grid center means the same batch seed
    recreates exactly the same per-job seed assignment without giving every
    docking the identical random stream.
    """
    token = f"{batch_seed}|{receptor_key}|{ligand_key}|{cx}|{cy}|{cz}".encode("utf-8")
    digest = hashlib.sha256(token).digest()
    seed = int.from_bytes(digest[:4], "big") & 0x7FFFFFFF
    return seed or 1


def _pose_retention_mode(energy_range: float) -> str:
    if energy_range >= FULL_ENSEMBLE_ENERGY_RANGE:
        return "full_ensemble"
    if abs(energy_range - DEFAULT_ENERGY_RANGE) < 1e-12:
        return "standard_vina"
    return "custom_energy_range"


# ---------- Docking worker ----------
def run_docking(job):
    vina_exe = job["vina_exe"]

    # Write every relevant Vina setting explicitly. This makes each individual
    # docking directory self-describing and reproducible.
    with open(job["config_path"], "w", encoding="utf-8") as cfg:
        cfg.write(f"receptor = {job['receptor_file']}\n")
        cfg.write(f"ligand = {job['ligand_file']}\n")
        cfg.write(f"center_x = {job['cx']}\n")
        cfg.write(f"center_y = {job['cy']}\n")
        cfg.write(f"center_z = {job['cz']}\n")
        cfg.write("size_x = 20\n")
        cfg.write("size_y = 20\n")
        cfg.write("size_z = 20\n")
        cfg.write(f"num_modes = {job['num_modes']}\n")
        cfg.write(f"exhaustiveness = {job['exhaustiveness']}\n")
        cfg.write(f"min_rmsd = {job['min_rmsd']}\n")
        cfg.write(f"energy_range = {job['energy_range']}\n")
        cfg.write(f"seed = {job['seed']}\n")
        cfg.write(f"out = {job['output_pdbqt']}\n")

    # Keep each Vina process to 1 thread; the Python pool controls total concurrency.
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"

    result = subprocess.run(
        [vina_exe, "--config", job["config_path"], "--cpu", "1"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )

    # Save whichever stream has content.
    with open(job["output_log"], "wb") as log_file:
        if result.stdout:
            log_file.write(result.stdout)
        if result.stderr:
            if result.stdout:
                log_file.write(b"\n--- STDERR ---\n")
            log_file.write(result.stderr)

    if result.returncode == 0:
        return True, f"✅ {job['ligand']} → {job['pdbid']}  seed={job['seed']}"

    err = result.stderr.decode(errors="ignore").strip()
    return False, f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {err}"


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


def get_vina_version(vina_exe: str) -> str:
    try:
        result = subprocess.run(
            [vina_exe, "--version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )
    except OSError as exc:
        return f"unknown ({exc})"

    text = (result.stdout or result.stderr or "").strip()
    return text.splitlines()[0].strip() if text else "unknown"


def build_jobs(
    results_dir,
    ligand_dir,
    receptor_dir,
    grid_rows,
    num_modes,
    vina_exe,
    exhaustiveness,
    min_rmsd,
    energy_range,
    batch_seed,
    receptor_file=None,
):
    jobs = []
    ligand_files = glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
    receptor_files = []
    for suffix in RECEPTOR_SUFFIXES:
        receptor_files.extend(glob.glob(os.path.join(receptor_dir, f"*{suffix}")))
    all_receptor_files = receptor_files.copy()

    if receptor_file:
        selected_receptor = os.path.abspath(receptor_file)
        receptor_files = [
            path for path in receptor_files
            if os.path.abspath(path) == selected_receptor
        ]

    ligand_names = [os.path.splitext(os.path.basename(f))[0] for f in ligand_files]

    manifest_map = {}
    for manifest_path in find_ligand_state_manifests([Path(ligand_dir), Path(ligand_dir).parent]):
        manifest_map.update(load_ligand_state_manifest(manifest_path))

    if not ligand_names:
        raise FileNotFoundError(f"No ligands (*.pdbqt) found in {ligand_dir}")
    if not receptor_files:
        raise FileNotFoundError(
            f"No receptors (*{', *'.join(RECEPTOR_SUFFIXES)}) found in {receptor_dir}"
        )

    for row in grid_rows:
        pdbid = row["PDB_ID"]
        cx, cy, cz = row["X"], row["Y"], row["Z"]
        # Resolve against the complete directory first.  With only one file in
        # the candidate list, the legacy fuzzy matching can otherwise assign a
        # neighboring center row to the selected receptor.
        resolved_receptor, _ = _resolve_receptor_file(pdbid, all_receptor_files)
        if not resolved_receptor or (
            receptor_file and os.path.abspath(resolved_receptor) != selected_receptor
        ):
            continue

        safe_pdbid = os.path.basename(resolved_receptor)
        safe_pdbid = safe_pdbid.replace(".pdbqt", "").replace(".pdb", "").replace(".mol2", "")

        for ligand in ligand_names:
            ligand_file = os.path.join(ligand_dir, f"{ligand}.pdbqt")
            if not os.path.exists(ligand_file):
                continue

            output_subdir = os.path.join(results_dir, safe_pdbid, ligand)
            os.makedirs(output_subdir, exist_ok=True)

            job_seed = _make_job_seed(batch_seed, safe_pdbid, ligand, cx, cy, cz)

            jobs.append(
                {
                    "ligand": ligand,
                    "ligand_file": ligand_file,
                    "receptor_file": resolved_receptor,
                    "output_pdbqt": os.path.join(output_subdir, "out.pdbqt"),
                    "output_log": os.path.join(output_subdir, "log.txt"),
                    "config_path": os.path.join(output_subdir, "config.txt"),
                    "cx": cx,
                    "cy": cy,
                    "cz": cz,
                    "pdbid": pdbid,
                    "resolved_receptor": safe_pdbid,
                    "num_modes": num_modes,
                    "exhaustiveness": exhaustiveness,
                    "min_rmsd": min_rmsd,
                    "energy_range": energy_range,
                    "seed": job_seed,
                    "vina_exe": vina_exe,
                    "ligand_metadata": merge_ligand_metadata(
                        ligand,
                        manifest_row=manifest_map.get(ligand),
                    ),
                }
            )

    return jobs


def write_job_manifest(results_dir: str, jobs: list[dict]) -> str:
    """Write a compact machine-readable index for downstream Step 8 analysis."""
    manifest_path = os.path.join(results_dir, "docking_jobs.csv")
    fieldnames = [
        "PDB_ID",
        "ResolvedReceptor",
        "Ligand",
        "ReceptorFile",
        "LigandFile",
        "CenterX",
        "CenterY",
        "CenterZ",
        "NumModes",
        "Exhaustiveness",
        "MinRMSD",
        "EnergyRange",
        "Seed",
        "ConfigPath",
        "OutputPDBQT",
        "OutputLog",
        "LigandMetadataJSON",
    ]

    root = Path(results_dir).resolve()

    def rel(path: str) -> str:
        try:
            return str(Path(path).resolve().relative_to(root))
        except ValueError:
            return str(Path(path).resolve())

    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for job in jobs:
            writer.writerow(
                {
                    "PDB_ID": job["pdbid"],
                    "ResolvedReceptor": job["resolved_receptor"],
                    "Ligand": job["ligand"],
                    "ReceptorFile": str(Path(job["receptor_file"]).resolve()),
                    "LigandFile": str(Path(job["ligand_file"]).resolve()),
                    "CenterX": job["cx"],
                    "CenterY": job["cy"],
                    "CenterZ": job["cz"],
                    "NumModes": job["num_modes"],
                    "Exhaustiveness": job["exhaustiveness"],
                    "MinRMSD": job["min_rmsd"],
                    "EnergyRange": job["energy_range"],
                    "Seed": job["seed"],
                    "ConfigPath": rel(job["config_path"]),
                    "OutputPDBQT": rel(job["output_pdbqt"]),
                    "OutputLog": rel(job["output_log"]),
                    "LigandMetadataJSON": json.dumps(
                        job.get("ligand_metadata", {}),
                        sort_keys=True,
                        default=str,
                    ),
                }
            )

    return manifest_path


def write_provenance(
    path: str,
    *,
    start_stamp: str,
    end_stamp: str | None,
    status: str,
    vina_exe: str,
    vina_version: str,
    receptor_dir: str,
    receptor_file: str | None,
    ligand_dir: str,
    vina_csv: str,
    results_dir: str,
    num_modes: int,
    exhaustiveness: int,
    min_rmsd: float,
    energy_range: float,
    batch_seed: int,
    total_jobs: int,
    successes: int | None = None,
    failures: int | None = None,
) -> None:
    payload = {
        "schema_version": 1,
        "status": status,
        "start_time": start_stamp,
        "end_time": end_stamp,
        "vina": {
            "executable": vina_exe,
            "version": vina_version,
        },
        "inputs": {
            "receptor_dir": receptor_dir,
            "receptor_file": receptor_file,
            "ligand_dir": ligand_dir,
            "centers_csv": vina_csv,
        },
        "outputs": {
            "results_dir": results_dir,
            "job_manifest": "docking_jobs.csv",
        },
        "docking_parameters": {
            "num_modes": num_modes,
            "exhaustiveness": exhaustiveness,
            "min_rmsd_angstrom": min_rmsd,
            "energy_range_kcal_mol": energy_range,
            "pose_retention_mode": _pose_retention_mode(energy_range),
            "box_size_angstrom": [20.0, 20.0, 20.0],
            "threads_per_vina_job": 1,
        },
        "random_seed": {
            "batch_seed": batch_seed,
            "strategy": (
                "deterministic SHA-256-derived per-job seed from "
                "batch_seed + receptor + ligand + grid center"
            ),
        },
        "jobs": {
            "total": total_jobs,
            "successes": successes,
            "failures": failures,
        },
    }

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def cpu_count_cgroup_aware():
    """Respect cgroup/LSF limits; fall back to os.cpu_count()."""
    try:
        n = len(os.sched_getaffinity(0))
    except Exception:
        n = os.cpu_count() or 1

    # LSF often sets this; honor it if present.
    lsf_n = os.environ.get("LSB_DJOB_NUMPROC")
    if lsf_n:
        try:
            n = min(n, int(lsf_n))
        except ValueError:
            pass
    return max(1, n)


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Batch AutoDock Vina runner (scheduler-friendly). Provide flags for "
            "non-interactive use; falls back to interactive input otherwise."
        )
    )
    ap.add_argument("--receptors", help="Path to receptor folder")
    ap.add_argument(
        "--receptor-file",
        help=(
            "Run only this receptor file from --receptors. Used by 3B_ServerDocks.py "
            "to submit one independent LSF job per receptor."
        ),
    )
    ap.add_argument("--ligands", help="Path to ligand folder with .pdbqt ligands")
    ap.add_argument("--centers_csv", help="Path to vina centers CSV with headers PDB_ID,X,Y,Z")
    ap.add_argument("--poses", type=int, help="num_modes per ligand (e.g., 9, 20, 64)")
    ap.add_argument("--vina-exe", help="Path to AutoDock Vina executable")
    ap.add_argument(
        "--reserve_cores",
        type=int,
        default=1,
        help="How many cores to reserve for system/IO (default: 1)",
    )
    ap.add_argument(
        "--exhaustiveness",
        type=int,
        default=DEFAULT_EXHAUSTIVENESS,
        help=f"Vina global-search exhaustiveness (default: {DEFAULT_EXHAUSTIVENESS})",
    )
    ap.add_argument(
        "--min-rmsd",
        type=float,
        default=DEFAULT_MIN_RMSD,
        help=f"Minimum RMSD between retained output modes in Angstrom (default: {DEFAULT_MIN_RMSD})",
    )
    ap.add_argument(
        "--energy-range",
        type=float,
        default=None,
        help=(
            "Maximum kcal/mol above the best Vina mode to write. If omitted in "
            "non-interactive mode, Vina's standard 3.0 kcal/mol behavior is used."
        ),
    )
    ap.add_argument(
        "--full-pose-ensemble",
        action="store_true",
        help=(
            f"Use a wide {FULL_ENSEMBLE_ENERGY_RANGE:g} kcal/mol output window so "
            "downstream ensemble analyses can see essentially all retained modes."
        ),
    )
    ap.add_argument(
        "--seed",
        type=int,
        default=None,
        help=(
            "Reproducible batch seed. A deterministic unique seed is derived for "
            "each receptor-ligand job. Omit or use 0 for AUTO."
        ),
    )
    return ap.parse_args()


# ---------- Main ----------
def main():
    args = parse_args()

    if args.full_pose_ensemble and args.energy_range is not None:
        raise SystemExit("Use either --full-pose-ensemble OR --energy-range, not both.")

    start_time = time.time()
    start_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    timestamp_tag = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    cwd = os.getcwd()

    vina_exe = resolve_vina_executable(args.vina_exe)
    vina_version = get_vina_version(vina_exe)

    # Resolve inputs (flags preferred; else interactive).
    non_interactive = bool(args.receptors and args.ligands and args.centers_csv and args.poses)

    if non_interactive:
        receptor_dir = os.path.abspath(args.receptors)
        receptor_file = os.path.abspath(args.receptor_file) if args.receptor_file else None
        ligand_dir = os.path.abspath(args.ligands)
        vina_csv = os.path.abspath(args.centers_csv)
        num_modes = int(args.poses)
    else:
        dirs = [d for d in os.listdir(cwd) if os.path.isdir(d)]
        receptor_choice = choose("📁 Select receptor folder:", dirs)
        ligand_choice = choose("📁 Select ligand folder:", dirs)
        csv_files = [f for f in os.listdir(cwd) if f.endswith(".csv")]
        csv_choice = choose("📄 Select vina_centers.csv file:", csv_files)
        num_modes = int(input("\n🔢 How many poses per ligand? (e.g., 9, 20, 64): "))

        receptor_dir = os.path.join(cwd, receptor_choice)
        receptor_file = None
        ligand_dir = os.path.join(cwd, ligand_choice)
        vina_csv = os.path.join(cwd, csv_choice)

    exhaustiveness = int(args.exhaustiveness)
    min_rmsd = float(args.min_rmsd)

    if args.full_pose_ensemble:
        energy_range = FULL_ENSEMBLE_ENERGY_RANGE
    elif args.energy_range is not None:
        energy_range = float(args.energy_range)
    elif non_interactive:
        # Backward-compatible direct CLI behavior.
        energy_range = DEFAULT_ENERGY_RANGE
    else:
        retention_choice = input(
            "\nPose retention for out.pdbqt?\n"
            f" [1] Standard Vina ({DEFAULT_ENERGY_RANGE:g} kcal/mol)\n"
            f" [2] Full pose ensemble / PROTACability ({FULL_ENSEMBLE_ENERGY_RANGE:g} kcal/mol)\n"
            "Choose 1/2 [1]: "
        ).strip() or "1"
        if retention_choice == "2":
            energy_range = FULL_ENSEMBLE_ENERGY_RANGE
        else:
            energy_range = DEFAULT_ENERGY_RANGE

        exhaustiveness = int(input_default("Vina exhaustiveness", exhaustiveness))
        min_rmsd = float(input_default("Minimum RMSD between retained modes (Angstrom)", min_rmsd))

    _validate_vina_parameters(num_modes, exhaustiveness, min_rmsd, energy_range)
    batch_seed = _normalize_batch_seed(args.seed)

    # --- Read grid centers FIRST (so we can validate before building names/dirs) ---
    with open(vina_csv, "r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        grid_rows = list(reader)
        fieldnames = reader.fieldnames or []

    required_cols = {"PDB_ID", "X", "Y", "Z"}
    if not required_cols.issubset(set(fieldnames)):
        raise ValueError(f"CSV {vina_csv} must have headers: {sorted(required_cols)}")

    # --- Build naming tags from folder names only (no centers tag) ---
    rec_tag = _receptor_tag_from_dir(receptor_dir)
    if receptor_file:
        rec_tag = f"{rec_tag}_{_sanitize_tag(_strip_known_suffix(Path(receptor_file).name))}"
    lig_tag = _ligand_tag_from_dir(ligand_dir)

    # --- Outputs with receptor/ligand/poses/timestamp ---
    results_dir = os.path.join(
        cwd,
        f"Docking_Results_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}",
    )
    os.makedirs(results_dir, exist_ok=True)

    # Run log (append all results; keeps terminal quiet) — tagged.
    run_log_path = os.path.join(
        cwd,
        f"run_log_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt",
    )
    run_log = open(run_log_path, "a", encoding="utf-8")

    jobs = build_jobs(
        results_dir,
        ligand_dir,
        receptor_dir,
        grid_rows,
        num_modes,
        vina_exe,
        exhaustiveness,
        min_rmsd,
        energy_range,
        batch_seed,
        receptor_file,
    )

    receptor_files = []
    for suffix in RECEPTOR_SUFFIXES:
        receptor_files.extend(glob.glob(os.path.join(receptor_dir, f"*{suffix}")))
    if receptor_file:
        receptor_files = [
            path for path in receptor_files
            if os.path.abspath(path) == receptor_file
        ]

    missing_receptor_msgs = []
    for row in grid_rows:
        resolved_receptor, match_kind = _resolve_receptor_file(row["PDB_ID"], receptor_files)
        if not resolved_receptor:
            missing_receptor_msgs.append(
                f"❌ Unmatched receptor center: {row['PDB_ID']} (match={match_kind})"
            )

    for msg in missing_receptor_msgs:
        run_log.write(msg + "\n")
    run_log.flush()

    ncpus = cpu_count_cgroup_aware()
    reserve = max(0, int(args.reserve_cores))
    max_workers = max(1, ncpus - reserve)

    total_jobs = len(jobs)
    if total_jobs == 0:
        run_log.close()
        available = ", ".join(sorted(os.path.basename(path) for path in receptor_files)[:10]) or "(none found)"
        details = "; ".join(missing_receptor_msgs[:5]) if missing_receptor_msgs else "no receptor-center rows matched"
        raise RuntimeError(
            "No jobs prepared. Check receptors/ligands/CSV inputs. "
            f"Available receptors: {available}. Details: {details}"
        )

    # Machine-readable provenance for future 8_Protacability.py.
    job_manifest_path = write_job_manifest(results_dir, jobs)
    provenance_path = os.path.join(results_dir, "docking_parameters.json")
    write_provenance(
        provenance_path,
        start_stamp=start_stamp,
        end_stamp=None,
        status="prepared",
        vina_exe=vina_exe,
        vina_version=vina_version,
        receptor_dir=receptor_dir,
        receptor_file=receptor_file,
        ligand_dir=ligand_dir,
        vina_csv=vina_csv,
        results_dir=results_dir,
        num_modes=num_modes,
        exhaustiveness=exhaustiveness,
        min_rmsd=min_rmsd,
        energy_range=energy_range,
        batch_seed=batch_seed,
        total_jobs=total_jobs,
    )

    print(
        f"\n🧠 Detected {ncpus} schedulable cores. "
        f"Running up to {max_workers} concurrent Vina jobs (reserve {reserve}).",
        flush=True,
    )
    print(f"🗂  Results dir      : {results_dir}", flush=True)
    print(f"🧪 Vina executable  : {vina_exe}", flush=True)
    print(f"🧾 Vina version     : {vina_version}", flush=True)
    print(f"🎯 Requested modes  : {num_modes}", flush=True)
    print(f"🔎 Exhaustiveness   : {exhaustiveness}", flush=True)
    print(f"📐 Minimum RMSD     : {min_rmsd:g} Å", flush=True)
    print(f"⚡ Energy range     : {energy_range:g} kcal/mol ({_pose_retention_mode(energy_range)})", flush=True)
    print(f"🎲 Batch seed       : {batch_seed}", flush=True)
    print(f"📋 Job manifest     : {job_manifest_path}", flush=True)
    print(f"🧬 Dock provenance  : {provenance_path}", flush=True)

    successes, failures = 0, 0
    results_log = []

    progress_bar(0, total_jobs, successes=0, failures=0)

    # Concurrency limited by max_workers; quiet terminal, log to file.
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

    # Update provenance after completion.
    write_provenance(
        provenance_path,
        start_stamp=start_stamp,
        end_stamp=end_stamp,
        status="complete" if failures == 0 else "complete_with_failures",
        vina_exe=vina_exe,
        vina_version=vina_version,
        receptor_dir=receptor_dir,
        receptor_file=receptor_file,
        ligand_dir=ligand_dir,
        vina_csv=vina_csv,
        results_dir=results_dir,
        num_modes=num_modes,
        exhaustiveness=exhaustiveness,
        min_rmsd=min_rmsd,
        energy_range=energy_range,
        batch_seed=batch_seed,
        total_jobs=total_jobs,
        successes=successes,
        failures=failures,
    )

    # Duration file in working directory (tagged).
    duration_txt = os.path.join(
        cwd,
        f"job_duration_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt",
    )
    with open(duration_txt, "w", encoding="utf-8") as df:
        df.write(f"Start Time : {start_stamp}\n")
        df.write(f"End Time   : {end_stamp}\n")
        df.write(f"Duration   : {duration_min} minutes\n")
        df.write(f"Duration_s : {duration_sec} seconds\n")
        df.write(f"Jobs Total : {total_jobs}\n")
        df.write(f"Successes  : {successes}\n")
        df.write(f"Failures   : {failures}\n")
        df.write(f"Batch Seed : {batch_seed}\n")

    # Summary file (tagged).
    summary_path = os.path.join(
        cwd,
        f"docking_summary_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt",
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
        f.write(" - Threads per Vina job: 1 (via --cpu 1 / OMP_NUM_THREADS=1)\n\n")

        f.write("Vina Parameters:\n")
        f.write(f" - Vina executable : {vina_exe}\n")
        f.write(f" - Vina version    : {vina_version}\n")
        f.write(f" - num_modes       : {num_modes}\n")
        f.write(f" - exhaustiveness  : {exhaustiveness}\n")
        f.write(f" - min_rmsd        : {min_rmsd} Angstrom\n")
        f.write(f" - energy_range    : {energy_range} kcal/mol\n")
        f.write(f" - retention_mode  : {_pose_retention_mode(energy_range)}\n")
        f.write(f" - batch_seed      : {batch_seed}\n")
        f.write(" - per-job seeds   : deterministic SHA-256 derivation from batch seed + job identity\n\n")

        f.write("Inputs/Outputs:\n")
        f.write(f" - Receptor dir    : {receptor_dir}\n")
        f.write(f" - Ligand dir      : {ligand_dir}\n")
        f.write(f" - Centers CSV     : {vina_csv}\n")
        f.write(f" - Results dir     : {results_dir}\n")
        f.write(f" - Job manifest    : {job_manifest_path}\n")
        f.write(f" - Parameters JSON : {provenance_path}\n\n")

        f.write("Results:\n--------\n")
        for line in results_log:
            f.write(line + "\n")

    print(f"\n⏱️  Duration file : {duration_txt}", flush=True)
    print(f"🧾 Run log       : {run_log_path}", flush=True)
    print(f"📄 Summary saved : {summary_path}", flush=True)
    print(f"🧬 Parameters    : {provenance_path}", flush=True)
    print(f"📋 Job manifest : {job_manifest_path}", flush=True)
    print(f"✅ Done. Results in: {results_dir}", flush=True)


if __name__ == "__main__":
    main()
