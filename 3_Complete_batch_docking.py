#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
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
import tempfile
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
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
VINA_RESULT_RE = re.compile(
    r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)\s+(\S+)\s+(\S+)", re.MULTILINE
)


def _write_vina_config(job: dict, config_path: str, output_pdbqt: str) -> None:
    """Write the unchanged Vina configuration, optionally to compact scratch."""
    with open(config_path, "w", encoding="utf-8") as cfg:
        cfg.write(f"receptor = {job['receptor_file']}\n")
        cfg.write(f"ligand = {job['ligand_file']}\n")
        cfg.write(f"center_x = {job['cx']}\n")
        cfg.write(f"center_y = {job['cy']}\n")
        cfg.write(f"center_z = {job['cz']}\n")
        cfg.write("size_x = 20\nsize_y = 20\nsize_z = 20\n")
        cfg.write(f"num_modes = {job['num_modes']}\nexhaustiveness = {job['exhaustiveness']}\n")
        cfg.write(f"min_rmsd = {job['min_rmsd']}\nenergy_range = {job['energy_range']}\n")
        cfg.write(f"seed = {job['seed']}\nout = {output_pdbqt}\n")


def _compact_rows(job: dict, output_pdbqt: str) -> tuple[list[dict], bytes]:
    """Extract scores before the worker's private scratch directory is removed."""
    pose_bytes = Path(output_pdbqt).read_bytes()
    text = pose_bytes.decode("utf-8", errors="replace")
    rows = []
    for pose, match in enumerate(VINA_RESULT_RE.finditer(text), start=1):
        affinity, rmsd_lb, rmsd_ub = match.groups()
        rows.append({
            "PDB_ID": job["pdbid"], "ResolvedReceptor": job["resolved_receptor"],
            "Ligand": job["ligand"], "Pose": pose, "Binding_Affinity": affinity,
            "RMSD_LB": rmsd_lb, "RMSD_UB": rmsd_ub, "Seed": job["seed"],
            "LigandMetadataJSON": json.dumps(job.get("ligand_metadata", {}), sort_keys=True, default=str),
            "Status": "success",
        })
    return rows, pose_bytes


def run_docking(job):
    """Run one ligand safely; no single worker error may abort a whole batch."""
    try:
        vina_exe = job["vina_exe"]

        if job.get("compact_scratch_root"):
            # Parent removes this private directory only after its CSV record is
            # fsync'd, so an interruption cannot mark a cleaned result complete.
            work = tempfile.mkdtemp(prefix="vina_job_", dir=job["compact_scratch_root"])
            try:
                config_path = os.path.join(work, "config.txt")
                output_pdbqt = os.path.join(work, "out.pdbqt")
                _write_vina_config(job, config_path, output_pdbqt)
                ok, line, rows, pose_bytes = _run_vina(job, config_path, output_pdbqt, compact=True)
                return ok, line, rows, work, pose_bytes
            except Exception as exc:
                return False, f"❌ {job['ligand']} vs {job['pdbid']}: {exc}", [], work, b""
        _write_vina_config(job, job["config_path"], job["output_pdbqt"])
        return _run_vina(job, job["config_path"], job["output_pdbqt"], compact=False)
    except Exception as exc:
        message = f"runner exception: {type(exc).__name__}: {exc}"
        if job.get("compact_scratch_root"):
            return False, f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {message}", [], "", b""
        return False, f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {message}"


def _run_vina(job: dict, config_path: str, output_pdbqt: str, *, compact: bool):
    """Run Vina once; compact callers receive score rows rather than pose paths."""
    try:
        vina_exe = job["vina_exe"]

        # Keep each Vina process to 1 thread; the Python pool controls total concurrency.
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        result = subprocess.run(
            [vina_exe, "--config", config_path, "--cpu", "1"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
        )

        if result.returncode == 0:
            line = f"✅ {job['ligand']} → {job['pdbid']}  seed={job['seed']}"
            if compact:
                rows, pose_bytes = _compact_rows(job, output_pdbqt)
                if not rows:
                    return False, f"❌ {job['ligand']} vs {job['pdbid']}: Vina wrote no score rows", [], b""
                return True, line, rows, pose_bytes
            with open(job["output_log"], "wb") as log_file:
                log_file.write(result.stdout)
                if result.stderr:
                    log_file.write(b"\n--- STDERR ---\n" + result.stderr)
            return True, line

        err = result.stderr.decode(errors="ignore").strip()
        line = f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {err}"
        if compact:
            return False, line, [], b""
        with open(job["output_log"], "wb") as log_file:
            log_file.write(result.stdout)
            if result.stderr:
                log_file.write(b"\n--- STDERR ---\n" + result.stderr)
        return False, line
    except Exception as exc:
        # In particular, fork_exec can raise ENOMEM. Record it for retry while
        # allowing all other queued ligands to finish.
        message = f"runner exception: {type(exc).__name__}: {exc}"
        if compact:
            return False, f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {message}", [], b""
        try:
            with open(job["output_log"], "ab") as log_file:
                log_file.write((message + "\n").encode("utf-8", errors="replace"))
        except OSError:
            pass
        return False, f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: {message}"


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
    exhaustiveness=DEFAULT_EXHAUSTIVENESS,
    min_rmsd=DEFAULT_MIN_RMSD,
    energy_range=DEFAULT_ENERGY_RANGE,
    batch_seed=1,
    receptor_file=None,
    ligand_files=None,
    compact_scratch_root: str | None = None,
):
    jobs = []
    ligand_files = list(ligand_files) if ligand_files is not None else glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
    ligand_files = [os.path.abspath(path) for path in ligand_files]
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

    ligand_by_name = {os.path.splitext(os.path.basename(path))[0]: path for path in ligand_files}
    if len(ligand_by_name) != len(ligand_files):
        raise ValueError("Ligand file names must be unique when using --ligand-list")
    ligand_names = sorted(ligand_by_name)

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
            ligand_file = ligand_by_name[ligand]

            output_subdir = os.path.join(results_dir, safe_pdbid, ligand)
            # Legacy mode deliberately retains one directory per ligand. Compact
            # mode never creates this tree: each worker uses private scratch.
            if not compact_scratch_root:
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
                    "compact_scratch_root": compact_scratch_root,
                    "ligand_metadata": merge_ligand_metadata(
                        ligand,
                        manifest_row=manifest_map.get(ligand),
                    ),
                }
            )

    return jobs


def write_job_manifest(results_dir: str, jobs: list[dict], *, compact_mode: bool = False) -> str:
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
                    "ConfigPath": "temporary scratch" if compact_mode else rel(job["config_path"]),
                    "OutputPDBQT": "not retained (compact results)" if compact_mode else rel(job["output_pdbqt"]),
                    "OutputLog": "compact_failures.csv" if compact_mode else rel(job["output_log"]),
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
    compact: dict | None = None,
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
    if compact is not None:
        payload["compact_results"] = compact

    _atomic_json(Path(path), payload)


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


COMPACT_RESULT_COLUMNS = [
    "PDB_ID", "ResolvedReceptor", "Ligand", "Pose", "Binding_Affinity", "RMSD_LB",
    "RMSD_UB", "Seed", "LigandMetadataJSON", "Status",
]
COMPACT_FAILURE_COLUMNS = ["PDB_ID", "ResolvedReceptor", "Ligand", "Seed", "Error", "Status"]


def resolve_compact_scratch_root(cli_root: str | None, results_dir: str) -> tuple[str, str]:
    """Apply CLI > environment/profile export > results-root fallback precedence."""
    if cli_root:
        return os.path.abspath(cli_root), "cli"
    env_root = os.environ.get("VINA_SCRATCH_ROOT", "").strip()
    if env_root:
        return os.path.abspath(os.path.expandvars(os.path.expanduser(env_root))), "environment/profile"
    return os.path.abspath(results_dir), "results-root fallback (temporary; not verified local NVMe)"


def stage_compact_inputs(
    receptor_dir: str,
    receptor_file: str | None,
    ligand_files: list[str],
    centers_csv: str,
    scratch_root: str,
) -> tuple[str, str | None, list[str], str, str]:
    """Stage one bounded ligand shard and its small shared inputs to local SSD."""
    workspace = tempfile.mkdtemp(prefix="vina_inputs_", dir=scratch_root)
    staged_receptors = Path(workspace) / "Receptors"
    staged_ligands = Path(workspace) / "Ligands"
    staged_receptors.mkdir()
    staged_ligands.mkdir()

    for path in Path(receptor_dir).iterdir():
        if path.is_file() and path.suffix.lower() in RECEPTOR_SUFFIXES:
            shutil.copy2(path, staged_receptors / path.name)

    staged_ligand_files: list[str] = []
    copied_names: set[str] = set()
    for source in ligand_files:
        src = Path(source)
        if src.name in copied_names:
            raise ValueError(f"Duplicate ligand filename in shard: {src.name}")
        if not src.is_file():
            raise FileNotFoundError(f"Ligand listed for staging is missing: {src}")
        target = staged_ligands / src.name
        shutil.copy2(src, target)
        copied_names.add(src.name)
        staged_ligand_files.append(str(target))

    staged_centers = Path(workspace) / Path(centers_csv).name
    shutil.copy2(centers_csv, staged_centers)
    staged_receptor_file = str(staged_receptors / Path(receptor_file).name) if receptor_file else None
    if staged_receptor_file and not Path(staged_receptor_file).is_file():
        raise FileNotFoundError(f"Selected receptor listed for staging is missing: {receptor_file}")
    return str(staged_receptors), staged_receptor_file, staged_ligand_files, str(staged_centers), workspace


def _job_key(job: dict) -> tuple[str, str, str]:
    return str(job["resolved_receptor"]), str(job["ligand"]), str(job["seed"])


def _completed_compact_jobs(path: Path) -> set[tuple[str, str, str]]:
    if not path.is_file():
        return set()
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row.get("ResolvedReceptor", ""), row.get("Ligand", ""), row.get("Seed", ""))
            for row in csv.DictReader(handle)
            if row.get("Status") == "success"
        }


def _append_csv_rows(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    """Durably append final records; a record is complete only after fsync."""
    if not rows:
        return
    new_file = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if new_file:
            writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())


def _atomic_json(path: Path, payload: dict) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)


def _write_retained_poses(results_dir: Path, candidates: dict[str, list[tuple[float, dict, bytes]]]) -> list[dict]:
    """Atomically retain only ranked compact pose files, never a full job tree."""
    retained: list[dict] = []
    root = results_dir / "retained_poses"
    for receptor, entries in sorted(candidates.items()):
        receptor_dir = root / _sanitize_tag(receptor)
        receptor_dir.mkdir(parents=True, exist_ok=True)
        for rank, (score, job, pose_bytes) in enumerate(sorted(entries, key=lambda item: (item[0], item[1]["ligand"], item[1]["seed"])), start=1):
            name = f"top_{rank:03d}_{_sanitize_tag(job['ligand'])}_seed{job['seed']}.pdbqt"
            path = receptor_dir / name
            tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
            tmp.write_bytes(pose_bytes)
            os.replace(tmp, path)
            retained.append({
                "rank": rank, "receptor": receptor, "ligand": job["ligand"], "seed": job["seed"],
                "binding_affinity": score, "path": str(path.relative_to(results_dir)),
            })
    _atomic_json(results_dir / "compact_retained_poses.json", {"schema_version": 1, "retained": retained})
    return retained


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
        "--ligand-list",
        help="Optional text file containing one ligand .pdbqt path per line; enables sharded/resume runs.",
    )
    ap.add_argument(
        "--results-dir",
        help="Optional explicit result directory. Required by array wrappers to isolate each shard.",
    )
    ap.add_argument(
        "--compact-results",
        action="store_true",
        help="Collect scores into shard CSVs and discard per-ligand pose/config/log files after extraction.",
    )
    ap.add_argument(
        "--compact-scratch-root",
        help="Temporary Vina workspace. Overrides VINA_SCRATCH_ROOT; defaults to --results-dir temporarily.",
    )
    ap.add_argument(
        "--stage-inputs-to-scratch",
        action="store_true",
        help="Copy this bounded ligand shard, receptors, and centers CSV to compact scratch before docking.",
    )
    ap.add_argument(
        "--compact-retain-top-poses", type=int, default=0, metavar="N",
        help="In compact mode, retain full PDBQT poses only for the top N ligands per receptor (default: 0).",
    )
    ap.add_argument(
        "--workers",
        type=int,
        help="Maximum concurrent Vina processes. Defaults to schedulable cores minus --reserve_cores.",
    )
    ap.add_argument(
        "--max-pending",
        type=int,
        help="Maximum submitted-but-unfinished worker tasks; defaults to twice --workers.",
    )
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
    if args.compact_retain_top_poses < 0:
        raise SystemExit("--compact-retain-top-poses must be zero or greater")
    if args.compact_retain_top_poses and not args.compact_results:
        raise SystemExit("--compact-retain-top-poses requires --compact-results")
    if args.stage_inputs_to_scratch and not args.compact_results:
        raise SystemExit("--stage-inputs-to-scratch requires --compact-results")

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
    results_dir = os.path.abspath(args.results_dir) if args.results_dir else os.path.join(
        cwd, f"Docking_Results_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}"
    )
    os.makedirs(results_dir, exist_ok=True)
    compact_scratch_root = None
    scratch_source = None
    if args.compact_results:
        compact_scratch_root, scratch_source = resolve_compact_scratch_root(args.compact_scratch_root, results_dir)
        os.makedirs(compact_scratch_root, exist_ok=True)

    # Run log (append all results; keeps terminal quiet) — tagged.
    run_log_path = os.path.join(
        cwd,
        f"run_log_{rec_tag}_{lig_tag}_{num_modes}Poses_{timestamp_tag}.txt",
    )
    run_log = open(run_log_path, "a", encoding="utf-8")

    ligand_files = None
    if args.ligand_list:
        list_path = Path(args.ligand_list).resolve()
        if not list_path.is_file():
            raise FileNotFoundError(f"Ligand list does not exist: {list_path}")
        ligand_files = []
        for raw_line in list_path.read_text(encoding="utf-8").splitlines():
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            path = Path(value)
            path = path if path.is_absolute() else list_path.parent / path
            if path.suffix.lower() != ".pdbqt" or not path.is_file():
                raise FileNotFoundError(f"Invalid ligand path in {list_path}: {value}")
            ligand_files.append(str(path.resolve()))

    staged_input_workspace = None
    if args.stage_inputs_to_scratch:
        ligand_files_for_staging = ligand_files or glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
        receptor_dir, receptor_file, ligand_files, vina_csv, staged_input_workspace = stage_compact_inputs(
            receptor_dir, receptor_file, ligand_files_for_staging, vina_csv, compact_scratch_root
        )
        atexit.register(shutil.rmtree, staged_input_workspace, ignore_errors=True)

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
        ligand_files,
        compact_scratch_root,
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
    max_workers = int(args.workers) if args.workers is not None else max(1, ncpus - reserve)
    if max_workers < 1:
        raise ValueError("--workers must be greater than zero")
    if max_workers > ncpus:
        raise ValueError(f"--workers ({max_workers}) exceeds schedulable cores ({ncpus})")
    max_pending = int(args.max_pending) if args.max_pending is not None else max_workers * 2
    if max_pending < max_workers:
        raise ValueError("--max-pending must be at least --workers")

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
    job_manifest_path = write_job_manifest(results_dir, jobs, compact_mode=args.compact_results)
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
        compact=({
            "enabled": True, "results_csv": "compact_results.csv", "failures_csv": "compact_failures.csv",
            "completion_metadata": "compact_completion.json", "scratch_root": compact_scratch_root,
            "scratch_source": scratch_source, "poses_retained": bool(args.compact_retain_top_poses),
            "retain_top_poses_per_receptor": args.compact_retain_top_poses,
        } if args.compact_results else None),
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
    if args.compact_results:
        pose_note = "no pose PDBQT files retained" if not args.compact_retain_top_poses else f"top {args.compact_retain_top_poses} pose PDBQT files retained per receptor"
        print(f"🗜  Compact results : enabled ({pose_note})", flush=True)
        print(f"🧹 Scratch root     : {compact_scratch_root} ({scratch_source})", flush=True)

    successes, failures = 0, 0
    results_log = []
    retained_candidates: dict[str, list[tuple[float, dict, bytes]]] = {}

    compact_results_path = Path(results_dir) / "compact_results.csv"
    compact_failures_path = Path(results_dir) / "compact_failures.csv"
    compact_completion_path = Path(results_dir) / "compact_completion.json"
    retained_manifest_path = Path(results_dir) / "compact_retained_poses.json"
    completed_jobs = _completed_compact_jobs(compact_results_path) if args.compact_results else set()
    replay_completed_for_retention = bool(
        args.compact_results and args.compact_retain_top_poses and completed_jobs and not retained_manifest_path.exists()
    )
    if completed_jobs and not replay_completed_for_retention:
        jobs = [job for job in jobs if _job_key(job) not in completed_jobs]
        print(f"♻️  Compact resume: skipping {len(completed_jobs)} completed docking jobs.", flush=True)
    elif replay_completed_for_retention:
        print("♻️  Compact resume: replaying completed dockings to reconstruct requested retained poses.", flush=True)

    progress_bar(0, len(jobs), successes=0, failures=0)

    # Keep only a bounded task buffer. Submitting tens of thousands of futures
    # consumes unnecessary parent-process memory and makes an OOM more likely.
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        job_iter = iter(jobs)
        futures = {}
        for _ in range(min(max_pending, len(jobs))):
            job = next(job_iter)
            futures[executor.submit(run_docking, job)] = job
        done_count = 0
        while futures:
            completed, _ = wait(futures, return_when=FIRST_COMPLETED)
            for future in completed:
                job = futures.pop(future)
                try:
                    outcome = future.result()
                    if args.compact_results:
                        ok, line, score_rows, scratch_job_dir, pose_bytes = outcome
                    else:
                        ok, line = outcome
                        score_rows = []
                        scratch_job_dir = ""
                        pose_bytes = b""
                except Exception as exc:
                    ok = False
                    score_rows = []
                    scratch_job_dir = ""
                    pose_bytes = b""
                    line = (
                        f"❌ {job['ligand']} vs {job['pdbid']}  seed={job['seed']}: "
                        f"worker future exception: {type(exc).__name__}: {exc}"
                    )
                results_log.append(line)
                run_log.write(line + "\n")
                run_log.flush()

                if ok:
                    successes += 1
                    if args.compact_results and _job_key(job) not in completed_jobs:
                        _append_csv_rows(compact_results_path, COMPACT_RESULT_COLUMNS, score_rows)
                    if args.compact_retain_top_poses and score_rows:
                        candidates = retained_candidates.setdefault(job["resolved_receptor"], [])
                        candidates.append((min(float(row["Binding_Affinity"]) for row in score_rows), job, pose_bytes))
                        candidates.sort(key=lambda item: (item[0], item[1]["ligand"], item[1]["seed"]))
                        del candidates[args.compact_retain_top_poses:]
                else:
                    failures += 1
                    if args.compact_results:
                        _append_csv_rows(compact_failures_path, COMPACT_FAILURE_COLUMNS, [{
                            "PDB_ID": job["pdbid"], "ResolvedReceptor": job["resolved_receptor"],
                            "Ligand": job["ligand"], "Seed": job["seed"], "Error": line, "Status": "failed",
                        }])
                if scratch_job_dir:
                    # The directory was made by this run with mkdtemp; never
                    # remove a configured scratch root or existing results dir.
                    shutil.rmtree(scratch_job_dir, ignore_errors=True)

                done_count += 1
                progress_bar(done_count, len(jobs), successes=successes, failures=failures)

                try:
                    next_job = next(job_iter)
                except StopIteration:
                    continue
                futures[executor.submit(run_docking, next_job)] = next_job

    end_time = time.time()
    duration_min = round((end_time - start_time) / 60, 2)
    duration_sec = round(end_time - start_time, 2)
    end_stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_log.close()

    if args.compact_results:
        retained = _write_retained_poses(Path(results_dir), retained_candidates) if args.compact_retain_top_poses else []
        _atomic_json(compact_completion_path, {
            "schema_version": 1,
            "status": "complete" if failures == 0 else "complete_with_failures",
            "completed_at": end_stamp,
            "results_csv": compact_results_path.name,
            "failures_csv": compact_failures_path.name,
            "successful_jobs_this_run": successes,
            "failed_jobs_this_run": failures,
            "resumed_completed_jobs": len(completed_jobs),
            "scratch_root": compact_scratch_root,
            "scratch_source": scratch_source,
            "poses_retained": bool(args.compact_retain_top_poses),
            "retained_top_poses_per_receptor": args.compact_retain_top_poses,
            "retained_pose_files": retained,
        })

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
        compact=({
            "enabled": True, "results_csv": "compact_results.csv", "failures_csv": "compact_failures.csv",
            "completion_metadata": "compact_completion.json", "scratch_root": compact_scratch_root,
            "scratch_source": scratch_source, "poses_retained": bool(args.compact_retain_top_poses),
            "retain_top_poses_per_receptor": args.compact_retain_top_poses,
        } if args.compact_results else None),
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
