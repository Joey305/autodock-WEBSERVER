"""Non-interactive campaign assembly for AutoDock/Vina package uploads."""
from __future__ import annotations

import csv
import io
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from typing import Any


LIGAND_SUFFIXES = {".csv", ".sdf", ".smi", ".smiles"}
RUNTIME_FILES = ("0_LIGSPLIT.py", "1_ConformerGeneration.py", "3_Complete_batch_docking.py", "ligand_manifest.py", "ligand_naming.py", "campaign_status.py")


def safe_name(value: str, fallback: str = "campaign-001") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", (value or "").strip()).strip(".-_")
    return cleaned[:80] or fallback


def _safe_member(name: str) -> Path | None:
    rel = Path((name or "").replace("\\", "/").strip("/"))
    if not rel.parts or rel.is_absolute() or ".." in rel.parts:
        return None
    if any(part in {"__MACOSX", ".DS_Store"} or part.startswith(".") for part in rel.parts):
        return None
    return rel


def extract_package(upload, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    # Werkzeug's upload stream is not guaranteed to expose ``seekable`` on
    # every supported Python version; ZipFile needs a fully seekable source.
    with zipfile.ZipFile(io.BytesIO(upload.read())) as archive:
        for info in archive.infolist():
            if info.is_dir() or (info.external_attr >> 16) & 0o170000 == 0o120000:
                continue
            rel = _safe_member(info.filename)
            if rel is None:
                raise ValueError(f"Unsafe ZIP entry: {info.filename}")
            out = (destination / rel).resolve()
            if destination.resolve() not in out.parents:
                raise ValueError(f"Unsafe ZIP entry: {info.filename}")
            out.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as source, out.open("wb") as target:
                shutil.copyfileobj(source, target)


def _file_list(root: Path, suffixes: set[str]) -> list[Path]:
    return sorted((p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in suffixes and not any(part.startswith(".") for part in p.relative_to(root).parts)), key=lambda p: p.as_posix().lower())


def _count_input(path: Path) -> int:
    if path.suffix.lower() == ".csv":
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            next(reader, None)
            return sum(1 for row in reader if any(cell.strip() for cell in row))
    if path.suffix.lower() == ".sdf":
        return sum(1 for line in path.open(encoding="utf-8", errors="replace") if line.strip() == "$$$$")
    return sum(1 for line in path.open(encoding="utf-8", errors="replace") if line.strip())


def _count_center_rows(path: Path) -> int:
    with path.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        return sum(1 for row in csv.DictReader(handle) if any((value or "").strip() for value in row.values()))


def inspect_package(root: Path) -> dict[str, Any]:
    ligands = _file_list(root, LIGAND_SUFFIXES)
    # Ignore result artefacts; users should choose only original ligand inputs.
    ligands = [
        p for p in ligands
        if "Ligands_PDBQT" not in p.as_posix()
        and "Docking_Results" not in p.as_posix()
        and not p.name.lower().startswith("vina_centers")
    ]
    receptors = _file_list(root, {".pdbqt"})
    centers = sorted(root.rglob("vina_centers*.csv"), key=lambda p: p.as_posix().lower())
    candidates = []
    for item in ligands:
        rel = item.relative_to(root).as_posix()
        entry: dict[str, Any] = {"path": rel, "type": item.suffix.lower().lstrip("."), "ligand_count": _count_input(item)}
        if item.suffix.lower() == ".csv":
            with item.open(newline="", encoding="utf-8-sig") as handle:
                entry["headers"] = next(csv.reader(handle), [])
        candidates.append(entry)
    return {
        "ligand_inputs": candidates,
        "receptors": [p.relative_to(root).as_posix() for p in receptors],
        "centers": [p.relative_to(root).as_posix() for p in centers],
        "required_scripts_present": [name for name in RUNTIME_FILES if (root / name).is_file()],
    }


def _quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _lsf_header(
    job: str, queue: str, project: str, workers: int, memory: int, walltime: str,
    *, exclusive_host: bool = False, memory_limit_mb: int = 0,
) -> str:
    lines = [
        "#!/bin/bash", f"#BSUB -J {job}", f"#BSUB -q {queue}", f"#BSUB -P {project}",
        f"#BSUB -n {workers}", f'#BSUB -R "rusage[mem={memory}]"', f"#BSUB -W {walltime}",
    ]
    if exclusive_host:
        lines.append("#BSUB -x")
    if memory_limit_mb:
        lines.append(f"#BSUB -M {memory_limit_mb}")
    lines.extend((f"#BSUB -oo logs/{job}.%J.out", f"#BSUB -eo logs/{job}.%J.err", "set -euo pipefail", ""))
    return "\n".join(lines)


def build_campaign(
    *, source_root: Path, app_root: Path, output_root: Path, campaign_name: str, ligand_rel: str,
    split_mode: str, split_value: int, batch_prefix: str, settings: dict[str, Any],
) -> dict[str, Any]:
    campaign_name = safe_name(campaign_name)
    ligand = (source_root / ligand_rel).resolve()
    if source_root.resolve() not in ligand.parents or not ligand.is_file() or ligand.suffix.lower() not in LIGAND_SUFFIXES:
        raise ValueError("Choose a valid ligand input found in the uploaded ZIP.")
    if split_mode not in {"group_size", "num_groups"} or split_value < 1:
        raise ValueError("Choose a positive batch size or number of batches.")
    batch_prefix = safe_name(batch_prefix, "batch")
    campaign = output_root / campaign_name
    if campaign.exists():
        shutil.rmtree(campaign)
    shutil.copytree(source_root, campaign)
    for filename in RUNTIME_FILES:
        source = app_root / filename
        target = campaign / filename
        if not target.exists() and source.is_file():
            shutil.copy2(source, target)
    missing = [name for name in ("0_LIGSPLIT.py", "1_ConformerGeneration.py", "3_Complete_batch_docking.py") if not (campaign / name).is_file()]
    if missing:
        raise ValueError("Package is missing required workflow scripts: " + ", ".join(missing))

    rel_in_campaign = ligand.relative_to(source_root.resolve())
    batches = campaign / "batches"
    command = [sys.executable, str(campaign / "0_LIGSPLIT.py"), "--input", str(campaign / rel_in_campaign), "--output", str(batches), f"--{split_mode.replace('_', '-')}", str(split_value), "--prefix", batch_prefix, "--yes"]
    subprocess.run(command, cwd=campaign, check=True, capture_output=True, text=True)
    shards = sorted(p for p in batches.glob(f"{batch_prefix}_*") if p.is_dir())
    if not shards:
        raise ValueError("The ligand splitter did not create any batches.")
    receptors = [p for p in campaign.rglob("*.pdbqt") if "Ligands_PDBQT" not in p.as_posix()]
    centers = sorted(campaign.rglob("vina_centers*.csv"))
    if not receptors or not centers:
        raise ValueError("The ZIP must include prepared receptor .pdbqt file(s) and vina_centers.csv.")

    queue, project = settings.get("queue", "normal"), settings.get("project", "brd")
    conf_workers, vina_workers = int(settings.get("confgen_workers", 16)), int(settings.get("vina_workers", 8))
    memory = int(settings.get("mem_per_core", 8000))
    vina_exclusive_host = bool(settings.get("vina_exclusive_host", False))
    vina_memory_limit_mb = int(settings.get("vina_memory_limit_mb", 0))
    vina_max_rss_mb = int(settings.get("vina_max_rss_mb", 0))
    vina_timeout_seconds = int(settings.get("vina_timeout_seconds", 0))
    conf_poses, vina_poses = int(settings.get("confgen_poses", 64)), int(settings.get("vina_poses", 20))
    common_setup = settings.get("setup", "")
    state_flags = " ".join(settings.get("state_flags", []))
    csv_smiles, csv_id = settings.get("csv_smiles_col", ""), settings.get("csv_id_col", "")
    center_rel = centers[0].relative_to(campaign).as_posix()
    receptor_dir = receptors[0].parent.relative_to(campaign).as_posix()
    batch_progress = []
    state_multiplier = int(settings.get("max_protomers", 1)) * int(settings.get("max_tautomers", 1)) if settings.get("enumerate_protomers") else 1
    docking_targets = _count_center_rows(centers[0])
    for index, shard in enumerate(shards, start=1):
        shard_rel = shard.relative_to(campaign).as_posix()
        job_suffix = f"{index:03d}"
        input_files = sorted(p for p in shard.iterdir() if p.is_file() and p.suffix.lower() in LIGAND_SUFFIXES)
        input_file = input_files[0]
        input_ligands = _count_input(input_file)
        pdbqt_maximum = input_ligands * conf_poses * state_multiplier
        if input_file.suffix.lower() == ".csv":
            conf_command = f'"$PYBIN" "$CAMPAIGN_ROOT/1_ConformerGeneration.py" --mode 1 --csv "{input_file.name}" --smiles-col {_quote(csv_smiles)}'
            if csv_id:
                conf_command += f" --id-col {_quote(csv_id)}"
        elif input_file.suffix.lower() == ".sdf":
            conf_command = f'"$PYBIN" "$CAMPAIGN_ROOT/1_ConformerGeneration.py" --mode 2 --folder . --filetype sdf'
        else:
            conf_command = f'"$PYBIN" "$CAMPAIGN_ROOT/1_ConformerGeneration.py" --mode 2 --folder . --filetype smiles'
        conf_command += f" --num-confs {conf_poses} --workers {conf_workers} {state_flags}".rstrip() + "\n"
        # The board reads this fixed-size value rather than globbing a large
        # PDBQT directory on every scheduler poll.
        conf_command += 'LIGAND_DIR="$(find . -maxdepth 1 -type d -name "*_Ligands_PDBQT_*" | sort | tail -n 1)"\n[ -n "$LIGAND_DIR" ] || { echo "ConfGen PDBQT output was not found"; exit 2; }\nfind "$LIGAND_DIR" -maxdepth 1 -type f -name "*.pdbqt" -printf . | wc -c > confgen_pdbqt_count.txt\n'
        prelude = f'CAMPAIGN_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"\ncd "{shard_rel}"\nmkdir -p logs\n{common_setup}\nPYBIN="${{PYBIN:-python}}"\n'
        (shard / "run_confgen.lsf").write_text(_lsf_header(f"{campaign_name}-conf-{job_suffix}", queue, project, conf_workers, memory, settings.get("confgen_walltime", "48:00")) + prelude + conf_command)
        vina = _lsf_header(
            f"{campaign_name}-vina-{job_suffix}", queue, project, vina_workers, memory,
            settings.get("vina_walltime", "240:00"), exclusive_host=vina_exclusive_host,
            memory_limit_mb=vina_memory_limit_mb,
        ) + prelude
        vina += 'LIGAND_DIR="$(find . -maxdepth 1 -type d -name "*_Ligands_PDBQT_*" | sort | tail -n 1)"\n[ -n "$LIGAND_DIR" ] || { echo "ConfGen PDBQT output was not found"; exit 2; }\n'
        vina_options = (
            f"--exhaustiveness {int(settings.get('exhaustiveness', 8))} "
            f"--min-rmsd {float(settings.get('min_rmsd', 1.0))} "
            f"--compact-retain-top-poses {int(settings.get('retain_top_poses', 0))}"
        )
        if settings.get("energy_range") not in (None, ""):
            vina_options += f" --energy-range {float(settings['energy_range'])}"
        if settings.get("seed") not in (None, "", "0", 0):
            vina_options += f" --seed {int(settings['seed'])}"
        safety_options = (
            f' --max-rss-mb {vina_max_rss_mb}' if vina_max_rss_mb else ""
        ) + (
            f' --timeout-seconds {vina_timeout_seconds}' if vina_timeout_seconds else ""
        ) + f' --memory-trace "$CAMPAIGN_ROOT/results/{shard.name}/vina_memory_trace.tsv"'
        vina += f'"$PYBIN" "$CAMPAIGN_ROOT/3_Complete_batch_docking.py" --receptors "$CAMPAIGN_ROOT/{receptor_dir}" --ligands "$LIGAND_DIR" --centers_csv "$CAMPAIGN_ROOT/{center_rel}" --poses {vina_poses} --workers {vina_workers} --compact-results --results-dir "$CAMPAIGN_ROOT/results/{shard.name}" {vina_options}{safety_options}\n'
        (shard / "run_vina.lsf").write_text(vina)
        batch_progress.append({"name": shard.name, "input_ligands": input_ligands, "pdbqt_maximum": pdbqt_maximum, "vina_maximum": pdbqt_maximum * docking_targets, "docking_targets": docking_targets})

    manifest = {"campaign": campaign_name, "ligand_input": ligand_rel, "batch_count": len(shards), "batches": batch_progress, "settings": settings}
    (campaign / "campaign.json").write_text(json.dumps(manifest, indent=2) + "\n")
    submit = campaign / "submit_campaign.sh"
    submit.write_text(f'''#!/bin/bash
set -euo pipefail
CAMPAIGN_ROOT="$(cd "$(dirname "$0")" && pwd)"
MAX_PARALLEL="${{MAX_PARALLEL:-{int(settings.get("max_parallel", 3))}}}"
mkdir -p "$CAMPAIGN_ROOT/results"
submit_limited() {{
  while [ "$(bjobs -noheader -J "{campaign_name}-conf-*" 2>/dev/null | wc -l | tr -d ' ')" -ge "$MAX_PARALLEL" ]; do sleep 20; done
  bsub < "$1"
}}
for batch in "$CAMPAIGN_ROOT"/batches/{batch_prefix}_*; do
  [ -d "$batch" ] || continue
  conf_out="$(submit_limited "$batch/run_confgen.lsf")"
  conf_id="$(printf '%s' "$conf_out" | sed -n 's/.*<\\([0-9][0-9]*\\)>.*/\\1/p')"
  [ -n "$conf_id" ] || {{ echo "Could not parse ConfGen job ID: $conf_out" >&2; exit 2; }}
  bsub -w "done($conf_id)" < "$batch/run_vina.lsf"
done
''')
    submit.chmod(0o755)
    monitor = campaign / "monitor_campaign.sh"
    monitor.write_text('''#!/bin/bash
set -euo pipefail
CAMPAIGN_ROOT="$(cd "$(dirname "$0")" && pwd)"
exec "${PYBIN:-python}" "$CAMPAIGN_ROOT/campaign_status.py" --root "$CAMPAIGN_ROOT" --watch "${1:-15}" --write
''')
    monitor.chmod(0o755)
    return manifest | {"path": campaign}
