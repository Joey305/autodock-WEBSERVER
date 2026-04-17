# ==============================
# lsf_templates.py
# ==============================
from __future__ import annotations
from pathlib import Path
from datetime import datetime
import re
import os

DEFAULT_EMAIL = "jxs794@miami.edu"
DEFAULT_QUEUE = "hihg"
DEFAULT_PROJECT = "brd"
DEFAULT_WORKERS = 16
DEFAULT_MEM_PER_CORE = 2000

# -------- utils --------
def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)

def _chmod_executable(path: Path):
    """Ensure a file is executable."""
    try:
        mode = os.stat(path).st_mode
        os.chmod(path, mode | 0o111)
    except Exception:
        pass

HEADER = """#!/bin/bash
# Auto-generated: {timestamp}
#BSUB -J {jobname}
#BSUB -P {project}
#BSUB -o logs/{log_prefix}_%J.out
#BSUB -e logs/{log_prefix}_%J.err
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

{env_block}{vina_pin}
"""

# Body snippets
CONFGEN_BODY = (
    'PYBIN="$(command -v python3 || command -v python)"\n'
    'echo "Using Python: $PYBIN"\n\n'
    '"$PYBIN" 1_ConformerGeneration.py {flags} --num-confs {poses} --workers {workers}\n'
)

VINA_BODY = (
    'python 3_Complete_batch_docking.py \\\n'
    '  --receptors "{receptors}" \\\n'
    '  --ligands   "{ligands}" \\\n'
    '  --centers_csv "{centers_csv}" \\\n'
    '  --poses {poses}\n'
)

# -------- ConfGen builder --------
def build_confgen_lsfs(jobroot: Path,
                       lsf_dir: Path,
                       *,
                       poses: int,
                       workers: int,
                       queue: str,
                       project: str,
                       walltime: str,
                       mem_per_core: int,
                       email: str,
                       env_line: str,
                       lig_mode,                 # "1" CSV, "2" folder, "3" single SDF
                       lig_filetype,            # None | "sdf" | "smiles"
                       csv_smiles_col,
                       csv_id_col,
                       single_sdf_rel):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)

    # Default env block if UI left it blank
    if env_line and env_line.strip():
        env_block = f"# Conda env (from UI)\n{env_line.strip()}\n"
    else:
        env_block = (
            "# Activate environment\n"
            "source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh && conda activate vina_env\n"
        )

    jobname = sanitize_name(f"confgen_{jobroot.name}")
    ts = _ts()
    header = HEADER.format(
        timestamp=ts,
        jobname=jobname,
        log_prefix=f"confgen_{jobroot.name}",
        project=project or DEFAULT_PROJECT,
        walltime=walltime or "108:00",
        queue=queue or DEFAULT_QUEUE,
        workers=workers or DEFAULT_WORKERS,
        mem_per_core=mem_per_core or DEFAULT_MEM_PER_CORE,
        email=email or DEFAULT_EMAIL,
        env_block=env_block,
        vina_pin=""  # none for confgen
    )

    # Build flags by mode
    if lig_mode == "1":
        # CSV mode: try workspace and Ligands
        csvs = list(jobroot.parent.glob("*.csv")) + list((jobroot / "Ligands").glob("*.csv"))
        if not csvs:
            flags = '--mode 1 --csv "MISSING.csv"'
        else:
            csv_rel = csvs[0].relative_to(jobroot.parent)
            flags = f'--mode 1 --csv "{csv_rel}"'
            if csv_smiles_col:
                flags += f" --smiles-col {csv_smiles_col}"
            if csv_id_col:
                flags += f" --id-col {csv_id_col}"

    elif lig_mode == "3":
        # Single SDF path relative to workspace root
        srel = single_sdf_rel or "Ligands/example.sdf"
        flags = f'--mode 3 --sdf "{srel}"'

    else:
        # Mode 2: folder; force filetype to sdf if not given
        ft = (lig_filetype or "sdf").lower()
        # Default folder name inside job is "Ligands"
        flags = f'--mode 2 --folder "Ligands" --filetype {ft}'

    body = CONFGEN_BODY.format(flags=flags, poses=poses, workers=workers)

    out = lsf_dir / "run_confgen_job.lsf"
    out.write_text(header + body)
    _chmod_executable(out)

    submit = lsf_dir / "submit_all_confgen.sh"
    submit.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    _chmod_executable(submit)

# -------- Vina builder --------
def build_vina_lsfs(jobroot: Path,
                    lsf_dir: Path,
                    poses: int,
                    workers: int,
                    queue: str,
                    project: str,
                    walltime: str,
                    mem_per_core: int,
                    email: str,
                    env_line: str,
                    vina_path: str | None = None):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)

    # Default env block if UI left it blank (multi-line)
    if env_line and env_line.strip():
        env_block = f"# Conda env (from UI)\n{env_line.strip()}\n"
    else:
        env_block = (
            "# Conda env with vina\n"
            "source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh\n"
            "conda activate vina_env\n"
        )

    # Default VINA_EXE export if not provided
    if vina_path and str(vina_path).strip():
        vina_pin = f'export VINA_EXE="{str(vina_path).strip()}"\n'
    else:
        vina_pin = 'export VINA_EXE="$HOME/miniconda3/envs/vina_env/bin/vina"\n'

    rec_dir = (jobroot / "Receptors").name
    lig_dir = (jobroot / "Ligands").name
    centers_csv = next((p.name for p in jobroot.glob("vina_centers*.csv")), "vina_centers.csv")

    jobtag = sanitize_name(f"vina_{jobroot.name}")
    ts = _ts()
    header = HEADER.format(
        timestamp=ts,
        jobname=jobtag,
        log_prefix=f"vina_{jobroot.name}",
        project=project or DEFAULT_PROJECT,
        walltime=walltime or "96:00",
        queue=queue or DEFAULT_QUEUE,
        workers=workers or DEFAULT_WORKERS,
        mem_per_core=mem_per_core or DEFAULT_MEM_PER_CORE,
        email=email or DEFAULT_EMAIL,
        env_block=env_block,
        vina_pin=vina_pin
    )

    body = VINA_BODY.format(
        receptors=rec_dir,
        ligands=lig_dir,
        centers_csv=centers_csv,
        poses=poses
    )

    out = lsf_dir / "run_vina_job.lsf"
    out.write_text(header + body)
    _chmod_executable(out)

    submit = lsf_dir / "submit_all_vina.sh"
    submit.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    _chmod_executable(submit)
