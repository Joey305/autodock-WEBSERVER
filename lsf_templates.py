# ==============================
# lsf_templates.py
# ==============================
from __future__ import annotations
from pathlib import Path
from datetime import datetime
import re

def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)

HEADER = """#!/bin/bash
# Auto-generated: {timestamp}
#BSUB -J {jobname}
#BSUB -P {project}
#BSUB -o logs/{jobname}_%J.out
#BSUB -e logs/{jobname}_%J.err
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

# Conda env
{env_line}

{vina_pin}
"""

CONFGEN_BODY = (
    'PYBIN="$(command -v python3 || command -v python)"\n'
    'echo "Using Python: $PYBIN"\n\n'
    '"$PYBIN" 1_ConformerGeneration.py --mode {mode} {mode_flags} --num-confs {poses} --workers {workers}\n'
)

VINA_BODY = (
    'python 3_Complete_batch_docking.py \\\n  --receptors "{receptors}" \\\n  --ligands   "{ligands}" \\\n  --centers_csv "{centers_csv}" \\\n  --poses {poses}\n'
)


# lsf_templates.py
from pathlib import Path

def build_confgen_lsfs(jobroot: Path, lsf_dir: Path, *, poses:int, workers:int, queue:str,
                       project:str, walltime:str, mem_per_core:int, email:str, env_line:str,
                       lig_mode, lig_filetype, csv_smiles_col, csv_id_col, single_sdf_rel):
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)   # <-- required
    out = lsf_dir / "run_confgen_job.lsf"
    header = "...your header..."
    body   = "...your body..."
    out.write_text(header + body)

    """
    lig_mode:
      "1" CSV of SMILES in workspace root or Ligands/
      "2" Folder of ligands under Ligands/ (filetype sdf|smiles)
      "3" Single SDF relative path inside workspace (e.g., Ligands/foo.sdf)
    """
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    jobname = sanitize_name(f"confgen_{jobroot.name}")
    header = HEADER.format(
        timestamp=ts, jobname=jobname, project=project, walltime=walltime,
        queue=queue, workers=workers, mem_per_core=mem_per_core, email=email,
        env_line=env_line, vina_pin=""
    )

    if lig_mode == "1":
        # find CSV(s) under workspace/Ligands or ws root
        csvs = list((jobroot.parent).glob("*.csv")) + list((jobroot / "Ligands").glob("*.csv"))
        if not csvs:
            body = 'echo "No CSV found for mode=1"; exit 2\n'
        else:
            # Use the first CSV by default; users can edit later if needed
            csv_rel = csvs[0].relative_to(jobroot.parent)
            flags = f'--mode 1 --csv "{csv_rel}"'
            if csv_smiles_col: flags += f" --smiles-col {csv_smiles_col}"
            if csv_id_col: flags += f" --id-col {csv_id_col}"
            body = (
                'PYBIN="$(command -v python3 || command -v python)"\n'
                'echo "Using Python: $PYBIN"\n'
                f'"$PYBIN" 1_ConformerGeneration.py {flags} --num-confs {poses} --workers {workers}\n'
            )
    elif lig_mode == "3":
        srel = single_sdf_rel or "Ligands/example.sdf"
        body = (
            'PYBIN="$(command -v python3 || command -v python)"\n'
            'echo "Using Python: $PYBIN"\n'
            f'"$PYBIN" 1_ConformerGeneration.py --mode 3 --sdf "{srel}" --num-confs {poses} --workers {workers}\n'
        )
    else:
        # mode 2 = folder
        body = (
            'PYBIN="$(command -v python3 || command -v python)"\n'
            'echo "Using Python: $PYBIN"\n'
            f'"$PYBIN" 1_ConformerGeneration.py --mode 2 --folder Ligands --filetype {lig_filetype} '
            f'--num-confs {poses} --workers {workers}\n'
        )

    out = lsf_dir / f"run_{jobname}.lsf"
    out.write_text(header + body)
    sub = lsf_dir / "submit_all_confgen.sh"
    sub.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    sub.chmod(0o755)

from datetime import datetime

def build_vina_lsfs(jobroot: Path, lsf_dir: Path, poses: int, workers: int, queue: str,
                    project: str, walltime: str, mem_per_core: int, email: str,
                    env_line: str, vina_path: str | None = None):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)  # <-- Step 4b: ensure dir exists

    rec_dir = (jobroot / "Receptors").name
    lig_dir = (jobroot / "Ligands").name

    # Centers CSV default name (prefer vina_centers*.csv if present)
    centers_csv = next((p.name for p in jobroot.glob("vina_centers*.csv")), "vina_centers.csv")

    vina_pin = f'export VINA_EXE="{vina_path}"\n' if vina_path else ""

    jobtag = sanitize_name(f"vina_{jobroot.name}")
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header = HEADER.format(
        timestamp=ts, jobname=jobtag, project=project, walltime=walltime,
        queue=queue, workers=workers, mem_per_core=mem_per_core, email=email,
        env_line=env_line, vina_pin=vina_pin
    )
    body = VINA_BODY.format(
        receptors=rec_dir, ligands=lig_dir,
        centers_csv=centers_csv, poses=poses
    )

    out = lsf_dir / f"run_{jobtag}.lsf"
    out.write_text(header + body)

    submit = lsf_dir / "submit_all_vina.sh"
    submit.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    submit.chmod(0o755)