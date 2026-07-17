from __future__ import annotations

from datetime import datetime
from pathlib import Path
import os
import re

from hpc_profiles import HPCProfile, JOEY_LSF_PROFILE, render_lsf_header, render_setup_block, save_packaged_profile

DEFAULT_EMAIL = JOEY_LSF_PROFILE.email
DEFAULT_QUEUE = JOEY_LSF_PROFILE.queue
DEFAULT_PROJECT = JOEY_LSF_PROFILE.project
DEFAULT_WORKERS = JOEY_LSF_PROFILE.workers
DEFAULT_MEM_PER_CORE = JOEY_LSF_PROFILE.mem_per_core_mb


def _ts() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def sanitize_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def _chmod_executable(path: Path):
    try:
        mode = os.stat(path).st_mode
        os.chmod(path, mode | 0o111)
    except Exception:
        pass


def _header_with_timestamp(profile: HPCProfile, *, jobname: str, log_prefix: str, walltime: str) -> str:
    return f"#!/bin/bash\n# Auto-generated: {_ts()}\n" + render_lsf_header(
        profile=profile,
        jobname=jobname,
        log_prefix=log_prefix,
        walltime=walltime,
    ).split("\n", 1)[1]


def _python_export(profile: HPCProfile) -> str:
    return f'PYBIN="{profile.python_command}"\nif [ -z "$PYBIN" ]; then\n  echo "No Python command configured"; exit 127\nfi\necho "Using Python: $PYBIN"\n\n'


def _vina_export(profile: HPCProfile) -> str:
    if profile.vina_executable:
        return f'export VINA_EXE="{profile.vina_executable}"\n'
    return ""


CONFGEN_BODY = '"$PYBIN" 1_ConformerGeneration.py {flags} --num-confs {poses} --workers {workers}\n'

VINA_BODY = (
    '"$PYBIN" 3_Complete_batch_docking.py \\\n'
    '  --receptors "{receptors}" \\\n'
    '  --ligands   "{ligands}" \\\n'
    '  --centers_csv "{centers_csv}" \\\n'
    '  --poses {poses}\n'
)


def build_confgen_lsfs(
    jobroot: Path,
    lsf_dir: Path,
    *,
    profile: HPCProfile,
    poses: int,
    lig_mode,
    lig_filetype,
    csv_smiles_col,
    csv_id_col,
    single_sdf_rel,
):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)
    save_packaged_profile(jobroot, profile)

    jobname = sanitize_name(f"confgen_{jobroot.name}")
    header = _header_with_timestamp(
        profile,
        jobname=jobname,
        log_prefix=f"confgen_{jobroot.name}",
        walltime=profile.confgen_walltime,
    )
    setup_block = render_setup_block(profile)

    if lig_mode == "1":
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
        flags = f'--mode 3 --sdf "{single_sdf_rel or "Ligands/example.sdf"}"'
    else:
        flags = f'--mode 2 --folder "Ligands" --filetype {(lig_filetype or "sdf").lower()}'

    out = lsf_dir / "run_confgen_job.lsf"
    out.write_text(header + setup_block + _python_export(profile) + CONFGEN_BODY.format(flags=flags, poses=poses, workers=profile.workers))
    _chmod_executable(out)

    submit = lsf_dir / "submit_all_confgen.sh"
    submit.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    _chmod_executable(submit)


def build_vina_lsfs(
    jobroot: Path,
    lsf_dir: Path,
    *,
    profile: HPCProfile,
    poses: int,
):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)
    save_packaged_profile(jobroot, profile)

    rec_dir = (jobroot / "Receptors").name
    lig_dir = (jobroot / "Ligands").name
    centers_csv = next((p.name for p in jobroot.glob("vina_centers*.csv")), "vina_centers.csv")
    jobtag = sanitize_name(f"vina_{jobroot.name}")
    header = _header_with_timestamp(
        profile,
        jobname=jobtag,
        log_prefix=f"vina_{jobroot.name}",
        walltime=profile.vina_walltime,
    )

    body = VINA_BODY.format(
        receptors=rec_dir,
        ligands=lig_dir,
        centers_csv=centers_csv,
        poses=poses,
    )

    out = lsf_dir / "run_vina_job.lsf"
    out.write_text(header + render_setup_block(profile) + _vina_export(profile) + _python_export(profile) + body)
    _chmod_executable(out)

    submit = lsf_dir / "submit_all_vina.sh"
    submit.write_text(f"#!/bin/bash\nbsub < {out.name}\n")
    _chmod_executable(submit)
