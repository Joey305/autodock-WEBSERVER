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


CSV_RESOLVER_BODY = r'''CSV="${CSV:-}"
if [ -z "$CSV" ]; then
  CSV="$(ls -1t ALL_Docking_Results_with_provenance_*.csv *vina_docking_scores_sorted.csv 2>/dev/null | head -n 1 || true)"
fi
if [ -z "$CSV" ] || [ ! -f "$CSV" ]; then
  echo "No docking score CSV found. Set CSV=/path/to/scores.csv or run 4_ParseScores.py/4C_ConcatenateScores.py first."
  exit 2
fi
echo "Using score CSV: $CSV"
'''

PYMOL_BODY = r'''"$PYBIN" 5C_BuildPymolSesh.py \
  --csv "$CSV" \
  --mode "${PYMOL_MODE:-per_receptor}" \
  --top "${PYMOL_TOP:-5}" \
  --receptor-roots ${RECEPTOR_ROOTS:-Receptors Receptors_PDBQT .} \
  --outdir "${PYMOL_OUTDIR:-.}" \
  --hydrogen-mode "${HYDROGEN_MODE:-none}" \
  --non-interactive ${PYMOL_EXTRA_ARGS:-}
'''

COMPACTED_SDF_HTML_BODY = r'''"$PYBIN" 5_COMPACTED_SDF_HTML.py \
  --csv "$CSV" \
  --outdir "${COMPACTED_OUTDIR:-.}" \
  --receptor-roots ${RECEPTOR_ROOTS:-Receptors Receptors_PDBQT .} \
  --top-ligands "${COMPACTED_TOP_LIGANDS:-25}" \
  --top-poses "${COMPACTED_TOP_POSES:-25}" \
  --project-name "${COMPACTED_PROJECT_NAME:-Docking_HTML_Viz_Project_Compacted_SDF}" \
  --page-title "${COMPACTED_PAGE_TITLE:-Compacted SDF Docking Visualization Project}" \
  --hydrogen-mode "${HYDROGEN_MODE:-none}" ${COMPACTED_EXTRA_ARGS:-}
'''


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


def _write_output_submitter(lsf_dir: Path):
    candidates = [
        "run_pymol_job.lsf",
        "run_compacted_sdf_html_job.lsf",
    ]
    existing = [name for name in candidates if (lsf_dir / name).is_file()]
    if not existing:
        return
    submit = lsf_dir / "submit_all_outputs.sh"
    lines = [
        "#!/bin/bash",
        "set -euo pipefail",
        f'echo "Submitting {len(existing)} output builder job(s)..."',
    ]
    for name in existing:
        lines.append(f"bsub < {name}")
    submit.write_text("\n".join(lines) + "\n")
    _chmod_executable(submit)


def build_pymol_lsf(
    jobroot: Path,
    lsf_dir: Path,
    *,
    profile: HPCProfile,
    mode: str = "per_receptor",
    top: int = 5,
):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)
    save_packaged_profile(jobroot, profile)

    jobtag = sanitize_name(f"pymol_{jobroot.name}")
    header = _header_with_timestamp(
        profile,
        jobname=jobtag,
        log_prefix=f"pymol_{jobroot.name}",
        walltime=profile.vina_walltime,
    )
    body = (
        f'export PYMOL_MODE="${{PYMOL_MODE:-{mode}}}"\n'
        f'export PYMOL_TOP="${{PYMOL_TOP:-{int(top)}}}"\n'
        + CSV_RESOLVER_BODY
        + PYMOL_BODY
    )

    out = lsf_dir / "run_pymol_job.lsf"
    out.write_text(header + render_setup_block(profile) + _python_export(profile) + body)
    _chmod_executable(out)
    _write_output_submitter(lsf_dir)


def build_compacted_sdf_html_lsf(
    jobroot: Path,
    lsf_dir: Path,
    *,
    profile: HPCProfile,
    top_ligands: int = 25,
    top_poses: int = 25,
):
    jobroot = Path(jobroot)
    lsf_dir = Path(lsf_dir)
    lsf_dir.mkdir(parents=True, exist_ok=True)
    save_packaged_profile(jobroot, profile)

    jobtag = sanitize_name(f"compact_sdf_{jobroot.name}")
    header = _header_with_timestamp(
        profile,
        jobname=jobtag,
        log_prefix=f"compact_sdf_{jobroot.name}",
        walltime=profile.vina_walltime,
    )
    body = (
        f'export COMPACTED_TOP_LIGANDS="${{COMPACTED_TOP_LIGANDS:-{int(top_ligands)}}}"\n'
        f'export COMPACTED_TOP_POSES="${{COMPACTED_TOP_POSES:-{int(top_poses)}}}"\n'
        + CSV_RESOLVER_BODY
        + COMPACTED_SDF_HTML_BODY
    )

    out = lsf_dir / "run_compacted_sdf_html_job.lsf"
    out.write_text(header + render_setup_block(profile) + _python_export(profile) + body)
    _chmod_executable(out)
    _write_output_submitter(lsf_dir)
