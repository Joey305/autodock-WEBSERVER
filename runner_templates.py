from __future__ import annotations

import os
from pathlib import Path


PORTABLE_README = """# Local Run Guide

This package was built in portable mode for local or generic environment use.

Requirements:
- `python3`
- `obabel` available on `PATH` or set with `OBABEL_BIN`
- `vina` available on `PATH` or set with `VINA_EXE`
- the Python dependencies required by the bundled scripts

Environment setup:
- Use `./create_vina_env.sh` to create a conda environment named `vina_env` from the bundled `docking.yaml`
- If the solver fails while installing `vina`, the script retries by creating the environment first and installing `vina` separately

Recommended workflow:
0. Optional: run `python 0_LIGSPLIT.py` first if you want to break a very large ligand library into batch-sized folders.
1. Review `run_confgen_local.sh` and adjust the variables near the top if needed.
2. Run `./run_confgen_local.sh`.
3. Review `run_vina_local.sh` and adjust the variables if needed.
4. Run `./run_vina_local.sh`.
5. Parse or combine docking scores with `4_ParseScores.py` and `4C_ConcatenateScores.py` as needed.

Naming and provenance:
- For SDF and SMILES file uploads, ligand identity is taken from the input filename
- Internal SDF titles are kept only as metadata such as `SourceRecordName`
- Ligand ZIP and folder uploads are normalized so supported files land directly under `Ligands/`
- If you manually create a nested `Ligands/Ligands/` layout, `1_ConformerGeneration.py` now falls back to recursive ligand discovery
- `1_ConformerGeneration.py` preserves generated names like `<LigandBase>_p01_t01_c001`
- `1_ConformerGeneration.py` also writes `ligand_state_manifest.csv` beside the generated PDBQT files
- `3_Complete_batch_docking.py` keeps the full ligand variant as the docking output folder name
- `4_ParseScores.py` adds `LigandBase`, `LigandVariant`, state tags, conformer/protomer/tautomer indexes, and chemical metadata from the manifest when available
- `4C_ConcatenateScores.py` preserves those fields in the combined CSV
- `7_Graphs.py --write-state-summaries` writes ligand, tautomer, protomer, state, and variant summary CSVs
- `7_Graphs.py` also supports richer docking analytics such as best-score, median-score, top-k mean, hit-rate, outlier-gap, and consensus-rank summaries
- Use `7_Graphs.py --write-ranking-plots --write-hit-rate-plots --write-state-heatmaps --write-pdf-report` to generate ranking charts, heatmaps, summary tables, a run manifest JSON, and an optional PDF review report
- `5C_BuildPymolSesh.py` now requires an exact generated SDF for each `LigandVariant` by default, prefers `Ligands_TMP_SDF_*` sources, excludes previous `11_AA_*` outputs from reference lookup, and writes corrected docked SDF exports plus an audit CSV
- `5C_BuildPymolSesh.py` applies a hydrogen policy before PyMOL/export; the default is `--hydrogen-mode nonpolar`
- Use `--hydrogen-mode none` for fully dehydrogenated display/export molecules, or `--hydrogen-mode all` to keep all hydrogens
- Base/source SDF fallback is disabled unless you pass `--allow-reference-fallback`
- Partial MCS coordinate transfer is disabled unless you pass `--allow-partial-mcs`
- Recommended PyMOL preflight: `python 5C_BuildPymolSesh.py --audit-reference-resolution --strict-reference-sdf`
- Recommended full export: `python 5C_BuildPymolSesh.py --strict-reference-sdf --write-atom-map-audit --hydrogen-mode nonpolar`
- Complex PDB exports are coordinate views only; preserved ligand chemistry lives in `11_AA_Ligands/*_corrected_docked.sdf`
- `Pose` means the Vina output pose number, not the RDKit conformer index

Vina setup:
- Prefer making `vina --version` work in your shell
- Or set `VINA_EXE=/path/to/vina`
- The runtime scripts treat Vina as an executable dependency; they do not require the Python `vina` package

`run_all_local.sh` runs both steps sequentially.
"""


RUN_CONFGEN = """#!/bin/bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
OBABEL_BIN="${OBABEL_BIN:-obabel}"
POSES="${POSES:-64}"
WORKERS="${WORKERS:-4}"
LIG_MODE="${LIG_MODE:-2}"
LIG_FOLDER="${LIG_FOLDER:-Ligands}"
LIG_FILETYPE="${LIG_FILETYPE:-sdf}"
CSV_FILE="${CSV_FILE:-Ligands/input.csv}"
CSV_SMILES_COL="${CSV_SMILES_COL:-}"
CSV_ID_COL="${CSV_ID_COL:-}"
SINGLE_SDF="${SINGLE_SDF:-Ligands/input.sdf}"

export OBABEL_BIN

case "$LIG_MODE" in
  1)
    ARGS=(--mode 1 --csv "$CSV_FILE")
    if [[ -n "$CSV_SMILES_COL" ]]; then ARGS+=(--smiles-col "$CSV_SMILES_COL"); fi
    if [[ -n "$CSV_ID_COL" ]]; then ARGS+=(--id-col "$CSV_ID_COL"); fi
    ;;
  3)
    ARGS=(--mode 3 --sdf "$SINGLE_SDF")
    ;;
  *)
    ARGS=(--mode 2 --folder "$LIG_FOLDER" --filetype "$LIG_FILETYPE")
    ;;
esac

"$PYTHON_BIN" 1_ConformerGeneration.py "${ARGS[@]}" --num-confs "$POSES" --workers "$WORKERS"
"""


RUN_VINA = """#!/bin/bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python3}"
VINA_EXE="${VINA_EXE:-vina}"
POSES="${POSES:-9}"
RECEPTORS_DIR="${RECEPTORS_DIR:-Receptors}"
LIGANDS_DIR="${LIGANDS_DIR:-Ligands}"
CENTERS_CSV="${CENTERS_CSV:-vina_centers.csv}"

export VINA_EXE

"$PYTHON_BIN" 3_Complete_batch_docking.py \\
  --receptors "$RECEPTORS_DIR" \\
  --ligands "$LIGANDS_DIR" \\
  --centers_csv "$CENTERS_CSV" \\
  --poses "$POSES" \\
  --vina-exe "$VINA_EXE"
"""


RUN_ALL = """#!/bin/bash
set -euo pipefail

./run_confgen_local.sh
./run_vina_local.sh
"""


def _write_executable(path: Path, content: str):
    path.write_text(content)
    mode = os.stat(path).st_mode
    os.chmod(path, mode | 0o111)


def build_portable_runners(jobroot: Path):
    jobroot = Path(jobroot)
    _write_executable(jobroot / "run_confgen_local.sh", RUN_CONFGEN)
    _write_executable(jobroot / "run_vina_local.sh", RUN_VINA)
    _write_executable(jobroot / "run_all_local.sh", RUN_ALL)
    (jobroot / "README_RUN_LOCAL.md").write_text(PORTABLE_README)
