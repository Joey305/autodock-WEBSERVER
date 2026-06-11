# AutoDock-Vina Prep Portal

AutoDock-Vina Prep Portal is a public-access Flask application for preparing and packaging molecular docking jobs for local execution or HPC submission. It provides a guided workflow for receptor intake, docking box definition, ligand intake, and final job assembly, producing either a portable local-use ZIP or an optional Joey/LSF lab package with additional HPC scripts.

The application is designed to reduce repetitive setup work and make AutoDock Vina workflows more consistent, portable, and easier to hand off between users and environments.

---

## What it does

The portal helps you build a complete docking job package by walking through the major setup steps in order:

- Upload or fetch receptor structures
- Define docking centers in an interactive 3D viewer
- Stage receptor cleanup and preparation settings
- Upload ligand inputs
- Generate a final ZIP bundle containing inputs and execution scripts

The resulting package is suitable for workflows that use:

- AutoDock Vina
- Vina-compatible docking pipelines
- LSF-based HPC submission
- Local command-line execution of the included scripts

---

## Key features

### Guided receptor workflow
The portal supports both single-receptor and multi-receptor job setup. Receptors can be provided as:

- Single structure files
- ZIP archives of receptor folders
- Browser folder uploads
- PDB IDs fetched directly from RCSB

Supported receptor formats include:

- `.pdb`
- `.cif`
- `.mmcif`
- `.ent`
- `.pdbqt`

### Interactive box generation
Docking centers can be defined directly inside the built-in 3D viewer. Users can:

- Click one or more atoms and average their positions
- Paste explicit XYZ coordinates
- Use residue-based context while inspecting the structure

Centers are saved into a canonical CSV format:

```text
PDB_ID,X,Y,Z,SIZE
```

This file becomes part of the final packaged job.

### Protein staging and preparation
Before packaging, receptors can be staged with cleanup parameters such as:

- HETATM removal
- chain removal
- alternate location handling

Prepared receptor outputs are then used for downstream job assembly.

### Flexible ligand intake
Ligands can be uploaded as:

- a single `.sdf`
- a single `.smiles`
- a single `.smi`
- a single `.csv`
- a ZIP archive containing ligand files
- a browser folder upload containing ligand files

The downstream conformer-generation workflow supports multiple ligand input modes, including CSV-based SMILES input, folder-based SDF/SMILES collections, and single multi-record SDF files.
ZIP and folder uploads are normalized so supported ligand files land directly under `Ligands/` in the workspace and final package, even if the archive contains a top-level `Ligands/` folder or macOS metadata folders.

### Job bundle generation
The build step now supports two package modes:

- `portable` for public/local use
- `lsf` for Joey's Docking Server Mode / lab HPC use

Portable builds include the normalized runtime files, receptor and ligand inputs, center CSV, and local runner scripts. LSF builds include all portable content plus the extra LSF/batch scripts used in the lab environment.

### Public access mode
The application now runs without account registration or login. Ordinary users land directly on the main workflow and can create workspaces, upload inputs, prepare receptors, and build job bundles without credentials.

Workspaces are still isolated per generated job directory, but they are no longer keyed to authenticated user accounts. A minimal public identity is used internally only where the legacy code expects a user-like email value.

---

## Typical workflow

### 1. Protein Prep
Start by creating a workspace and loading receptor structures.

Available intake methods:

- Upload a single receptor file
- Upload a ZIP archive of receptors
- Upload a receptor folder from the browser
- Fetch a structure from the Protein Data Bank by PDB ID

For fetched structures, an optional chain filter can be provided.

### 2. Box Generation
Open a receptor in the viewer and define the docking center.

You can do this by:

- selecting atoms interactively in the 3D viewer
- averaging multiple picked atom positions
- pasting an explicit center as `x,y,z`

Each saved center is written into the centers CSV.

### 3. Protein Staging
Choose receptor cleanup options before preparation.

Examples include:

- removing all HET groups
- removing selected chains
- collapsing alternate locations

These settings are saved as part of the receptor preparation workflow.

### 4. Protein Preparation
Run receptor preparation so the docking-ready receptor files are available for packaging.

This stage is designed to fit naturally between center selection and final job assembly.

### 5. Ligand Prep
Upload the ligand input set for the job.

Accepted inputs include:

- single ligand files
- ligand ZIP archives
- ligand folder uploads
- CSV-based ligand collections

When the single ligand input is a CSV, the UI now reads the header row and offers dropdowns for:

- SMILES column
- ligand ID column

The SMILES column is required for CSV workflows. The ID column is optional and can be left blank if the source CSV does not have one.
For ZIP and folder uploads, the portal reports accepted ligand files and any ignored junk such as `__MACOSX/`, `.DS_Store`, `._*`, or unsupported file types.

The packaging step later passes the appropriate settings into the ligand conformer-generation workflow.

### 6. Build Job Package
Once receptors, centers, and ligands are in place, the portal builds the final job package.

In portable mode, the generated ZIP includes generic local runner scripts such as:

- `run_confgen_local.sh`
- `run_vina_local.sh`
- `run_all_local.sh`
- `README_RUN_LOCAL.md`

In Joey/LSF mode, the UI reveals the lab/HPC parameters such as:

- queue
- project
- walltime
- core count
- memory per core
- optional Vina executable path
- environment activation line

The generated ZIP can then be downloaded and transferred to the target execution environment.

---

## What the ZIP package contains

The final archive is built from the prepared workspace and includes the materials required to launch a docking run outside the web app.

Typical portable contents include:

- `Receptors/`
- `Receptors_PDB/`
- `Ligands/`
- `vina_centers.csv`
- `1_ConformerGeneration.py`
- `3_Complete_batch_docking.py`
- `4_ParseScores.py`
- `4C_ConcatenateScores.py`
- `5C_BuildPymolSesh.py`
- `6_MDpymacs.py`
- `7_Graphs.py`
- `ligand_naming.py`
- `ligand_manifest.py`
- `AutoDockTools_py3/`
- local runner scripts and local run README

LSF packages include all of the portable contents plus:

- `4B_LSFbatch.py`
- `3B_ServerDocks.py`
- `1B_confgen_batch.py`
- `3a_PDB2PDBQTbatch.py`
- `lsf_templates.py`
- `runDOCKING-tmux.sh`
- generated `.lsf` submission scripts

In practice, the package is intended to be execution-ready once moved into the appropriate environment with the expected software available.

---

## Included execution scripts

### Ligand conformer and PDBQT generation
The packaged ligand preparation utility supports three input modes:

1. CSV input containing SMILES plus an optional ligand ID column
2. One or more folders containing `.sdf` or `.smiles` files
3. A single multi-record SDF

The script uses RDKit and Open Babel to generate conformers and convert ligands into docking-ready PDBQT outputs. It also supports configurable pose counts, worker counts, and optional cleanup of temporary intermediate files.

### Batch AutoDock Vina runner
The packaged Vina runner is designed for batch docking across receptor and ligand sets. It expects:

- a receptor directory
- a ligand directory containing `.pdbqt` ligands
- a centers CSV with headers `PDB_ID,X,Y,Z`
- a specified number of Vina output poses

The runner writes:

- Vina config files
- docking result folders
- run logs
- summary files
- duration reports

It is written to work cleanly in scheduler-driven environments and limits each Vina process to one thread while using process-level parallelism across jobs.
It now checks the Vina executable up front before starting worker processes. The recommended setup is to make `vina --version` work in the target shell, or to set `VINA_EXE=/path/to/vina`.

### Post-docking score workflow
The cleaned post-docking parser workflow now uses:

- `4B_LSFbatch.py`
- `4_ParseScores.py`
- `4C_ConcatenateScores.py`

Retired parser variants are no longer part of the active workflow:

- `4A_Parse_VinaResults.py`
- `4A_PARSEvinascores.py`
- `4B_Parse_VinaResults.py`
- `4CMULTIVINA.py`
- `4F_Combine_PerDirCSVs_with_Provenance.py`
- `4_SERVERBATCH.py`

### Ligand naming and provenance
Generated ligand variants keep their full names through the workflow, for example:

- `obj01_p01_t01_c001`
- `compoundA_p02_t03_c032`
- `legacyLig_pose38`

The runtime scripts now preserve the full `LigandVariant` and add derived metadata columns such as:

- `LigandBase`
- `LigandVariant`
- `ProtomerTag`
- `TautomerTag`
- `ConformerTag`
- `LegacyPoseTag`
- `StateTag`

`Pose` remains the AutoDock Vina output pose index and is intentionally kept separate from the generated RDKit conformer tag such as `c001`.

### Ligand-state manifest
During conformer generation, `1_ConformerGeneration.py` writes `ligand_state_manifest.csv` into the generated ligand PDBQT folder. Each row describes one generated ligand variant and includes:

- `LigandBase`
- `LigandVariant`
- `ProtomerTag`
- `TautomerTag`
- `StateTag`
- `ConformerTag`
- canonical and isomeric SMILES
- state-level SMILES
- formula, charge, atom counts, rotatable bonds, and exact molecular weight
- source input tracking and generated file references

This manifest is optional for legacy workflows. If downstream scripts do not find it, they fall back to parsing the ligand filename so older `_pose#` and newer `_p##_t##_c###` naming schemes still work.

For SDF and SMILES file uploads, the primary ligand identity comes from the source filename. Internal SDF molecule titles are stored only as metadata such as `SourceRecordName`, `SDFTitle`, and `OriginalMolName`. This prevents files like `DR7.sdf`, `DR8.sdf`, and `DR9.sdf` from all collapsing into a shared internal title like `obj01`.
The same normalization also keeps ligand inputs directly under `Ligands/` in generated packages. If a package is manually edited into a nested `Ligands/Ligands/` layout, `1_ConformerGeneration.py --mode 2 --folder Ligands` now falls back to recursive ligand discovery so the folder still runs.

### Data model
- `LigandBase` = original molecule identity
- `ProtomerTag` = protonation/ionization state, for example `p01`
- `TautomerTag` = tautomer state, for example `t02`
- `StateTag` = protomer + tautomer, for example `p01_t02`
- `ConformerTag` = generated RDKit conformer, for example `c005`
- `LigandVariant` = full generated structure name, for example `obj01_p01_t02_c005`
- `Pose` = Vina docking pose index from `REMARK VINA RESULT`
- `Binding_Affinity` = Vina score

The public default `--max-tautomers` remains `4`. Increasing tautomer enumeration broadens chemical-state coverage but also multiplies the conformer and docking workload, so the default is intended as a practical balance for public/server use.

### Docking analytics and reports

`7_Graph.py` still produces the original per-ligand histograms, per-ligand boxplots, per-receptor combined boxplots, pooled histograms, and overlay histograms. It now also writes richer analytics for ranking and robustness review.

Key ligand-level metrics now include:

- `best_affinity` for the single strongest Vina score
- `median_affinity` for typical performance
- `top3_mean_affinity`, `top5_mean_affinity`, and `top10_mean_affinity` for stronger multi-pose ranking
- hit counts and hit rates at thresholds such as `<= -8`, `<= -8.5`, `<= -9`, `<= -9.5`, and `<= -10`
- `outlier_gap = median_affinity - best_affinity` to show whether the best score is a lucky outlier or close to typical behavior
- `consensus_rank`, which combines best-score, median-score, top-k mean, hit-rate, and outlier-gap rankings

More negative Vina affinity values are better.

The expanded script can also write:

- ligand, variant, tautomer, protomer, and state summary CSVs
- ranking bar charts
- consensus-rank charts
- hit-rate heatmaps
- state heatmaps
- `Plots/graph_report_manifest.json`
- an optional PDF report for human review

Recommended command:

```bash
python 7_Graph.py \
  --csv auto \
  --top-n 10 \
  --write-state-summaries \
  --write-ranking-plots \
  --write-hit-rate-plots \
  --write-state-heatmaps \
  --write-pdf-report
```

The PDF report is a convenient review artifact. The summary CSVs remain the authoritative machine-readable outputs.

### PyMOL top-hit exports

`5C_BuildPymolSesh.py` now reconstructs each selected hit by combining:

- the selected receptor structure
- the selected Vina pose coordinates from the chosen `MODEL` inside `out.pdbqt`
- the exact generated SDF for the same `LigandVariant`

This preserves the ligand's bond order, tautomer/protomer identity, and conformer identity from the generated SDF while still using the docked coordinates from Vina. The corrected ligand chemistry is written to `11_AA_Ligands/*_corrected_docked.sdf`; complex PDB files are coordinate views and should not be treated as the bond-order source of truth.

The PyMOL builder also applies a hydrogen cleanup policy before writing corrected docked ligands, loading ligand objects, and saving complex exports. The default is:

- `--hydrogen-mode nonpolar`

This removes all hydrogens and adds back only nonpolar hydrogens (defined here as hydrogens bonded to carbon atoms). Other options are:

- `--hydrogen-mode none`
- `--hydrogen-mode all`

For the cleanest session/export behavior, the recommended command is:

```bash
python 5C_BuildPymolSesh.py --hydrogen-mode nonpolar
```

The PyMOL export workflow also writes an audit CSV that links every saved output back to:

- `LigandVariant`
- `StateTag`
- the selected Vina `Pose`
- the exact reference SDF path
- the coordinate-transfer method used

---

## Installation

### Requirements

- Python 3.9+
- Flask application dependencies from `requirements.txt`
- A writable temporary workspace location
- Optional: RDKit
- Optional: Open Babel
- Optional: Tailscale for remote private access
- Optional: access to an LSF-based HPC environment for scheduler submission

### Environment setup

```bash
conda create -n docking python=3.10 -y
conda activate docking
pip install -r requirements.txt
```

Optional RDKit installation:

```bash
conda install -c conda-forge rdkit
```

Optional Open Babel installation:

```bash
conda install -c conda-forge openbabel
```

---

## Configuration

The application uses environment variables for its core configuration.

```bash
export PORTAL_SECRET="replace-with-a-secure-secret"
export PORTAL_DB="sqlite:///portal.db"
export PORTAL_TMP="/tmp/autodock_prep"
export PORTAL_MAX_RECEPTORS="5"
export PORTAL_ENV_LINE="source /path/to/conda.sh && conda activate vina_env"
export PORTAL_PUBLIC_EMAIL="public@autodock.local"
export PUBLIC_MODE="true"
export ENABLE_AUTH="false"
export ENABLE_LSF_PACKAGE="true"
export DEFAULT_PACKAGE_MODE="portable"
```

### Configuration notes

- `PORTAL_SECRET`  
  Secret key used by Flask.

- `PORTAL_DB`  
  SQLAlchemy database URI. SQLite is suitable for local use; a production database is recommended for multi-user deployment.

- `PORTAL_TMP`  
  Root directory for user workspaces and generated artifacts.

- `PORTAL_MAX_RECEPTORS`  
  Maximum number of receptors accepted per workspace.

- `PORTAL_ENV_LINE`  
  Default environment activation line inserted into generated job scripts.

- `PORTAL_PUBLIC_EMAIL`  
  Optional public identity email used in places where the legacy workflow expected a user email. Defaults to `public@autodock.local`.

- `PUBLIC_MODE`
  Enables the public no-auth workflow. Defaults to `true`.

- `ENABLE_AUTH`
  Reserved compatibility flag for deployments that still need account-based auth behavior. The default public setup keeps this `false`.

- `ENABLE_LSF_PACKAGE`
  Enables Joey's Docking Server Mode / LSF packaging in the UI and build endpoint.

- `DEFAULT_PACKAGE_MODE`
  Chooses the default build mode shown in the UI. Supported values are `portable` and `lsf`.

### Example public configuration

```bash
PUBLIC_MODE=true
ENABLE_AUTH=false
ENABLE_LSF_PACKAGE=false
DEFAULT_PACKAGE_MODE=portable
```

### Example Joey/lab configuration

```bash
PUBLIC_MODE=true
ENABLE_AUTH=false
ENABLE_LSF_PACKAGE=true
DEFAULT_PACKAGE_MODE=lsf
```

---

## Running the application

Start the server with:

```bash
python app.py
```

Or run it with Flask directly:

```bash
flask run -h 0.0.0.0 -p 5050
```

Then open:

```text
http://localhost:5050
```

The app opens directly in public no-auth mode. No sign-in step is required.

Portable ZIPs require the receiving environment to provide its own Python, AutoDock Vina, Open Babel, and any relevant Python dependencies.

---

## Remote access

For private remote access, the portal can also be exposed through a Tailscale network.

Example:

```bash
sudo tailscale up
```

Then connect through the machine’s Tailscale IP and configured application port.

---

## Input expectations

### Receptors
Recommended for visualization:

- `.pdb`
- `.cif`
- `.mmcif`
- `.ent`

Accepted for docking-oriented workflows:

- `.pdbqt`

Note: PDBQT inputs may be usable for packaging, but structure visualization is centered around standard structure formats.

### Ligands
Supported ligand inputs include:

- `.sdf`
- `.smiles`
- `.csv`
- ZIP archives containing ligand files

### Centers CSV
The canonical center format is:

```text
PDB_ID,X,Y,Z,SIZE
```

Example:

```csv
PDB_ID,X,Y,Z,SIZE
receptor_A.pdbqt,12.4,18.2,5.9,20
receptor_B.pdbqt,-4.8,10.1,27.3,20
```

---

## Architecture overview

The portal is built with:

- Flask
- SQLAlchemy
- Bootstrap-based templates
- NGL Viewer for interactive structure inspection

At a high level, the application:

1. creates a per-job workspace
2. stores receptor and ligand inputs
3. records docking center definitions in a canonical CSV
4. assembles a job tree from the prepared workspace
5. generates either portable local runner scripts or optional LSF scripts
6. zips the complete job for download

---

## Production considerations

For shared or production deployment:

- run the app behind a production WSGI server such as Gunicorn or uWSGI
- place it behind a reverse proxy such as Nginx
- use a production-grade database
- store secrets securely
- ensure the temporary workspace directory is writable and appropriately managed
- add rate limits, upload limits, workspace cleanup, and anonymous job quotas for public deployments

---

## Troubleshooting

### Viewer does not show the structure
Use a structure format intended for visualization, such as `.pdb`, `.cif`, `.mmcif`, or `.ent`. Packaging may still work with `.pdbqt`, but viewer behavior is best with standard structure files.

### PDB fetch fails
Verify:

- the PDB ID is valid
- outbound network access is available
- any chain filter is provided in a comma-separated format such as `A,B`

### Build step says centers are incomplete
Every expected receptor must have a saved entry in the centers CSV before the final ZIP can be generated.

### Build step says ligands are missing
At least one ligand input must be uploaded before packaging can proceed.

### Open Babel is not found
Install Open Babel or make sure `obabel` is available on `PATH`. The ligand preparation scripts depend on it for file conversion.

### RDKit is not installed
Install RDKit if you plan to use the bundled ligand conformer-generation workflow.

---

## Scope

This application prepares and packages docking jobs. It does not replace the external software environment required to actually run those jobs. The final ZIP is intended to bridge the gap between web-based setup and command-line or scheduler-based execution.

---

## License

MIT
