# AutoDock-Vina Prep Portal

AutoDock-Vina Prep Portal is a login-protected Flask application for preparing and packaging molecular docking jobs for local execution or HPC submission. It provides a guided workflow for receptor intake, docking box definition, ligand intake, and final job assembly, producing a downloadable bundle that includes the required inputs, scheduler scripts, and batch docking utilities.

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
- a single `.csv`
- a ZIP archive containing ligand files

The downstream conformer-generation workflow supports multiple ligand input modes, including CSV-based SMILES input, folder-based SDF/SMILES collections, and single multi-record SDF files.

### Job bundle generation
The build step packages the working job into a downloadable ZIP that includes:

- receptor inputs
- ligand inputs
- center CSV
- generated LSF submission scripts
- ligand conformer/PDBQT generation utilities
- batch AutoDock Vina execution utilities

This allows the prepared job to be moved to another machine or submitted to an HPC environment with minimal manual editing.

### Authentication and per-user workspaces
The application uses authenticated sessions and creates isolated workspaces for each user job. This keeps uploaded inputs and generated artifacts separated across runs.

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
- CSV-based ligand collections

The packaging step later passes the appropriate settings into the ligand conformer-generation workflow.

### 6. Build LSF + ZIP
Once receptors, centers, and ligands are in place, the portal builds the final job package.

The UI exposes job-level parameters such as:

- queue
- project
- walltime
- core count
- memory per core
- environment activation line

The generated ZIP can then be downloaded and transferred to the target execution environment.

---

## What the ZIP package contains

The final archive is built from the prepared workspace and includes the materials required to launch a docking run outside the web app.

Typical contents include:

- `Receptors/`
- `Ligands/`
- `vina_centers*.csv`
- generated `.lsf` submission scripts
- ligand conformer-generation and PDBQT conversion script(s)
- batch AutoDock Vina runner script(s)

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

---

## Installation

### Requirements

- Python 3.10+
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

Sign in with an existing user account. For a fresh deployment, initialize your authentication flow according to your project’s user-creation method.

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
- Flask-Login
- SQLAlchemy
- Bootstrap-based templates
- NGL Viewer for interactive structure inspection

At a high level, the application:

1. creates a per-job workspace
2. stores receptor and ligand inputs
3. records docking center definitions in a canonical CSV
4. assembles a job tree from the prepared workspace
5. generates LSF scripts for conformer generation and Vina
6. zips the complete job for download

---

## Production considerations

For shared or production deployment:

- run the app behind a production WSGI server such as Gunicorn or uWSGI
- place it behind a reverse proxy such as Nginx
- use a production-grade database
- store secrets securely
- ensure the temporary workspace directory is writable and appropriately managed

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
