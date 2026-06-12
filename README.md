# AutoDock-Vina PrepServer

<p align="center">
  <strong>AutoDock-Vina PrepServer: A browser and API-driven workspace builder for reproducible AutoDock Vina docking packages</strong>
</p>

<p align="center">
  <em>Create receptor workspaces, define docking boxes, prepare PDBQT receptors, upload ligand libraries, and export ready-to-run AutoDock Vina job archives.</em>
</p>

<p align="center">
  <a href="https://autodockvina.com">
    <img src="https://img.shields.io/badge/Open%20Web%20Server-autodockvina.com-success?style=for-the-badge&logo=googlechrome" alt="Open AutoDock-Vina PrepServer">
  </a>
  <a href="https://github.com/Joey305/autodock-WEBSERVER">
    <img src="https://img.shields.io/badge/View%20Repository-GitHub-black?style=for-the-badge&logo=github" alt="View GitHub repository">
  </a>
  <a href="#environment-setup">
    <img src="https://img.shields.io/badge/Install-Conda%20Environment-yellow?style=for-the-badge&logo=anaconda" alt="Install conda environment">
  </a>
  <a href="#headless-api-tutorial">
    <img src="https://img.shields.io/badge/API%20Tutorial-curl%20Workflow-blue?style=for-the-badge&logo=gnubash" alt="API tutorial">
  </a>
  <a href="#real-example-walkthrough">
    <img src="https://img.shields.io/badge/Run%20Example-DR7%20Ligand-orange?style=for-the-badge&logo=readthedocs" alt="Run DR7 example">
  </a>
</p>

<p align="center">
  <a href="mailto:jmschulz@med.miami.edu?subject=AutoDock-Vina%20PrepServer%20Question">
    <img src="https://img.shields.io/badge/Contact-Joseph%20M.%20Schulz-blue?style=for-the-badge&logo=gmail" alt="Contact Joseph M. Schulz">
  </a>
  <a href="https://schurerlab.org">
    <img src="https://img.shields.io/badge/Sch%C3%BCrer%20Lab-Molecular%20Design-lightgrey?style=for-the-badge" alt="Schürer Lab">
  </a>
</p>

---

## Overview

**AutoDock-Vina PrepServer** is a Flask-based web server and headless API for assembling reproducible molecular docking packages.

The server helps users move from raw receptor and ligand inputs to a structured, downloadable AutoDock Vina job archive. It is designed for researchers, students, and computational chemistry workflows where the preparation steps should be visible, reproducible, and easier to debug than a single black-box endpoint.

The public deployment is available at:

```text
https://autodockvina.com
```

The web interface and API follow the same staged workflow:

```text
Create workspace
→ Add receptor
→ Save docking center
→ Prepare receptor PDBQT
→ Upload ligand
→ Build package
→ Download ZIP
→ Run package locally or on compute
```

---

## Why this exists

Docking workflows often fail for practical reasons before docking even begins:

* receptor files are inconsistently named
* docking centers are not mapped to the correct receptor
* PDBQT conversion fails silently
* ligand libraries are uploaded in mixed formats
* AutoDock Vina is not available in the execution environment
* HPC or local execution folders are assembled by hand
* output ZIPs are difficult to reproduce later

AutoDock-Vina PrepServer makes these preparation stages explicit. Each step has a visible state and a scriptable API endpoint so failures can be diagnosed before running docking.

---

## Key Features

* Browser-based docking package preparation
* Versioned headless API under `/api/v1`
* Workspace-based file organization
* Receptor upload or RCSB PDB fetch
* Docking box center saving by XYZ coordinates or resolver-supported methods
* Open Babel-based receptor PDBQT preparation
* Ligand upload from `.sdf`, `.smi`, `.smiles`, `.csv`, ZIP, or folder-style inputs
* Portable package generation for local/non-HPC workflows
* Optional LSF-oriented package generation for scheduler-based environments
* Downloadable ZIP artifacts
* JSON status endpoints for workflow automation
* Clear error responses for missing receptors, incomplete centers, missing ligands, missing Open Babel, or conversion failures
* Documented Vina CLI installation path for running downloaded packages

---

## Companion Molecular Design Tools

AutoDock-Vina PrepServer is part of a growing structure-guided molecular design and analysis ecosystem.

| Tool                     | Website                    | Focus                                                              |
| ------------------------ | -------------------------- | ------------------------------------------------------------------ |
| AutoDock-Vina PrepServer | https://autodockvina.com   | Receptor/ligand intake and docking package preparation             |
| Warhead Hunter           | https://warheadhunter.com  | Solvent-exposed ligand atom and warhead/linker follow-up           |
| PROTAC Builder           | https://protacbuilder.com  | Degrader design continuation and linker/recruiter/warhead assembly |
| E3 Ligandalyzer          | https://e3ligandalyzer.com | E3 recruiter and ligase-context exploration                        |
| V-LiSEMOD                | https://vlisemod.com       | Viral ligand solvent-exposed moiety and PROTACability-style triage |

---

## Repository Navigation

<p align="center">
  <a href="#quick-start-web-interface">
    <img src="https://img.shields.io/badge/Quick%20Start-Web%20Interface-success?style=for-the-badge&logo=googlechrome" alt="Web interface quick start">
  </a>
  <a href="#environment-setup">
    <img src="https://img.shields.io/badge/Install-docking.yaml-yellow?style=for-the-badge&logo=anaconda" alt="Install environment">
  </a>
  <a href="#manual-vina-binary-install">
    <img src="https://img.shields.io/badge/Vina-Manual%20Binary%20Install-purple?style=for-the-badge&logo=github" alt="Manual Vina binary install">
  </a>
  <a href="#headless-api-tutorial">
    <img src="https://img.shields.io/badge/Tutorial-Headless%20API-blue?style=for-the-badge&logo=gnubash" alt="Headless API tutorial">
  </a>
  <a href="#real-example-walkthrough">
    <img src="https://img.shields.io/badge/Example-DR7%20SDF-orange?style=for-the-badge&logo=readthedocs" alt="DR7 example">
  </a>
</p>

* [Quick Start: Web Interface](#quick-start-web-interface)
* [Repository Structure](#repository-structure)
* [Environment Setup](#environment-setup)
* [Manual Vina Binary Install](#manual-vina-binary-install)
* [Running Locally](#running-locally)
* [Headless API Tutorial](#headless-api-tutorial)
* [Real Example Walkthrough](#real-example-walkthrough)
* [What is inside the downloaded ZIP?](#what-is-inside-the-downloaded-zip)
* [Endpoint Reference](#endpoint-reference)
* [Accepted Input Types](#accepted-input-types)
* [Generated Package Modes](#generated-package-modes)
* [Deployment Notes](#deployment-notes)
* [Heroku + Open Babel Notes](#heroku--open-babel-notes)
* [Testing and Validation](#testing-and-validation)
* [Troubleshooting](#troubleshooting)
* [Security and Data Notes](#security-and-data-notes)
* [Citation and Acknowledgements](#citation-and-acknowledgements)
* [Contact](#contact)

---

## Repository Structure

```text
.
├── app.py                         # Flask application and API routes
├── packager.py                    # Workspace assembly and ZIP packaging helpers
├── center_resolver.py             # Docking center resolution helpers
├── runner_templates.py            # Portable runner script templates
├── lsf_templates.py               # LSF package-generation templates
├── ligand_manifest.py             # Ligand provenance and upload metadata helpers
├── requirements.txt               # Python dependencies for pip-style installs
├── docking.yaml                   # Recommended conda environment
├── Aptfile                        # Optional Heroku apt packages, including Open Babel
├── templates/                     # Flask HTML templates
├── static/                        # CSS, JS, images, favicon assets
├── tests/                         # Regression tests
├── Example/
│   └── Ligand_SDF/
│       └── DR7.sdf                # Example ligand used in the walkthrough below
├── AutoDockTools_py3/             # Bundled AutoDockTools source tree
├── 1_ConformerGeneration.py       # Ligand conformer/PDBQT generation
├── 2a_PDB2PDBQTbatch.py           # Receptor preparation utility
├── 3_Complete_batch_docking.py    # AutoDock Vina batch runner
├── 4_ParseScores.py               # Vina score parsing
├── 4B_LSFbatch.py                 # LSF-oriented batch helper
├── 4C_ConcatenateScores.py        # Score collation
├── 5C_BuildPymolSesh.py           # PyMOL/session export helper
├── 6_MDpymacs.py                  # Additional downstream analysis helper
├── 7_Graphs.py                    # Plotting/ranking summaries
└── README.md
```

---

<a id="quick-start-web-interface"></a>

## Quick Start: Web Interface

Open:

```text
https://autodockvina.com
```

Typical browser workflow:

1. Create a workspace.
2. Upload a receptor or fetch one by PDB ID.
3. Define a docking center.
4. Prepare the receptor.
5. Upload ligand input.
6. Build a package.
7. Download the generated ZIP.
8. Run the downloaded docking package locally or on your target compute environment.

The browser interface is recommended for first-time users because each workflow stage is visible.

---

<a id="environment-setup"></a>

## Environment Setup

This repository includes a conda environment file:

```text
docking.yaml
```

Create the environment:

```bash
conda env create -f docking.yaml
conda activate docking
```

Verify the main Python stack:

```bash
python --version
python -c "import flask, pandas, numpy, matplotlib; print('Core Python imports OK')"
python -c "from rdkit import Chem; print('RDKit OK')"
```

Verify Open Babel:

```bash
which obabel
obabel -V
```

Verify Vina:

```bash
which vina
vina --version
```

If `vina --version` works, the environment is ready to run downloaded docking packages.

If Vina does not install cleanly through conda, use the manual binary install below.

---

## Recommended `docking.yaml`

The repository environment can be created from this file:

```yaml
name: docking

channels:
  - conda-forge
  - defaults

dependencies:
  # Core language and Python tools
  - python=3.9
  - pip

  # Molecular modeling / structural biology binaries
  - openbabel=3.1.1
  - pymol-open-source=3.1.0
  - pyqt=5.15.11

  # AutoDock Vina CLI binary
  # If this causes solver issues on a specific machine, create the env first,
  # then install or replace Vina manually using the binary install section below.
  - vina

  # Scientific Python stack
  - numpy
  - pandas
  - matplotlib
  - rdkit

  # Flask web stack via pip
  - pip:
    - flask==3.0.3
    - flask-login==0.6.3
    - flask-wtf==1.2.1
    - wtforms==3.1.2
    - flask-sqlalchemy==3.1.1
    - sqlalchemy==2.0.32
    - email-validator==2.1.1
    - requests==2.32.3
    - gunicorn==23.0.0
```

### Why Vina has a separate fallback

The PrepServer web app prepares and packages docking jobs. The downloaded package expects an AutoDock Vina executable when the user runs docking locally or on a compute server.

In other words:

```text
PrepServer builds the docking package.
Vina runs the docking job after the package is downloaded.
```

On some machines, Vina may install cleanly from conda. On others, especially when solving older scientific environments or mixed architecture systems, it may be easier to download the official precompiled Vina executable and place it directly into the active conda environment.

---

<a id="manual-vina-binary-install"></a>

## Manual Vina Binary Install

Use this section if:

* `conda env create -f docking.yaml` fails because of `vina`
* `conda install -c conda-forge vina` fails
* `vina --version` does not work after environment creation
* a downloaded docking package fails with `vina: command not found`

The goal is to place the Vina executable directly into:

```text
$CONDA_PREFIX/bin/vina
```

That makes `vina` available whenever the `docking` environment is active.

### Step 1 — Activate the environment

```bash
conda activate docking
echo "$CONDA_PREFIX"
```

### Step 2 — Install the binary for your platform

#### macOS Apple Silicon / M1 / M2 / M3 / M4

```bash
conda activate docking

curl -L -o "$CONDA_PREFIX/bin/vina" \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64"

chmod +x "$CONDA_PREFIX/bin/vina"

which vina
vina --version
```

#### macOS Intel

```bash
conda activate docking

curl -L -o "$CONDA_PREFIX/bin/vina" \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_x86_64"

chmod +x "$CONDA_PREFIX/bin/vina"

which vina
vina --version
```

#### Linux x86_64

```bash
conda activate docking

curl -L -o "$CONDA_PREFIX/bin/vina" \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64"

chmod +x "$CONDA_PREFIX/bin/vina"

which vina
vina --version
```

Expected result:

```text
AutoDock Vina v1.2.7
```

### If the conda-installed `vina` exists but is broken

You can replace it:

```bash
conda activate docking

rm -f "$CONDA_PREFIX/bin/vina"

curl -L -o "$CONDA_PREFIX/bin/vina" \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64"

chmod +x "$CONDA_PREFIX/bin/vina"

which vina
vina --version
```

Use the matching download URL for your platform.

### Optional: install Vina through conda first

This may work on many systems:

```bash
conda activate docking
conda install -c conda-forge vina

which vina
vina --version
```

If that fails, use the manual binary install above.

---

<a id="running-locally"></a>

## Running Locally

### 1. Clone the repository

```bash
git clone https://github.com/Joey305/autodock-WEBSERVER.git
cd autodock-WEBSERVER
```

### 2. Create the conda environment

```bash
conda env create -f docking.yaml
conda activate docking
```

### 3. Verify required executables

```bash
which python
python --version

which obabel
obabel -V

which vina
vina --version
```

If `vina` is missing, follow [Manual Vina Binary Install](#manual-vina-binary-install).

### 4. Start the development server

```bash
python app.py
```

Then open:

```text
http://127.0.0.1:5050
```

Alternative Flask CLI:

```bash
flask --app app:create_app run --debug --port 5050
```

---

## Environment Variables

| Variable               | Purpose                                                       |
| ---------------------- | ------------------------------------------------------------- |
| `PORTAL_SECRET`        | Flask secret key                                              |
| `PORTAL_DB`            | SQLAlchemy database URI                                       |
| `PORTAL_TMP`           | Workspace root; defaults to `/tmp/autodock_prep`              |
| `PORTAL_MAX_RECEPTORS` | Maximum receptors allowed in one workspace                    |
| `PORTAL_ENV_LINE`      | Optional environment line inserted into generated LSF scripts |
| `PORTAL_PUBLIC_EMAIL`  | Placeholder identity for public mode                          |
| `PUBLIC_MODE`          | Enables public mode when true                                 |
| `ENABLE_AUTH`          | Enables login/auth behavior when true                         |
| `ENABLE_LSF_PACKAGE`   | Allows LSF package generation                                 |
| `DEFAULT_PACKAGE_MODE` | Default package mode: `portable` or `lsf`                     |
| `BABEL_LIBDIR`         | Open Babel plugin directory, important for some deployments   |

Example local environment:

```bash
export PORTAL_SECRET="change-me"
export PUBLIC_MODE=true
export ENABLE_AUTH=false
export PORTAL_TMP="/tmp/autodock_prep"
```

---

<a id="headless-api-tutorial"></a>

## Headless API Tutorial

This tutorial shows the full API workflow using:

* receptor: PDB `3EKY`, chain `A`
* ligand: a user-provided SDF file
* package mode: `portable`
* output: downloadable ZIP package

### Before starting

You need:

```bash
curl
```

and a ligand file in your current directory, for example:

```text
my_ligand.sdf
```

Set variables:

```bash
BASE="https://autodockvina.com"
JOB="demo-docking-$(date +%s)"
LIGAND="my_ligand.sdf"
```

Confirm your ligand exists locally:

```bash
ls -lh "$LIGAND"
```

---

### Step 1 — Health check

```bash
curl -sS "$BASE/api/v1/health"
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "service": "autodock-vina-prepserver",
    "api_version": "v1"
  }
}
```

---

### Step 2 — Create a workspace

```bash
curl -sS -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d "{\"workspace_name\":\"$JOB\"}"
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "jobname": "demo-docking-...",
    "workspace": "/tmp/autodock_prep/demo-docking-..."
  }
}
```

---

### Step 3 — Fetch a receptor by PDB ID

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/receptors/fetch" \
  -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":"A"}'
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "rel": "Receptors/3eky.pdb",
    "receptors": [
      {
        "display": "3eky.pdb",
        "rel": "Receptors/3eky.pdb",
        "status": "new"
      }
    ]
  }
}
```

Important: if the receptor was fetched as `Receptors/3eky.pdb`, use that same receptor path when saving the center.

---

### Step 4 — Save a docking center

This tutorial uses explicit XYZ coordinates.

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/centers/save" \
  -H "Content-Type: application/json" \
  -d '{
    "method":"xyz",
    "receptor":"Receptors/3eky.pdb",
    "center":[10.5, 22.1, -3.4],
    "size":20
  }'
echo ""
```

Check saved centers:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/centers"
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "centers": [
      {
        "receptor_pdbqt": "3eky.pdbqt",
        "center": [10.5, 22.1, -3.4],
        "size": 20.0
      }
    ]
  }
}
```

The docking center must be saved for the receptor PDBQT name expected by the preparation step.

---

### Step 5 — Prepare the receptor

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" \
  -H "Content-Type: application/json" \
  -d '{
    "remove_hets":"all",
    "remove_chains":"",
    "altloc":"collapse"
  }'
echo ""
```

Check preparation status:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/prep/status"
echo ""
```

Check workspace summary:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/summary"
echo ""
```

Success indicators:

```json
{
  "converted_count": 1,
  "converted_list": ["3eky.pdbqt"],
  "prep_done": true,
  "receptors_prepped": 1
}
```

---

### Step 6 — Upload a ligand

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@$LIGAND"
echo ""
```

Check ligand state:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/ligands"
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "ligands": [
      "my_ligand.sdf"
    ],
    "ligand_info": {
      "accepted_count": 1,
      "filetypes": [".sdf"],
      "upload_mode": "single"
    }
  }
}
```

---

### Step 7 — Build the package

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"portable",
    "vina_poses":20,
    "confgen_poses":64
  }'
echo ""
```

Expected response includes:

```json
{
  "ok": true,
  "data": {
    "package_mode": "portable",
    "zip": "/tmp/autodock_prep/.../job.zip",
    "download_url": "/download?path=/tmp/autodock_prep/.../job.zip"
  }
}
```

---

### Step 8 — List and download artifacts

List generated artifacts:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/artifacts"
echo ""
```

Download the latest ZIP into your current local directory:

```bash
curl -L -o "${JOB}_docking_package.zip" \
  "$BASE/api/v1/workspaces/$JOB/download"

ls -lh "${JOB}_docking_package.zip"
```

After this step, the ZIP should appear locally.

Important: workspace files are created on the server under `/tmp/autodock_prep/...`. They will not appear in your local folder until you download the ZIP.

---

## Full API Tutorial Script

Save the following as `run_autodock_api_tutorial.sh`:

```bash
#!/usr/bin/env bash
set -e

BASE="https://autodockvina.com"
JOB="demo-docking-$(date +%s)"
LIGAND="${1:-my_ligand.sdf}"

echo "Using BASE=$BASE"
echo "Using JOB=$JOB"
echo "Using LIGAND=$LIGAND"
echo ""

echo "0) Confirm ligand exists locally"
ls -lh "$LIGAND"
echo ""

echo "1) Health check"
curl -sS "$BASE/api/v1/health"
echo ""
echo ""

echo "2) Create workspace"
curl -sS -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d "{\"workspace_name\":\"$JOB\"}"
echo ""
echo ""

echo "3) Fetch receptor"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/receptors/fetch" \
  -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":"A"}'
echo ""
echo ""

echo "4) Save center"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/centers/save" \
  -H "Content-Type: application/json" \
  -d '{
    "method":"xyz",
    "receptor":"Receptors/3eky.pdb",
    "center":[10.5, 22.1, -3.4],
    "size":20
  }'
echo ""
echo ""

echo "5) Check centers"
curl -sS "$BASE/api/v1/workspaces/$JOB/centers"
echo ""
echo ""

echo "6) Start receptor preparation"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" \
  -H "Content-Type: application/json" \
  -d '{
    "remove_hets":"all",
    "remove_chains":"",
    "altloc":"collapse"
  }'
echo ""
echo ""

echo "7) Check preparation status"
curl -sS "$BASE/api/v1/workspaces/$JOB/prep/status"
echo ""
echo ""

echo "8) Upload ligand"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@$LIGAND"
echo ""
echo ""

echo "9) Check ligands"
curl -sS "$BASE/api/v1/workspaces/$JOB/ligands"
echo ""
echo ""

echo "10) Build package"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"portable",
    "vina_poses":20,
    "confgen_poses":64
  }'
echo ""
echo ""

echo "11) List artifacts"
curl -sS "$BASE/api/v1/workspaces/$JOB/artifacts"
echo ""
echo ""

echo "12) Download ZIP"
curl -L -o "${JOB}_docking_package.zip" \
  "$BASE/api/v1/workspaces/$JOB/download"
echo ""
echo ""

echo "13) Confirm local ZIP"
ls -lh "${JOB}_docking_package.zip"
echo ""

echo "Done."
echo "Downloaded package: ${JOB}_docking_package.zip"
```

Run it:

```bash
chmod +x run_autodock_api_tutorial.sh
./run_autodock_api_tutorial.sh my_ligand.sdf
```

---

<a id="real-example-walkthrough"></a>

## Real Example Walkthrough

This repository includes an example ligand:

```text
Example/Ligand_SDF/DR7.sdf
```

The walkthrough below uses the same staged API flow, but points to the repository example ligand.

### Option A — Run from the repository root

```bash
BASE="https://autodockvina.com"
JOB="dr7-example-$(date +%s)"
LIGAND="Example/Ligand_SDF/DR7.sdf"

echo "Using JOB=$JOB"
echo "Using LIGAND=$LIGAND"
ls -lh "$LIGAND"
```

### Option B — Download the ligand directly from GitHub

Use this if you are not inside the repository:

```bash
mkdir -p PrepServer_DR7_API_Demo
cd PrepServer_DR7_API_Demo

curl -L -o DR7.sdf \
  https://raw.githubusercontent.com/Joey305/autodock-WEBSERVER/main/Example/Ligand_SDF/DR7.sdf

BASE="https://autodockvina.com"
JOB="dr7-example-$(date +%s)"
LIGAND="DR7.sdf"

ls -lh "$LIGAND"
```

### Create the workspace

```bash
curl -sS -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d "{\"workspace_name\":\"$JOB\"}"
echo ""
```

### Fetch receptor `3EKY`, chain `A`

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/receptors/fetch" \
  -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":"A"}'
echo ""
```

### Save an example XYZ docking box

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/centers/save" \
  -H "Content-Type: application/json" \
  -d '{
    "method":"xyz",
    "receptor":"Receptors/3eky.pdb",
    "center":[10.5, 22.1, -3.4],
    "size":20
  }'
echo ""
```

### Prepare receptor

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" \
  -H "Content-Type: application/json" \
  -d '{
    "remove_hets":"all",
    "remove_chains":"",
    "altloc":"collapse"
  }'
echo ""
```

### Confirm receptor preparation

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/summary"
echo ""
```

Success should include:

```json
{
  "prep_done": true,
  "converted_count": 1,
  "converted_list": ["3eky.pdbqt"]
}
```

### Upload `DR7.sdf`

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@$LIGAND"
echo ""
```

Confirm ligand upload:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/ligands"
echo ""
```

### Build a portable package

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"portable",
    "vina_poses":20,
    "confgen_poses":64
  }'
echo ""
```

### Download the generated package

```bash
curl -L -o "${JOB}_docking_package.zip" \
  "$BASE/api/v1/workspaces/$JOB/download"

ls -lh "${JOB}_docking_package.zip"
```

You should now have a local file similar to:

```text
dr7-example-1781216080_docking_package.zip
```

---

## One-Shot DR7 Example Script

Save this as `run_dr7_example.sh`:

```bash
#!/usr/bin/env bash
set -e

BASE="https://autodockvina.com"
JOB="dr7-example-$(date +%s)"
LIGAND="${1:-Example/Ligand_SDF/DR7.sdf}"
DR7_URL="https://raw.githubusercontent.com/Joey305/autodock-WEBSERVER/main/Example/Ligand_SDF/DR7.sdf"

echo "======================================"
echo " AutoDock-Vina PrepServer DR7 Example"
echo "======================================"
echo ""
echo "BASE:   $BASE"
echo "JOB:    $JOB"
echo "LIGAND: $LIGAND"
echo ""

if [ ! -f "$LIGAND" ]; then
  echo "Ligand not found locally. Downloading DR7.sdf..."
  LIGAND="DR7.sdf"
  curl -L -o "$LIGAND" "$DR7_URL"
fi

echo "0) Confirm local example ligand"
ls -lh "$LIGAND"
echo ""

echo "1) Create workspace"
curl -sS -X POST "$BASE/api/v1/workspaces" \
  -H "Content-Type: application/json" \
  -d "{\"workspace_name\":\"$JOB\"}"
echo ""
echo ""

echo "2) Fetch receptor 3EKY chain A"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/receptors/fetch" \
  -H "Content-Type: application/json" \
  -d '{"pdb_id":"3EKY","chains":"A"}'
echo ""
echo ""

echo "3) Save example docking center"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/centers/save" \
  -H "Content-Type: application/json" \
  -d '{
    "method":"xyz",
    "receptor":"Receptors/3eky.pdb",
    "center":[10.5, 22.1, -3.4],
    "size":20
  }'
echo ""
echo ""

echo "4) Prepare receptor"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/prep/start" \
  -H "Content-Type: application/json" \
  -d '{
    "remove_hets":"all",
    "remove_chains":"",
    "altloc":"collapse"
  }'
echo ""
echo ""

echo "5) Check receptor preparation"
curl -sS "$BASE/api/v1/workspaces/$JOB/prep/status"
echo ""
echo ""

echo "6) Upload DR7 ligand"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@$LIGAND"
echo ""
echo ""

echo "7) Build portable docking package"
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"portable",
    "vina_poses":20,
    "confgen_poses":64
  }'
echo ""
echo ""

echo "8) List artifacts"
curl -sS "$BASE/api/v1/workspaces/$JOB/artifacts"
echo ""
echo ""

echo "9) Download package"
curl -L -o "${JOB}_docking_package.zip" \
  "$BASE/api/v1/workspaces/$JOB/download"
echo ""
echo ""

echo "10) Confirm local ZIP"
ls -lh "${JOB}_docking_package.zip"
echo ""

echo "Done."
echo "Downloaded package: ${JOB}_docking_package.zip"
```

Run it:

```bash
chmod +x run_dr7_example.sh
./run_dr7_example.sh
```

---

## Optional SMILES / SMI Ligand Example

If you do not have an SDF file, you can create a simple `.smi` ligand file and upload it through the same ligand endpoint.

```bash
cat > example_ligand.smi <<'EOF'
CC(=O)Oc1ccccc1C(=O)O aspirin_example
EOF

curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@example_ligand.smi"
echo ""

curl -sS "$BASE/api/v1/workspaces/$JOB/ligands"
echo ""
```

---

## What is inside the downloaded ZIP?

A generated portable package typically contains:

```text
job/
├── Receptors/
├── Receptors_PDBQT/
├── Ligands/
├── vina_centers.csv
├── run_confgen_local.sh
├── run_vina_local.sh
├── run_all_local.sh
├── README_RUN_LOCAL.md
└── bundled Python helper scripts
```

After download:

```bash
unzip dr7-example-*_docking_package.zip -d DR7_Docking_Package
cd DR7_Docking_Package
find . -maxdepth 2 -type f | sort | head -50
```

Read the package-specific instructions:

```bash
find . -name "README_RUN_LOCAL.md" -print
```

---

## Running a Downloaded Package

After downloading and unzipping a package from AutoDock-Vina PrepServer:

```bash
unzip my_docking_package.zip -d my_docking_package
cd my_docking_package

conda activate docking
bash run_all_local.sh
```

If ligand preparation is already complete and you only want to run docking:

```bash
conda activate docking
bash run_vina_local.sh
```

If this fails with `vina: command not found`, return to [Manual Vina Binary Install](#manual-vina-binary-install).

---

## Endpoint Reference

| Method | Endpoint                                        | Purpose                               |
| ------ | ----------------------------------------------- | ------------------------------------- |
| `GET`  | `/api/v1/health`                                | Check service health and API version  |
| `POST` | `/api/v1/workspaces`                            | Create a new workspace                |
| `GET`  | `/api/v1/workspaces/<jobname>`                  | Read workspace state                  |
| `GET`  | `/api/v1/workspaces/<jobname>/summary`          | Read workflow summary                 |
| `POST` | `/api/v1/workspaces/<jobname>/receptors/upload` | Upload receptor files                 |
| `POST` | `/api/v1/workspaces/<jobname>/receptors/fetch`  | Fetch receptor by PDB ID              |
| `GET`  | `/api/v1/workspaces/<jobname>/receptors`        | List workspace receptors              |
| `POST` | `/api/v1/workspaces/<jobname>/centers/resolve`  | Resolve a docking center              |
| `POST` | `/api/v1/workspaces/<jobname>/centers/save`     | Save receptor docking center          |
| `GET`  | `/api/v1/workspaces/<jobname>/centers`          | List saved centers                    |
| `POST` | `/api/v1/workspaces/<jobname>/prep/start`       | Prepare receptor PDBQT files          |
| `GET`  | `/api/v1/workspaces/<jobname>/prep/status`      | Read receptor preparation status      |
| `POST` | `/api/v1/workspaces/<jobname>/ligands/upload`   | Upload ligand files                   |
| `GET`  | `/api/v1/workspaces/<jobname>/ligands`          | List uploaded ligands                 |
| `POST` | `/api/v1/workspaces/<jobname>/build`            | Build downloadable package            |
| `GET`  | `/api/v1/workspaces/<jobname>/artifacts`        | List generated ZIP artifacts          |
| `GET`  | `/api/v1/workspaces/<jobname>/download`         | Download latest or requested artifact |

---

## Accepted Input Types

### Receptors

Supported receptor inputs include:

* `.pdb`
* `.pdbqt`
* `.cif`
* `.mmcif`
* `.ent`
* ZIP archives containing receptor files
* folder-style browser uploads
* PDB fetch by RCSB PDB ID

### Ligands

Supported ligand inputs include:

* `.sdf`
* `.smi`
* `.smiles`
* `.csv`
* ZIP archives containing supported ligand files
* folder-style browser uploads

The upload workflow ignores unsupported files and common macOS metadata such as:

```text
.DS_Store
._*
__MACOSX/
```

---

## Generated Package Modes

### Portable mode

Portable mode is the default choice for most users.

It generates a downloadable job archive with local runner helpers such as:

```text
run_confgen_local.sh
run_vina_local.sh
run_all_local.sh
README_RUN_LOCAL.md
```

Use this mode when the downloaded package will be run on a local workstation, server, or non-LSF environment.

Example build request:

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"portable",
    "vina_poses":20,
    "confgen_poses":64
  }'
```

### LSF mode

LSF mode is intended for lab environments that use an LSF scheduler.

Example build request:

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/build" \
  -H "Content-Type: application/json" \
  -d '{
    "package_mode":"lsf",
    "workers":16,
    "queue":"general",
    "project":"brd",
    "mem_per_core":2000,
    "vina_poses":20,
    "confgen_poses":64
  }'
```

Only use LSF mode in environments where the scheduler assumptions match your compute system.

---

<a id="deployment-notes"></a>

## Deployment Notes

AutoDock-Vina PrepServer can be deployed as a standard Flask application in environments that provide:

* Python
* writable temporary storage
* Open Babel for receptor conversion
* optional AutoDock Vina for local execution workflows
* a web dyno/process runner such as Gunicorn

Production-style deployments should:

* configure a strong `PORTAL_SECRET`
* set `PORTAL_TMP` to a writable path
* avoid storing confidential data on public deployments
* confirm `obabel` works before enabling receptor preparation
* gate public access if anonymous uploads are not intended
* clean old `/tmp` workspaces periodically if persistent storage is used

---

<a id="heroku--open-babel-notes"></a>

## Heroku + Open Babel Notes

For Heroku deployments, Open Babel is installed as a system package using an `Aptfile`.

### Aptfile

```text
openbabel
```

### Buildpacks

The apt buildpack should run before the Python buildpack:

```bash
heroku buildpacks:add --index 1 heroku-community/apt --app autodockvina
heroku buildpacks:add --index 2 heroku/python --app autodockvina
```

### Verify Open Babel

```bash
heroku run --no-tty "which obabel && obabel -V" --app autodockvina
```

If Open Babel is installed but cannot find plugins, locate the plugin directory:

```bash
heroku run --no-tty "find /app/.apt -type f \( -name 'pdbformat.so' -o -name 'pdbqtformat.so' -o -name 'plugindefines.txt' \) -print" --app autodockvina
```

Then set `BABEL_LIBDIR`, for example:

```bash
heroku config:set BABEL_LIBDIR=/app/.apt/usr/lib/x86_64-linux-gnu/openbabel/3.1.1 --app autodockvina
heroku restart --app autodockvina
```

Verify format loading:

```bash
heroku run --no-tty 'echo BABEL_LIBDIR=$BABEL_LIBDIR && which obabel && obabel -V && obabel -L formats | grep -E "pdb|pdbqt"' --app autodockvina
```

Expected output includes:

```text
pdb -- Protein Data Bank format
pdbqt -- AutoDock PDBQT format
```

---

## Testing and Validation

Run Python syntax checks:

```bash
python3 -m py_compile app.py packager.py runner_templates.py manage.py
```

Run tests:

```bash
python3 -m unittest discover -s tests -v
```

Run the API example:

```bash
chmod +x run_dr7_example.sh
./run_dr7_example.sh
```

Manual success checklist:

* `/api/v1/health` returns `"ok": true`
* workspace creation returns a valid `jobname`
* receptor fetch returns `Receptors/3eky.pdb`
* centers list contains `3eky.pdbqt`
* prep summary shows `prep_done: true`
* ligand upload shows `accepted_count: 1`
* build returns a `job.zip`
* download creates a local ZIP file
* `conda activate docking && vina --version` works before running the downloaded package

---

<a id="troubleshooting"></a>

## Troubleshooting

### `workspace_missing`

Example:

```json
{
  "ok": false,
  "error": "workspace_missing"
}
```

Possible causes:

* the workspace name is wrong
* the server restarted and `/tmp` was cleared
* the workspace was created on a different deployment or dyno state

Fix:

```bash
JOB="demo-docking-$(date +%s)"
```

Create a fresh workspace and rerun the staged workflow.

---

### `centers_incomplete`

Example:

```json
{
  "ok": false,
  "error": "centers_incomplete",
  "message": "Save a center for each receptor first."
}
```

Cause:

The receptor exists, but no docking center was saved for its expected PDBQT name.

If the receptor is:

```text
Receptors/3eky.pdb
```

then the center must map to:

```text
3eky.pdbqt
```

Check centers:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/centers"
```

---

### `obabel_missing`

Example:

```json
{
  "ok": false,
  "error": "obabel_missing"
}
```

Cause:

Open Babel is not installed or not on `PATH`.

Check:

```bash
which obabel
obabel -V
```

On Heroku:

```bash
heroku run --no-tty "which obabel && obabel -V" --app autodockvina
```

For local conda installs:

```bash
conda activate docking
conda install -c conda-forge openbabel
obabel -V
```

---

### Open Babel plugin error

Example log:

```text
Unable to find OpenBabel plugins. Try setting the BABEL_LIBDIR environment variable.
```

Fix on Heroku:

```bash
heroku run --no-tty "find /app/.apt -type f \( -name 'pdbformat.so' -o -name 'pdbqtformat.so' -o -name 'plugindefines.txt' \) -print" --app autodockvina
heroku config:set BABEL_LIBDIR=/app/.apt/usr/lib/x86_64-linux-gnu/openbabel/3.1.1 --app autodockvina
heroku restart --app autodockvina
```

---

### `receptor_conversion_failed`

Example:

```json
{
  "ok": false,
  "error": "receptor_conversion_failed"
}
```

Check prep status:

```bash
curl -sS "$BASE/api/v1/workspaces/$JOB/prep/status"
```

Look at the returned `log` text. Common causes include:

* Open Babel plugin path not configured
* malformed receptor input
* unsupported conversion behavior
* invalid or empty cleaned receptor file

---

### `vina: command not found`

If a generated docking package fails with:

```text
vina: command not found
```

activate the environment:

```bash
conda activate docking
```

Then check:

```bash
which vina
vina --version
```

If `which vina` returns nothing, install Vina manually:

```bash
curl -L -o "$CONDA_PREFIX/bin/vina" \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64"

chmod +x "$CONDA_PREFIX/bin/vina"

which vina
vina --version
```

Use the matching binary URL for your operating system.

---

### `Permission denied` when running Vina

If Vina exists but will not execute:

```bash
chmod +x "$CONDA_PREFIX/bin/vina"
vina --version
```

---

### Local `ls` does not show server files

This is expected.

The workspace exists on the server, for example:

```text
/tmp/autodock_prep/demo-docking-...
```

Your local folder will only contain output after you download the generated ZIP:

```bash
curl -L -o "${JOB}_docking_package.zip" \
  "$BASE/api/v1/workspaces/$JOB/download"
```

---

### `ligands_missing`

The package cannot be built until a ligand is uploaded.

Upload a ligand:

```bash
curl -sS -X POST "$BASE/api/v1/workspaces/$JOB/ligands/upload" \
  -F "mode=single" \
  -F "file=@my_ligand.sdf"
```

---

## Security and Data Notes

Do not upload confidential receptor structures, proprietary ligand libraries, or sensitive project data to a public deployment unless the instance is explicitly configured for that purpose.

Public deployments are useful for demonstrations and educational workflows, but private or proprietary research should use a controlled deployment.

Recommended practices:

* use private deployments for sensitive structures
* clean old workspaces regularly
* avoid committing generated ZIPs containing confidential molecules
* review downloaded packages before sharing
* use strong secrets and authentication for non-public instances

---

## Citation and Acknowledgements

AutoDock-Vina PrepServer builds on widely used open-source scientific software and molecular modeling tools.

Acknowledgements:

* Flask and related extensions for the web framework
* Open Babel for chemical file conversion
* AutoDock Vina for docking workflows
* RDKit for cheminformatics utilities
* AutoDockTools-derived utilities bundled in this repository
* Schürer Lab molecular design workflows and tooling ecosystem

If you use this tool in academic work, please cite the relevant upstream software packages used in your workflow, including AutoDock Vina and Open Babel.

---

## Contact

For questions, deployment feedback, or collaboration:

```text
jmschulz@med.miami.edu
```

Website:

```text
https://autodockvina.com
```

Repository:

```text
https://github.com/Joey305/autodock-WEBSERVER
```

Schürer Lab:

```text
https://schurerlab.org
```
