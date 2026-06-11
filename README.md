# AutoDock-Vina PrepServer

AutoDock-Vina PrepServer is a Flask-based web application for assembling molecular docking job packages. It provides a browser workflow for collecting receptor structures, defining docking boxes, preparing receptors, uploading ligand inputs, and exporting ready-to-run AutoDock Vina package archives for local execution or LSF-style cluster workflows.

The repository also includes the Python scripts bundled into generated packages for ligand conformer generation, batch docking, score parsing, result aggregation, plotting, and PyMOL-oriented export.

## Overview

This project is intended for teams that want a reproducible, browser-driven entry point for AutoDock Vina job preparation without hand-building directory trees for every run. It is especially useful when researchers, students, or shared lab users need a lightweight portal for packaging docking studies while still retaining access to the underlying runtime scripts.

## Who This Is For

- Researchers preparing AutoDock Vina screening jobs
- Lab members who want a guided browser workflow for receptor and ligand intake
- Users generating portable ZIP packages for downstream docking execution
- HPC-oriented groups packaging jobs for LSF-based environments
- Developers extending the Flask app or the packaged docking utilities
- Contributors evaluating the project structure and validation workflow

## Features

- Public, browser-based workflow by default
- Receptor intake from single-file uploads, ZIP uploads, folder uploads, or RCSB PDB fetches
- Interactive docking-center selection in the browser
- Versioned headless API for scripted workspace creation, coordinate-based box generation, receptor preparation, ligand upload, and package generation
- Receptor preparation into docking-ready outputs
- Ligand intake from `.sdf`, `.smiles`, `.smi`, `.csv`, ZIP archives, or folder uploads
- Ligand upload normalization that ignores common macOS metadata files and flattens nested ligand paths
- Package generation in two modes:
  - `portable` for generic local execution
  - `lsf` for environments that use LSF scheduler workflows
- Bundled post-processing scripts for score parsing, aggregation, plotting, and PyMOL-oriented export

## Tech Stack

- Backend: [Flask](https://flask.palletsprojects.com/), Flask-Login, Flask-WTF, Flask-SQLAlchemy
- Database: SQLite by default
- Language: Python 3
- Scientific Python libraries used in this repository: RDKit, pandas, NumPy, matplotlib
- External scientific executables used by the workflow:
  - Open Babel (`obabel`) for receptor and ligand conversion steps
  - AutoDock Vina (`vina`) for docking runs in generated packages
- Bundled tool source: `AutoDockTools_py3/`
- Frontend: server-rendered HTML templates with static CSS/JS assets
- Scheduler support: LSF helper script generation

## Repository Structure

```text
.
├── app.py                         # Flask app entrypoint
├── manage.py                      # User-management CLI for auth-enabled deployments
├── packager.py                    # Workspace assembly and ZIP packaging
├── runner_templates.py            # Portable runner script templates
├── requirements.txt               # Python dependencies used by the app and tests
├── templates/                     # Flask HTML templates
├── static/                        # CSS and images
├── tests/                         # unittest-based regression tests
├── AutoDockTools_py3/             # Bundled AutoDockTools source tree
├── 1_ConformerGeneration.py       # Ligand conformer/PDBQT generation
├── 2a_PDB2PDBQTbatch.py           # Receptor preparation utility
├── 3_Complete_batch_docking.py    # AutoDock Vina batch runner
├── 4_ParseScores.py               # Vina score parsing
├── 4C_ConcatenateScores.py        # Combined score collation
├── 5C_BuildPymolSesh.py           # PyMOL/session export workflow
├── 6_MDpymacs.py                  # Additional downstream analysis helper
└── 7_Graphs.py                    # Plotting and ranking summaries
```

## Prerequisites

### Required to run the web app

- Python 3.9 or newer is recommended
- `pip`

### Required for the full docking workflow

- Open Babel available as `obabel`
- AutoDock Vina available as `vina` if generated packages will be executed

### Optional or environment-specific

- RDKit for ligand conformer generation and downstream utilities
- `gemmi` and optionally `meeko` for the standalone receptor preparation script `2a_PDB2PDBQTbatch.py`
- An LSF environment if `lsf` package generation will be used in practice

Check local versions with:

```bash
python3 --version
obabel -V
vina --version
```

## Installation

### 1. Clone the repository

```bash
git clone <repo-url>
cd autodock-WEBSERVER
```

### 2. Create a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install Python dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

`requirements.txt` covers the Flask application and currently tracked test dependencies. External executables such as `obabel` and `vina` are installed separately.

## Environment Configuration

The repository includes a sample environment file at `.env.example`.

### 1. Copy the template

```bash
cp .env.example .env
```

### 2. Edit the values

Use deployment-appropriate values and keep `.env` out of version control.

### 3. Load the variables into your shell

The app does not auto-load `.env`, so load it before starting the server:

```bash
set -a
source .env
set +a
```

### Environment variable reference

- `PORTAL_SECRET`
  Required. Flask secret key for session signing.
- `PORTAL_DB`
  Optional. SQLAlchemy database URL. Defaults to SQLite if not set.
- `PORTAL_TMP`
  Optional. Workspace root for uploaded files and generated job packages.
- `PORTAL_MAX_RECEPTORS`
  Optional. UI limit for receptor count.
- `PORTAL_ENV_LINE`
  Optional. Shell activation line inserted into generated LSF scripts.
- `PORTAL_PUBLIC_EMAIL`
  Optional. Placeholder identity used in public mode.
- `PUBLIC_MODE`
  Optional. Defaults to `true`. Keeps the portal open without login.
- `ENABLE_AUTH`
  Optional. Defaults to `false`. Enables auth-dependent behavior when set to `true`.
- `ENABLE_LSF_PACKAGE`
  Optional. Enables `lsf` package generation in the UI.
- `DEFAULT_PACKAGE_MODE`
  Optional. Default package mode for builds. Supported values are `portable` and `lsf`.

## Running Locally

After installing dependencies and loading environment variables, start the Flask app with one of the following commands.

### Option 1: Run the module directly

```bash
python3 app.py
```

This starts the development server on `http://127.0.0.1:5050`.

### Option 2: Use Flask CLI

```bash
flask --app app:create_app run --debug --port 5050
```

### What to expect

- The homepage is public by default.
- Creating a workspace launches the upload and packaging workflow in the browser.
- Uploaded workspace data is stored under `PORTAL_TMP`.
- The local SQLite database is created wherever `PORTAL_DB` points.

## Typical Usage Workflow

1. Open the app in a browser.
2. Create a new workspace.
3. Add receptor structures by uploading files, uploading a ZIP or folder, or fetching an RCSB PDB entry.
4. Open each receptor in the viewer and save a docking center.
5. Prepare receptors so docking-ready files are generated.
6. Upload ligand input as a single file, ZIP archive, or folder.
7. If the ligand input is a CSV, map the SMILES column and optional ligand ID column.
8. Choose a package mode:
   - `portable` for local or generic environments
   - `lsf` for scheduler-oriented lab packaging
9. Build and download the generated ZIP archive.
10. Move the ZIP to the target execution environment and run the packaged scripts there.

## Headless API

PrepServer includes a versioned headless API under `/api/v1` for scripted workspace creation, receptor upload or fetch, residue/HETATM/XYZ docking-box generation, receptor preparation, ligand upload, and package generation.

Start with:

- `docs/HEADLESS_API.md` for the route table and workflow notes
- `docs/API_EXAMPLES.md` for curl and Python examples
- `docs/AGENT_WORKFLOW.md` for staged command-line or future agent usage
- `docs/examples/` for short `requests` scripts

The API follows the app's current deployment behavior and does not add authentication, rate limiting, quotas, or a production queue system.

## Accepted Input Types

### Receptors

- `.pdb`
- `.cif`
- `.mmcif`
- `.ent`
- `.pdbqt`
- ZIP archives containing receptor files
- Browser folder uploads containing receptor files

### Ligands

- `.sdf`
- `.smiles`
- `.smi`
- `.csv`
- ZIP archives containing supported ligand files
- Browser folder uploads containing supported ligand files

The ligand upload path intentionally ignores unsupported files and common macOS metadata such as `.DS_Store`, `._*`, and `__MACOSX/`.

## Generated Package Modes

### Portable mode

Portable builds include normalized runtime files plus helper scripts such as:

- `run_confgen_local.sh`
- `run_vina_local.sh`
- `run_all_local.sh`
- `README_RUN_LOCAL.md`

Use this mode for local execution or non-LSF environments.

### LSF mode

LSF builds include the portable content plus scheduler-oriented files such as:

- `1B_confgen_batch.py`
- `3B_ServerDocks.py`
- `4B_LSFbatch.py`
- `runDOCKING-tmux.sh`
- generated `.lsf` submission scripts

Use this mode only in environments that already provide the expected LSF scheduler workflow.

## Scripts and Commands

### App and support commands

- `python3 app.py`
  Starts the Flask development server on port `5050`.
- `flask --app app:create_app run --debug --port 5050`
  Starts the app via Flask CLI.
- `python3 manage.py`
  Opens the user-management CLI for auth-enabled deployments.

### Workflow scripts included in this repository

- `python3 1_ConformerGeneration.py --help`
  Ligand conformer generation and PDBQT conversion.
- `python3 2a_PDB2PDBQTbatch.py --help`
  Receptor preparation and PDBQT conversion.
- `python3 3_Complete_batch_docking.py --help`
  Batch AutoDock Vina execution.
- `python3 4_ParseScores.py --help`
  Parse Vina outputs into CSV summaries.
- `python3 5C_BuildPymolSesh.py --help`
  Build corrected export artifacts and PyMOL-oriented outputs.
- `python3 7_Graphs.py --help`
  Create score summaries, rankings, plots, and reports.

## Testing and Validation

This repository uses `unittest`-style tests under `tests/`.

Run the test suite with:

```bash
python3 -m unittest discover -s tests -v
```

Run a basic syntax check with:

```bash
python3 -m py_compile app.py packager.py runner_templates.py manage.py
```

### Manual validation checklist

If all scientific executables are not available locally, a lightweight browser workflow check is still useful:

1. Start the app.
2. Open `http://127.0.0.1:5050`.
3. Create a workspace.
4. Upload a small receptor file.
5. Save a center.
6. Upload a small ligand file.
7. Build a `portable` package and confirm that a ZIP archive is produced.

## Deployment Notes

AutoDock-Vina PrepServer can be deployed as a standard Flask application in environments that provide Python, writable storage, and the scientific executables required by the workflows users will run.

For production-style deployments:

- configure a strong `PORTAL_SECRET`
- set `PORTAL_TMP` to a writable non-repository path
- use a writable database location
- install `obabel` if receptor preparation will be available to users
- install `vina` if generated packages will be executed in the deployment environment
- disable or gate public access if open anonymous use is not intended

This repository does not include a production deployment manifest, reverse-proxy configuration, or container image definition.

## Docker

This repository does not currently include Docker configuration. Use the local Python setup instructions above.

## Troubleshooting

### Environment variables appear unchanged after editing `.env`

Confirm the file was loaded into the current shell:

```bash
set -a
source .env
set +a
env | rg '^PORTAL_|^PUBLIC_MODE|^ENABLE_'
```

### `obabel` is missing

Receptor preparation and several packaged workflow steps rely on Open Babel. Confirm availability with:

```bash
obabel -V
```

### `vina` is missing

The generated docking scripts expect a working AutoDock Vina executable. Confirm availability with:

```bash
vina --version
```

If needed:

```bash
export VINA_EXE=/path/to/vina
```

### RDKit installation is difficult on your machine

RDKit can be easier to install through Conda-based scientific environments than through a plain system Python setup. If the target execution environment already provides RDKit, use that environment for workflow-script execution.

### Port `5050` is already in use

Start the app on a different port:

```bash
flask --app app:create_app run --debug --port 5051
```

### The build button fails for `lsf` mode

Check whether `ENABLE_LSF_PACKAGE=true` is set in the environment. The UI falls back to `portable` mode when LSF packaging is disabled.

## Security and Data Notes

- Do not commit `.env`, credentials, or cluster-specific secrets.
- Use a strong `PORTAL_SECRET`, especially outside local development.
- Uploaded inputs and generated workspaces are stored on disk under `PORTAL_TMP`.
- Use a writable storage path outside source-controlled directories for `PORTAL_TMP` in shared or deployed environments.
- Review generated docking packages before sharing them; they may contain project-specific inputs, intermediate files, or derived outputs.
- Public mode should only be used in environments where open access is intended.

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for setup, validation, and pull request guidance.

## License

This repository does not currently include a license file. Use, modification, and redistribution rights are therefore not granted unless provided separately by the maintainers.

## Acknowledgements

- [Flask](https://flask.palletsprojects.com/) and related extensions used by the portal
- [RDKit](https://www.rdkit.org/) for cheminformatics processing
- [AutoDock Vina](https://vina.scripps.edu/) for docking workflows
- Open Babel for chemical file conversion
- The bundled `AutoDockTools_py3` source included in this repository
