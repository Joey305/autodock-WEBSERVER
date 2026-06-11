# AutoDock-Vina PrepServer

AutoDock-Vina PrepServer is a Flask-based web application for assembling molecular docking job packages. It helps a user upload or fetch receptor structures, define docking boxes, prepare receptors, upload ligand inputs, and export a ready-to-run package for either local execution or an LSF-style HPC workflow.

The repository also includes the Python scripts that the generated package uses for ligand conformer generation, AutoDock Vina batch docking, score parsing, PyMOL export, and plotting.

## Who this is for

- Researchers preparing AutoDock Vina screening jobs
- Lab members who want a browser-based setup workflow instead of hand-building folders
- Developers who need to maintain or extend the portal and its packaged runtime scripts

## Features

- Public, no-login workflow by default
- Receptor intake from single files, ZIP uploads, folder uploads, or RCSB PDB fetches
- Interactive docking-center selection in the browser
- Receptor cleanup and preparation into PDBQT-compatible outputs
- Ligand intake from `.sdf`, `.smiles`, `.smi`, `.csv`, ZIP archives, or folder uploads
- Ligand upload normalization that ignores macOS junk files and flattens nested ligand paths
- Build output in two modes:
  - `portable` for generic local execution
  - `lsf` for lab/HPC packaging with generated scheduler scripts
- Included post-processing scripts for docking score parsing, aggregation, plotting, and PyMOL-oriented export

## Tech Stack

- Backend: [Flask](https://flask.palletsprojects.com/), Flask-Login, Flask-WTF, Flask-SQLAlchemy
- Database: SQLite by default
- Language: Python 3
- Scientific Python libraries used in this repo: RDKit, pandas, NumPy, matplotlib
- External scientific executables used by the workflow:
  - Open Babel (`obabel`) for receptor and ligand conversion steps
  - AutoDock Vina (`vina`) for docking runs in generated packages
- Bundled tool source: `AutoDockTools_py3/`
- Frontend: server-rendered HTML templates plus static CSS/JS assets; no Node build pipeline is present
- Scheduler support: LSF helper script generation

## Repository Structure

```text
.
├── app.py                         # Flask app entrypoint
├── manage.py                      # Simple user-management CLI for auth-enabled mode
├── packager.py                    # Workspace assembly and ZIP packaging
├── runner_templates.py            # Portable runner script templates
├── requirements.txt               # Python dependencies used by the app/tests
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

- Python 3.9+ recommended
- `pip`

### Required for the full docking workflow

- Open Babel available as `obabel`
- AutoDock Vina available as `vina` if you plan to execute the generated docking package

### Optional or environment-specific

- RDKit: required by ligand conformer generation and some downstream utilities
- `gemmi` and optionally `meeko`: used by the standalone receptor preparation script `2a_PDB2PDBQTbatch.py`
- An LSF environment: only needed if you will generate and run `lsf` packages

If you are unsure which Python version you have, check with:

```bash
python3 --version
```

To confirm external executables:

```bash
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

Notes:

- `requirements.txt` covers the Flask app and test dependencies currently tracked in the repo.
- External tools such as `obabel` and `vina` are not installed by `pip install -r requirements.txt`.
- If you plan to run `2a_PDB2PDBQTbatch.py` directly, you may also need to install `gemmi` and optionally `meeko`.

## Environment Configuration

This project includes a sample environment file at `.env.example`.

### 1. Copy the template

```bash
cp .env.example .env
```

### 2. Edit the values

Use safe, local values only. Do not put production or shared secrets into version control.

### 3. Load the variables into your shell

The app does not currently auto-load `.env`, so load it before starting the server:

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
  Optional. Default `true`. Keeps the portal open without user login.
- `ENABLE_AUTH`
  Optional. Default `false`. Enables legacy auth-dependent behavior if set to `true`.
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

This starts the development server on:

```text
http://127.0.0.1:5050
```

### Option 2: Use Flask CLI

```bash
flask --app app:create_app run --debug --port 5050
```

### What to expect

- The homepage is public by default.
- Creating a workspace triggers the upload/build workflow in the browser.
- Uploaded workspace data is stored under `PORTAL_TMP`.
- The local SQLite database is created wherever your configured `PORTAL_DB` points.

## Docker

No `Dockerfile`, `docker-compose.yml`, or other container configuration is currently present in this repository.

If you want Docker-based setup before publishing broadly, add and validate container files separately rather than assuming a container workflow already exists.

## Typical Usage Workflow

1. Open the app in the browser.
2. Create a new workspace.
3. Add receptor structures by uploading files, uploading a ZIP/folder, or fetching a PDB entry.
4. Open each receptor in the viewer and save a docking center.
5. Prepare receptors so PDBQT-ready files are generated.
6. Upload ligand input as a single file, ZIP, or folder.
7. If the ligand input is a CSV, map the SMILES column and optional ligand ID column.
8. Choose a package mode:
   - `portable` for local or generic environments
   - `lsf` for scheduler-oriented lab packaging
9. Build and download the generated ZIP archive.
10. Move the ZIP to the target execution environment and run the packaged scripts there.

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

Portable builds include the normalized runtime files plus helper scripts such as:

- `run_confgen_local.sh`
- `run_vina_local.sh`
- `run_all_local.sh`
- `README_RUN_LOCAL.md`

Use this when you want a generic ZIP for local or non-LSF environments.

### LSF mode

LSF builds include the portable content plus scheduler-oriented files such as:

- `1B_confgen_batch.py`
- `3B_ServerDocks.py`
- `4B_LSFbatch.py`
- `runDOCKING-tmux.sh`
- generated `.lsf` submission scripts

Use this only if your lab or cluster environment already supports the expected LSF workflow.

## Scripts and Commands

### App and support commands

- `python3 app.py`
  Starts the Flask development server on port `5050`.
- `flask --app app:create_app run --debug --port 5050`
  Starts the app via Flask CLI.
- `python3 manage.py`
  Opens a simple interactive CLI for user management if auth mode is enabled.

### Workflow scripts included in this repo

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

This repository currently uses `unittest`-style tests under `tests/`.

### Run the test suite

```bash
python3 -m unittest discover -s tests -v
```

### Run a basic syntax check

```bash
python3 -m py_compile app.py packager.py runner_templates.py manage.py
```

### Manual validation checklist

If you do not have all scientific executables installed locally, at least verify the web app flow manually:

1. Start the app.
2. Open `http://127.0.0.1:5050`.
3. Create a workspace.
4. Upload a small receptor file.
5. Save a center.
6. Upload a small ligand file.
7. Build a `portable` package and confirm a ZIP is produced.

## Deployment Notes

This repository does not include a production deployment manifest, reverse-proxy config, or container image definition.

If you deploy it yourself:

- set a strong `PORTAL_SECRET`
- use a writable, non-repo path for `PORTAL_TMP`
- use a writable database location
- ensure `obabel` is installed if users will perform receptor preparation
- disable or gate public mode if the deployment should not be open to everyone

## GitHub Publishing Instructions

Review the repository carefully before the first public push, especially:

- `.env` and any other local environment files
- generated ZIPs and docking result folders
- local databases and logs
- any secrets that may have been committed previously

Important: `.gitignore` only affects new untracked files. This repository currently has tracked runtime artifacts such as `.pyc`, `.db`, and log files, so review and remove them from version control before the first public release commit if they are not meant to ship:

```bash
git rm --cached -r __pycache__
git rm --cached instance/*.db logs/server.out
git rm --cached .DS_Store AutoDockTools_py3/.DS_Store AutoDockTools_py3/AutoDockTools/.DS_Store
```

After review, the owner can publish manually with:

```bash
git init
git status
git add .
git commit -m "Initial public release"
git remote add origin <repo-url>
git branch -M main
git push -u origin main
```

If this repository is already initialized locally, skip `git init` and keep the existing history as appropriate.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for a lightweight contribution workflow.

At a minimum:

1. Fork the repository.
2. Create a branch.
3. Make focused changes.
4. Run validation commands.
5. Open a pull request with a clear summary.

## Troubleshooting

### `source .env` works but the app still uses defaults

Confirm you loaded the file into the current shell:

```bash
set -a
source .env
set +a
env | rg '^PORTAL_|^PUBLIC_MODE|^ENABLE_'
```

### `obabel` is missing

Receptor preparation and several packaged workflow steps rely on Open Babel. Install it and confirm:

```bash
obabel -V
```

### `vina` is missing

The generated docking scripts expect a working AutoDock Vina executable. Confirm:

```bash
vina --version
```

If needed, set:

```bash
export VINA_EXE=/path/to/vina
```

### RDKit installation is difficult on your machine

RDKit can be easier to install through Conda-based scientific environments than through a plain system Python setup. If your target environment already provides RDKit, use that environment for workflow-script execution.

### Port `5050` is already in use

Start the app with a different port:

```bash
flask --app app:create_app run --debug --port 5051
```

### The build button fails for `lsf` mode

Check whether `ENABLE_LSF_PACKAGE=true` is set in your environment. The UI intentionally falls back to `portable` mode when LSF packaging is disabled.

## Security and Data Notes

- Do not commit `.env`, credentials, or cluster-specific secrets.
- Uploaded inputs and generated workspaces are stored on disk under `PORTAL_TMP`.
- Local SQLite databases, logs, caches, and generated job outputs should not be committed.
- If a secret was ever committed previously, rotate it before making the repository public.
- Review generated docking result folders and ZIP packages before sharing them publicly; they may contain project-specific inputs or derived outputs.

## License

No repository-level license file is currently present.

Before publishing publicly, the repository owner should choose and add a license so users know how they may use, modify, and redistribute the project.

## Acknowledgements

- [Flask](https://flask.palletsprojects.com/) and related extensions used by the portal
- [RDKit](https://www.rdkit.org/) for cheminformatics processing
- [AutoDock Vina](https://vina.scripps.edu/) for docking workflows
- Open Babel for chemical file conversion
- The bundled `AutoDockTools_py3` source included in this repository
