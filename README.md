# AutoDocking Prep Portal

A login‑gated Flask web app that streamlines **end‑to‑end molecular docking job prep** for HPC clusters (LSF) and local experimentation. It focuses on fast, repeatable packaging of everything you need to run AutoDock Vina (or Vina‑compatible tools): receptor prep, box/center selection, ligand prep, and final job ZIPs that include executable submit scripts and templates.


---

## Highlights

* **Three guided modules** in one page, revealed as you go:

  1. **Protein Prep** – upload/fetch receptor(s), optional chain filtering, solvent/ion cleanup.
  2. **Box Generation** – pick centers from a bound ligand or by clicking atoms/residues in a 3D viewer; writes a consolidated `vina_centers*.csv`.
  3. **Ligand Prep** – single file, CSV of SMILES, or a folder ZIP (SDF/SMILES). Autodetects and passes the right flags to your confgen step.
* **Batch‑friendly**: process up to 5 receptors per job; the viewer walks you through each to ensure a center for all.
* **LSF ready**: generates fully parameterized `.lsf` scripts (queue, project, walltime, cores, mem/core, poses, etc.) and a submitter.
* **Portable job bundles**: produces a single ZIP with `Receptors/`, `Ligands/`, centers CSV, and placeholder (or your real) Python pipelines.
* **Viewer‑first UX**: integrated NGL 3D canvas; click atoms or type selections (e.g., `resname LIG`, `chain A and resi 82`).
* **Security**: email/password auth (Flask-Login), CSRF protection.

---

## Quick start

### 1) Install system deps

* Python 3.10+
* Conda or venv (recommended)
* (Optional) RDKit for local cheminformatics
* (Optional) Tailscale if you want VPN access across machines

### 2) Create an environment

```bash
conda create -n docking python=3.10 -y
conda activate docking
pip install -r requirements.txt
```

If you use RDKit:

```bash
conda install -c conda-forge rdkit
```

### 3) Initialize the app

```bash
export PORTAL_SECRET="change-me"
export PORTAL_DB="sqlite:///portal.db"
export PORTAL_TMP="/tmp/autodock_prep"
# Default conda env line used in LSF templates (customizable in UI)
export PORTAL_ENV_LINE="source /nethome/<user>/miniconda3/etc/profile.d/conda.sh && conda activate vina_env"
python app.py  # or: flask run -h 0.0.0.0 -p 5050
```

Visit `http://localhost:5050`. Create a user with the included auth blueprint (first login path typically exposes a basic sign‑up if seeded; otherwise seed via shell/SQLite).

### 4) (Optional) Make it reachable on your tailnet

If running under WSL/Linux:

```bash
sudo tailscale up
# Access from another device on the tailnet at: http://<tailscale-ip>:5050/
```

On Windows host you can also map paths via `tailscale serve` (see Tailscale docs).

---

## Workflow

### A) Protein Prep

* **Upload** a single `.pdb`/`.pdbqt` or **ZIP a Receptors/** folder; or **Fetch by PDB ID** (optional chain filter `A,B`).
* The app lists non‑standard residues (waters/ions/ligands) and lets you apply removal. A cleaned PDB is written alongside the original.

### B) Box / Center Generation

* Use the integrated viewer to:

  * **Ligand‑driven**: type `resname <3‑letter>` to select a bound ligand and auto‑center.
  * **Click‑driven**: click any atoms to build a selection; you can clear/re‑pick until satisfied.
  * **Residue form**: alternatively enter `resname`, `resi`, `chain` and we’ll compute the center of that selection.
* Press **Save center** to append one row per receptor to `vina_centers.csv` (later renamed to include receptor tags).

### C) Ligand Prep

* Provide **one** of:

  * **Single file**: `.sdf` / `.smiles` / `.csv` (with SMILES + optional ID column)
  * **Folder ZIP**: `Ligands/` containing SDF or SMILES files
* The backend autodetects the input mode for the confgen stage.

### D) Build ZIP (LSF + assets)

* Fill out queue, project, walltime, cores, mem/core, number of conformers and docking poses; optionally set a **VINA\_EXE** and **env line**.
* Click **Build Zip**. You’ll get a downloadable job bundle with:

  * `Receptors/`, `Ligands/`
  * `vina_centers*.csv`
  * `LSF/` for ConfGen and Vina
  * Your pipeline scripts (placeholders unless you replace them)

---

## File conventions

* **Receptors**: `.pdb` for visualization and cleanup; `.pdbqt` for docking. If you upload only PDBQT, you can still build jobs but the 3D viewer won’t render.
* **Centers CSV**: `PDB_ID,X,Y,Z,SIZE` — one row per receptor. During packaging it’s renamed to include receptor tags, e.g., `vina_centers_HDGF_… .csv`.
* **Ligands**: SDF/SMILES files; or CSV with SMILES and an optional ID column.

---

## LSF templates

Templates include the usual headers (queue, project, walltime, core count, mem/core, span-on-host, email notify) and a customizable **env block**:

```bash
# Conda env with vina
source /nethome/<user>/miniconda3/etc/profile.d/conda.sh
conda activate vina_env
export VINA_EXE="$HOME/miniconda3/envs/vina_env/bin/vina"  # overridable in UI
```

You can swap the env block per job via the UI. Scripts submit ConfGen first, then Vina with the generated centers file.

---

## Architecture (short)

* **Flask** app with **Flask‑Login** and **SQLAlchemy** (SQLite by default)
* **Templates**: Bootstrap 5, NGL for 3D viewer
* **Endpoints** (selection):

  * `POST /api/workspace` – create job workspace under `$PORTAL_TMP`
  * `POST /api/upload` – upload single/ZIP for receptors/ligands; autodetects content
  * `POST /api/fetch_pdb` – fetch from RCSB, optional `chains=A,B`
  * `GET  /api/list_residues` – list removable residue names
  * `POST /api/clean_receptor` – rewrite cleaned PDB
  * `POST /api/center` – append center row
  * `POST /api/build` – generate LSFs & bundle ZIP
  * `GET  /download?path=...` – download artifacts

---

## Security & operations

* Login required for all build actions.
* Each job lives in a per‑user **workspace** path under `$PORTAL_TMP`.
* CSRF is enabled via WTForms/Flask‑WTF in auth views.
* For multi‑user or production, run behind a proper WSGI server (gunicorn/uwsgi) + reverse proxy (nginx) and use a real database.

---

## Troubleshooting

* **Viewer shows nothing**: you likely loaded only PDBQT. Upload/fetch a PDB to visualize.
* **Cannot fetch PDB**: verify internet egress and PDB ID. Use `chains` like `A,B` (no spaces).
* **Tailscale access**: ensure the service binds `0.0.0.0`, `tailscale ip -4` shows an address, and the port is listening (`ss -ltnp | grep :5050`).
* **RDKit missing**: install from conda‑forge; the app runs without it, but certain ligand conveniences will be limited.

---

## Roadmap

* Optional batch >5 receptors with paging
* Built‑in PDB→PDBQT conversion (AutoDockTools) and protonation options
* Queue templates for SLURM/PBS
* Job status callbacks & minimal results browser

---

## License

MIT (see `LICENSE`).
