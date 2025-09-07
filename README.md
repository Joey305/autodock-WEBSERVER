
# ==============================
# README.md
# ==============================
# PWWP Web Prep Portal (minimal)

## Quick start
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export FLASK_APP=app.py
python -c "from app import create_app; app=create_app();\nfrom models import db,User;\napp.app_context().push(); db.drop_all(); db.create_all(); u=User(email='jxs794@miami.edu'); u.set_password('Pass!'); db.session.add(u); db.session.commit(); print('admin ok')"
flask run -p 5050
```

## Workflow
1. **New Workspace** → creates `/tmp/pwwp_prep/<MM-DD-YYYY-HH-MM-SS-USER>/`
2. **Upload ZIPs** for `Receptors/` and `Ligands/` (server unzips in workspace)
3. Optionally **Fetch PDB** (writes to Receptors/Ligands)
4. **Open Viewer** with a receptor file → set selection string → *Use selection center* → appends to `vina_centers_<receptor>.csv`
5. **Build LSF + Zip** → builds `LSF/run_*.lsf`, submitters, and bundles a job zip under `/tmp`
6. Download the zip or manually move it to Pegasus.

## Customize
- Default conda env line (`PORTAL_ENV_LINE`) can be overridden via env.
- Optional `VINA_EXE` path can be typed in the UI.
- Replace placeholders for your real `1_ConformerGeneration.py` and `3_Complete_batch_docking.py` when assembling.
