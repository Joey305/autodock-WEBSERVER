# ==============================
# app.py  — drop-in
# ==============================
from __future__ import annotations

import os
import time
from pathlib import Path

from flask import (
    Flask, render_template, request, redirect, url_for, flash,
    send_file, jsonify, current_app
)
from flask_login import LoginManager, login_required, current_user

from models import db, User
from auth import auth_bp

from packager import (
    make_workspace, save_uploaded_zip, ensure_subdir, write_centers_csv_row,
    assemble_job_tree, zip_job_tree, fetch_pdb_and_prep, move_dir,
    rename_centers_with_tags
)
from lsf_templates import build_confgen_lsfs, build_vina_lsfs


# ---------- config ----------

class Config:
    SECRET_KEY = os.getenv("PORTAL_SECRET", "dev-secret-change-me")
    SQLALCHEMY_DATABASE_URI = os.getenv("PORTAL_DB", "sqlite:///portal.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    TMP_ROOT = os.getenv("PORTAL_TMP", "/tmp/pwwp_prep")
    DEFAULT_EMAIL_DOMAIN = os.getenv("PORTAL_EMAIL_DOMAIN", "@miami.edu")
    DEFAULT_ENV_LINE = os.getenv(
        "PORTAL_ENV_LINE",
        "source /nethome/jxs794/miniconda3/etc/profile.d/conda.sh && conda activate vina_env",
    )


login_manager = LoginManager()
login_manager.login_view = "auth.login"


def create_app() -> Flask:
    app = Flask(__name__)
    app.config.from_object(Config)

    Path(app.config["TMP_ROOT"]).mkdir(parents=True, exist_ok=True)

    db.init_app(app)
    with app.app_context():
        db.create_all()

    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(uid):
        return db.session.get(User, int(uid))

    app.register_blueprint(auth_bp)

    # ---------- routes ----------

    @app.route("/")
    def home():
        if current_user.is_authenticated:
            return render_template("home.html")
        return redirect(url_for("auth.login"))

    # create workspace
    @app.post("/api/workspace")
    @login_required
    def api_workspace():
        stamp = time.strftime("%m-%d-%Y-%H-%M-%S")
        uname = current_user.email.split("@")[0]
        jobname = f"{stamp}-{uname}"
        ws = make_workspace(Path(current_app.config["TMP_ROOT"]) / jobname)
        return jsonify({"jobname": jobname, "workspace": str(ws)})

    # upload endpoints
    @app.post("/api/upload")
    @login_required
    def api_upload():
        kind = request.form.get("kind")  # receptors|receptors_file|ligands|ligands_file|centers
        jobname = request.form.get("jobname")
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)
        file = request.files.get("file")
        if not file:
            return ("no file", 400)

        if kind == "receptors":
            rec_dir = ensure_subdir(ws, "Receptors")
            saved = save_uploaded_zip(file, rec_dir)
            first_rel = None
            pdbs = sorted(rec_dir.rglob("*.pdb"))
            pdbqt = sorted(rec_dir.rglob("*.pdbqt"))
            if pdbs:
                first_rel = str(Path("Receptors") / pdbs[0].relative_to(rec_dir))
            elif pdbqt:
                first_rel = str(Path("Receptors") / pdbqt[0].relative_to(rec_dir))
            return jsonify({"saved": saved, "firstRel": first_rel, "detect": "receptors_zip"})

        if kind == "receptors_file":
            rec_dir = ensure_subdir(ws, "Receptors")
            out = rec_dir / Path(file.filename).name
            file.save(out)
            ext = out.suffix.lower().lstrip(".")
            rel = str(Path("Receptors") / out.name)
            return jsonify({"saved": out.name, "ext": ext, "relpath": rel})

        if kind == "ligands":
            lig_dir = ensure_subdir(ws, "Ligands")
            saved = save_uploaded_zip(file, lig_dir)
            detect = "folder"
            if any(lig_dir.rglob("*.sdf")):
                detect = "folder_sdf"
            elif any(lig_dir.rglob("*.smiles")):
                detect = "folder_smiles"
            elif any(lig_dir.rglob("*.csv")):
                detect = "csv"
            return jsonify({"saved": saved, "detect": detect})

        if kind == "ligands_file":
            lig_dir = ensure_subdir(ws, "Ligands")
            out = lig_dir / Path(file.filename).name
            file.save(out)
            ext = out.suffix.lower().lstrip(".")
            return jsonify({"saved": out.name, "detect": f"single_{ext}"})

        if kind == "centers":
            out = ws / "vina_centers.csv"
            file.save(out)
            return jsonify({"saved": out.name})

        return ("bad kind", 400)

    # RCSB fetch to Receptors/
    @app.post("/api/fetch_pdb")
    @login_required
    def api_fetch_pdb():
        jobname = request.form.get("jobname")
        pdb_code = (request.form.get("pdb") or "").strip()
        chains = (request.form.get("chains") or "").strip()
        if not pdb_code:
            return ("missing pdb", 400)
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)
        rec_dir = ensure_subdir(ws, "Receptors")
        out = fetch_pdb_and_prep(pdb_code, rec_dir, chains=chains)
        rel = str(Path("Receptors") / Path(out["pdb_path"]).name)
        return jsonify({"pdb_rel": rel, "pdbid": pdb_code.upper()})

    # file server for NGL
    @app.get("/api/wsfile")
    @login_required
    def api_wsfile():
        jobname = request.args.get("jobname", "")
        rel = request.args.get("rel", "")
        if not jobname or not rel:
            return ("missing", 400)
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        p = (ws / rel).resolve()
        if not str(p).startswith(str(ws)) or not p.exists() or not p.is_file():
            return ("not found", 404)
        return send_file(p)

    # residue list (for removal panel)
    @app.get("/api/list_residues")
    @login_required
    def api_list_residues():
        jobname = request.args.get("jobname")
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)
        rec_dir = ensure_subdir(ws, "Receptors")

        # prefer cleaned, else latest pdb (search subdirs too)
        cand = sorted(rec_dir.rglob("*.cleaned.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not cand:
            cand = sorted(rec_dir.rglob("*.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not cand:
            return jsonify({"resnames": []})

        std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
        res = set()
        with open(cand[0]) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    rn = line[17:20].strip().upper()
                    if rn not in std:
                        res.add(rn)
        return jsonify({"resnames": sorted(res)})

    # receptor cleanup (rewrite PDB)
    @app.post("/api/clean_receptor")
    @login_required
    def api_clean_receptor():
        data = request.get_json(force=True)
        jobname = data.get("jobname")
        keep_res = set(data.get("keep_resnames") or [])
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        # pick latest pdb
        srcs = sorted(rec_dir.rglob("*.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not srcs:
            return ("no pdb", 400)
        src = srcs[0]
        out = rec_dir / (src.stem + ".cleaned.pdb")

        std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

        with open(src) as fin, open(out, "w") as fout:
            for line in fin:
                if not line.startswith(("ATOM", "HETATM")):
                    fout.write(line)
                    continue
                resname = line[17:20].strip().upper()
                if resname in std:
                    fout.write(line)  # keep protein always
                elif resname in keep_res:
                    fout.write(line)  # user explicitly kept this residue
                # else drop it

        rel = str(Path("Receptors") / out.name)
        return jsonify({"out_rel": rel})

    # center capture (append to vina_centers.csv at workspace root)
    @app.post("/api/center")
    @login_required
    def api_center():
        data = request.get_json(force=True)
        jobname = data.get("jobname")
        receptor_label = data.get("receptor_label") or "receptor"
        center = data.get("center")
        size = float(data.get("size", 20))
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)
        centers_csv = ws / "vina_centers.csv"
        write_centers_csv_row(centers_csv, receptor_label, center, size)
        return jsonify({"wrote": centers_csv.name})

    # build zip (confgen + vina LSFs)
    @app.post("/api/build")
    @login_required
    def api_build():
        f = request.form
        jobname = f.get("jobname")
        ws = Path(current_app.config["TMP_ROOT"]) / jobname
        if not ws.exists():
            return ("workspace missing", 400)

        email = current_user.email
        queue = f.get("queue", "hihg")
        project = f.get("project", "brd")
        walltime = f.get("walltime", "96:00")
        workers = int(f.get("workers", 16))
        mem_per_core = int(f.get("mem_per_core", 2000))
        poses_conf = int(f.get("poses_conf", 64))
        poses_vina = int(f.get("poses_vina", 9))
        env_line = f.get("env_line") or current_app.config["DEFAULT_ENV_LINE"]
        vina_path = f.get("vina_path") or None

        # optional CSV columns (autodetect mode lives inside your confgen template)
        csv_smiles_col = f.get("csv_smiles_col", "")
        csv_id_col = f.get("csv_id_col", "")

        jobroot = assemble_job_tree(ws, ws / "Receptors", ws / "Ligands")
        lsf_dir = jobroot / "LSF"
        lsf_dir.mkdir(exist_ok=True)

        build_confgen_lsfs(
            jobroot,
            lsf_dir,
            poses=poses_conf,
            workers=workers,
            queue=queue,
            project=project,
            walltime=walltime,
            mem_per_core=mem_per_core,
            email=email,
            env_line=env_line,
            # pass through hints for autodetection
            lig_mode=None,
            lig_filetype=None,
            csv_smiles_col=csv_smiles_col,
            csv_id_col=csv_id_col,
            single_sdf_rel=None,
        )

        build_vina_lsfs(
            jobroot,
            lsf_dir,
            poses=poses_vina,
            workers=workers,
            queue=queue,
            project=project,
            walltime=walltime,
            mem_per_core=mem_per_core,
            email=email,
            env_line=env_line,
            vina_path=vina_path,
        )

        # Add a tag-expanded copy of centers if present
        rename_centers_with_tags(jobroot)

        zip_path = zip_job_tree(jobroot)
        return jsonify({"zip": str(zip_path)})

    # download
    @app.get("/download")
    @login_required
    def download():
        path = request.args.get("path")
        if not path:
            return ("missing path", 400)
        p = Path(path)
        if not p.exists():
            return ("not found", 404)
        return send_file(p, as_attachment=True, download_name=p.name)

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5050, debug=True)
