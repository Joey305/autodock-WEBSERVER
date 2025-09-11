# ==============================
# app.py — AutoDock-Vina PrepServer (finalized: correct converted detection, status flow)
# ==============================
from __future__ import annotations

import os, json, time, csv, subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

from flask import (
    Flask, render_template, request, redirect, url_for, send_file, jsonify, current_app
)
from flask_login import LoginManager, login_required, current_user

from models import db, User
from auth import auth_bp
from lsf_templates import build_confgen_lsfs, build_vina_lsfs
from packager import (
    make_workspace, ensure_subdir, save_uploaded_zip,
    assemble_job_tree, zip_job_tree, fetch_pdb_and_prep, rename_centers_with_tags
)

# ---------- Config ----------
class Config:
    SECRET_KEY = os.getenv("PORTAL_SECRET", "dev-secret-change-me")
    SQLALCHEMY_DATABASE_URI = os.getenv("PORTAL_DB", "sqlite:///portal.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    TMP_ROOT = os.getenv("PORTAL_TMP", "/tmp/autodock_prep")
    MAX_RECEPTORS = int(os.getenv("PORTAL_MAX_RECEPTORS", "5"))
    DEFAULT_ENV_LINE = os.getenv("PORTAL_ENV_LINE", "")

login_manager = LoginManager()
login_manager.login_view = "auth.login"

APP_ROOT = Path(__file__).resolve().parent

def create_app() -> Flask:
    app = Flask(__name__)
    app.config.from_object(Config)
    Path(app.config["TMP_ROOT"]).mkdir(parents=True, exist_ok=True)

    db.init_app(app)
    with app.app_context():
        db.create_all()

    login_manager.init_app(app)

    @login_manager.user_loader
    def _load(uid):
        return db.session.get(User, int(uid))

    app.register_blueprint(auth_bp)

    # ---------- PAGES ----------
    @app.get("/")
    def home():
        if not current_user.is_authenticated:
            return redirect(url_for("auth.login"))
        return render_template("home.html", max_receptors=app.config["MAX_RECEPTORS"])

    # ---------- STATE HELPERS ----------
    def _ws(jobname: str) -> Path:
        return Path(current_app.config["TMP_ROOT"]) / jobname

    def _load_state(ws: Path) -> Dict[str, Any]:
        s = ws / "_state.json"
        if not s.exists():
            return {"receptors": [], "centers_csv": "vina_centers.csv",
                    "ligands_uploaded": False, "prep_job": None}
        return json.loads(s.read_text())

    def _save_state(ws: Path, obj: Dict[str, Any]):
        (ws / "_state.json").write_text(json.dumps(obj, indent=2))

    # ---------- CENTERS CSV HELPERS ----------
    def _centers_csv_path(ws: Path, st: Dict[str, Any]) -> Path:
        name = st.get("centers_csv") or "vina_centers.csv"
        p = ws / name
        if not p.exists():
            with p.open("w", newline="") as f:
                csv.writer(f).writerow(["receptor_pdbqt", "center_x", "center_y", "center_z", "size"])
        return p

    def _read_centers(ws: Path, st: Dict[str, Any]) -> Dict[str, Tuple[float, float, float, float]]:
        p = _centers_csv_path(ws, st)
        out: Dict[str, Tuple[float, float, float, float]] = {}
        with p.open() as f:
            r = csv.DictReader(f)
            for row in r:
                try:
                    out[row["receptor_pdbqt"]] = (
                        float(row["center_x"]), float(row["center_y"]),
                        float(row["center_z"]), float(row["size"])
                    )
                except Exception:
                    continue
        return out

    def _write_centers(ws: Path, st: Dict[str, Any], mapping: Dict[str, Tuple[float, float, float, float]]):
        p = _centers_csv_path(ws, st)
        with p.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["receptor_pdbqt", "center_x", "center_y", "center_z", "size"])
            for k, (x, y, z, s) in mapping.items():
                w.writerow([k, x, y, z, s])

    def _upsert_center_row(ws: Path, st: Dict[str, Any], receptor_pdbqt: str,
                           center: Tuple[float, float, float], size: float):
        mapping = _read_centers(ws, st)
        mapping[receptor_pdbqt] = (float(center[0]), float(center[1]),
                                   float(center[2]), float(size))
        _write_centers(ws, st, mapping)

    # ---------- WORKSPACE ----------
    @app.post("/api/workspace")
    @login_required
    def api_workspace():
        stamp = time.strftime("%m-%d-%Y-%H-%M-%S")
        uname = current_user.email.split("@")[0]
        jobname = f"{stamp}-{uname}"
        ws = make_workspace(_ws(jobname))
        _save_state(ws, {"receptors": [], "centers_csv": "vina_centers.csv",
                         "ligands_uploaded": False, "prep_job": None})
        return jsonify({"jobname": jobname, "workspace": str(ws)})

    # ---------- RECEPTOR INGEST ----------
    @app.post("/api/receptors/upload")
    @login_required
    def api_receptors_upload():
        jobname = request.form.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        mode = request.form.get("mode")
        f = request.files.get("file")
        if not f:
            return ("no file", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        added: List[str] = []
        if mode == "single":
            out = rec_dir / Path(f.filename).name
            out.parent.mkdir(parents=True, exist_ok=True)
            f.save(out)
            if out.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                added.append(str(Path("Receptors") / out.name))
        elif mode == "zip":
            save_uploaded_zip(f, rec_dir)
            for p in sorted(rec_dir.rglob("*")):
                if p.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                    added.append(str(Path("Receptors") / p.relative_to(rec_dir)))
        else:
            return ("bad mode", 400)

        st = _load_state(ws)
        have = {r["rel"] for r in st["receptors"]}
        for rel in added:
            if rel not in have and len(st["receptors"]) < current_app.config["MAX_RECEPTORS"]:
                st["receptors"].append({"rel": rel, "display": Path(rel).name, "status": "new"})
        _save_state(ws, st)
        return jsonify({"count": len(st["receptors"]), "receptors": st["receptors"]})

    @app.post("/api/receptors/fetch")
    @login_required
    def api_receptors_fetch():
        jobname = request.form.get("jobname", "")
        pdbid = (request.form.get("pdb") or "").strip()
        chains = (request.form.get("chains") or "").strip()
        if not pdbid:
            return ("missing pdb", 400)
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        out = fetch_pdb_and_prep(pdbid, rec_dir, chains=chains)
        rel = str(Path("Receptors") / Path(out["pdb_path"]).name)

        st = _load_state(ws)
        if len(st["receptors"]) < current_app.config["MAX_RECEPTORS"]:
            st["receptors"].append({"rel": rel, "display": Path(rel).name, "status": "new"})
        _save_state(ws, st)
        return jsonify({"rel": rel, "count": len(st["receptors"])})

    @app.get("/api/receptors/list")
    @login_required
    def api_receptors_list():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        return jsonify(_load_state(ws))

    # Serve a file from workspace (for NGL)
    @app.get("/api/wsfile")
    @login_required
    def api_wsfile():
        jobname = request.args.get("jobname", "")
        rel = request.args.get("rel", "")
        ws = _ws(jobname)
        if not jobname or not rel or not ws.exists():
            return ("missing", 400)
        p = (ws / rel).resolve()
        if not str(p).startswith(str(ws)) or not p.exists():
            return ("not found", 404)
        return send_file(p)

    # ---------- HET COUNTS ----------
    @app.get("/api/het_counts")
    @login_required
    def api_het_counts():
        jobname = request.args.get("jobname","")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        st = _load_state(ws)

        std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
               "SEC","PYL","MSE"}
        counts: Dict[str,int] = {}
        for rec in st.get("receptors", []):
            rel = rec.get("rel","")
            p = (ws / rel).resolve()
            if not p.exists() or not p.is_file():
                continue
            with open(p) as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        rn = line[17:20].strip().upper()
                        if rn not in std and line.startswith("HETATM"):
                            counts[rn] = counts.get(rn, 0) + 1
        return jsonify({"het_counts": counts})

    # ---------- CENTER ----------
    @app.post("/api/receptor/center")
    @login_required
    def api_receptor_center():
        data = request.get_json(force=True)
        jobname = data.get("jobname")
        rel = data.get("rel")
        size = float(data.get("size") or 20.0)
        method = (data.get("method") or "").lower()
        ws = _ws(jobname)
        src = (ws / rel).resolve()
        if not src.exists():
            return ("not found", 404)

        def _avg(pts: List[Tuple[float, float, float]]):
            n = max(1, len(pts))
            return (sum(p[0] for p in pts) / n,
                    sum(p[1] for p in pts) / n,
                    sum(p[2] for p in pts) / n)

        center: Optional[Tuple[float, float, float]] = None

        if method == "xyz":
            arr = data.get("center") or []
            if len(arr) == 3:
                center = (float(arr[0]), float(arr[1]), float(arr[2]))

        elif method == "selection_atoms":
            pts = [(float(x), float(y), float(z)) for x, y, z in (data.get("atoms") or [])]
            if pts:
                center = _avg(pts)

        elif method in ("ligand", "residue"):
            lig_code = (data.get("lig_resname") or "").upper()
            resname = (data.get("resname") or "").upper()
            resi = (data.get("resi") or "").strip()
            chain = (data.get("chain") or "").upper()

            def match(line: str) -> bool:
                if not line.startswith(("ATOM", "HETATM")):
                    return False
                rn = line[17:20].strip().upper()
                ch = line[21].upper()
                idx = line[22:26].strip()
                if method == "ligand" and lig_code:
                    ok = (rn == lig_code)
                    if resi:  ok = ok and (idx == resi)
                    if chain: ok = ok and (ch == chain)
                    return ok
                if method == "residue" and resname and resi:
                    ok = (rn == resname and idx == resi)
                    if chain: ok = ok and (ch == chain)
                    return ok
                return False

            pts = []
            with open(src) as f:
                for line in f:
                    if match(line):
                        try:
                            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                            pts.append((x, y, z))
                        except Exception:
                            pass
            if pts:
                center = _avg(pts)

        if center is None:
            return ("could not compute center", 400)

        pdbqt_name = Path(rel).name
        pdbqt_name = (pdbqt_name.replace(".pdb", ".pdbqt")
                                  .replace(".cif", ".pdbqt")
                                  .replace(".mmcif", ".pdbqt")
                                  .replace(".ent", ".pdbqt"))
        st = _load_state(ws)
        _upsert_center_row(ws, st, pdbqt_name, center, size)

        # status becomes CENTERED here; only becomes PREPPED after 3a
        for r in st["receptors"]:
            if r["rel"] == rel:
                r["status"] = "centered"
        _save_state(ws, st)

        # all receptors centered?
        expected = []
        for r in st.get("receptors", []):
            nm = Path(r["rel"]).name
            nm = nm.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt")
            expected.append(nm)
        csv_map = _read_centers(ws, st)
        all_centered = (len(expected) > 0 and all(n in csv_map for n in expected))

        return jsonify({"center": center, "size": size,
                        "csv": _centers_csv_path(ws, st).name,
                        "all_centered": all_centered})

    # ---------- helpers for converted detection ----------
    def _stem(p: Path) -> str:
        # basename without suffix; "HDGF_7hg9.converted.pdbqt" -> "HDGF_7hg9.converted"
        return p.name[:-len(p.suffix)] if p.suffix else p.name

    def _expected_stem(rel_name: str) -> str:
        # "HDGF_7hg9.pdb" -> "HDGF_7hg9"
        return Path(rel_name).with_suffix("").name

    def _file_matches_receptor(receptor_rel: str, pdbqt_file: Path) -> bool:
        """
        True if pdbqt filename looks like it was generated from the receptor.
        Accepts NAME.pdbqt, NAME.converted.pdbqt, NAME.anything.pdbqt
        where NAME is the receptor base (without original extension).
        """
        base = _expected_stem(receptor_rel)
        return _stem(pdbqt_file).startswith(base)

    def _mark_prepped_from_output(ws: Path, st: Dict[str, Any], out_dir: Path):
        """Update in-memory state if corresponding PDBQT exists."""
        changed = False
        files = list(out_dir.glob("*.pdbqt"))
        for r in st.get("receptors", []):
            if r.get("status") == "prepped":
                continue
            for f in files:
                if _file_matches_receptor(Path(r["rel"]).name, f):
                    r["status"] = "prepped"
                    changed = True
                    break
        if changed:
            _save_state(ws, st)

    # ---------- SUMMARY ----------
    @app.get("/api/summary")
    @login_required
    def api_summary():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        st = _load_state(ws)

        out_dir = (ensure_subdir(ws, "Receptors").parent / "Receptors_PDBQT").resolve()
        if out_dir.exists():
            _mark_prepped_from_output(ws, st, out_dir)
            st = _load_state(ws)  # re-read after potential update

        total = len(st.get("receptors", []))
        centered = sum(1 for r in st.get("receptors", []) if r.get("status") in ("centered","prepped"))
        prepped  = sum(1 for r in st.get("receptors", []) if r.get("status") == "prepped")
        csv_map = _read_centers(ws, st)
        csv_rows = len(csv_map)

        expected = []
        for r in st.get("receptors", []):
            nm = Path(r["rel"]).name
            nm = nm.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt")
            expected.append(nm)

        converted = [p.name for p in sorted(out_dir.glob("*.pdbqt"))] if out_dir.exists() else []

        all_receptors_centered = (total > 0 and centered == total)
        have_rows_for_all = all(n in csv_map for n in expected) if expected else False

        return jsonify({
            "jobname": jobname,
            "receptors_total": total,
            "receptors_centered": centered,
            "receptors_prepped": prepped,
            "centers_csv": str(_centers_csv_path(ws, st)),
            "centers_rows": csv_rows,
            "expected_pdbqt": expected,
            "have_rows_for_all": have_rows_for_all,
            "conversion_out_dir": str(out_dir),
            "converted_count": len(converted),
            "converted_list": converted,
            "all_receptors_centered": all_receptors_centered
        })

    # ---------- 3a CONVERSION (batch) ----------
    @app.post("/api/prep/start")
    @login_required
    def api_prep_start():
        f = request.form
        jobname = f.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        rec_dir = ensure_subdir(ws, "Receptors")

        st = _load_state(ws)
        recs = st.get("receptors", [])
        if not recs:
            return ("Add receptors first.", 400)
        csv_map = _read_centers(ws, st)
        expected = [Path(r["rel"]).name for r in recs]
        expected = [n.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt") for n in expected]
        if not all(n in csv_map for n in expected):
            return ("Centers CSV missing one or more receptors. Save a center for each receptor first.", 400)

        py = (APP_ROOT / "3a_PDB2PDBQTbatch.py").resolve()
        if not py.exists():
            return ("3a_PDB2PDBQTbatch.py not found at " + str(py), 500)

        out_dir = (rec_dir.parent / "Receptors_PDBQT").resolve()
        log_path = (ws / "prep3a.log").resolve()

        cmd = [
            "python", str(py),
            "--folder", str(rec_dir.resolve()),
            "--headless",
            "--mode", "batch",
            "--backend", f.get("backend", "auto"),
            "--altloc", f.get("altloc", "collapse"),
            "--output-dir", str(out_dir),
        ]

        rm = (f.get("remove_het_csv", "").strip() or "all")
        if rm.lower() == "all":
            cmd += ["--remove-het", "all"]
        elif rm:
            cmd += ["--remove-het", rm]

        ch = f.get("remove_chains_csv", "").strip()
        if ch:
            cmd += ["--remove-chains", ch]

        pdbqt_only = (f.get("pdbqt_only", "true").lower() != "false")
        cmd += ["--pdbqt-only"] if pdbqt_only else ["--no-pdbqt-only"]
        if f.get("keep_clean_pdb", "false").lower() == "true":
            cmd += ["--keep-clean-pdb"]

        with open(log_path, "w") as logf:
            proc = subprocess.Popen(cmd, cwd=str(ws), stdout=logf, stderr=subprocess.STDOUT)

        st["prep_job"] = {"pid": proc.pid, "log": str(log_path), "out_dir": str(out_dir)}
        _save_state(ws, st)
        return jsonify({"pid": proc.pid})

    # ---------- 3a CONVERSION (single) ----------
    @app.post("/api/prep/start_one")
    @login_required
    def api_prep_start_one():
        f = request.form
        jobname = f.get("jobname", "")
        rel = f.get("rel", "")
        ws = _ws(jobname)
        if not ws.exists() or not rel:
            return ("workspace or rel missing", 400)
        rec_dir = ensure_subdir(ws, "Receptors")

        st = _load_state(ws)
        csv_map = _read_centers(ws, st)
        name = Path(rel).name
        expected = name.replace(".pdb",".pdbqt").replace(".cif",".pdbqt").replace(".mmcif",".pdbqt").replace(".ent",".pdbqt")
        if expected not in csv_map:
            return ("Save a center for this receptor first.", 400)

        py = (APP_ROOT / "3a_PDB2PDBQTbatch.py").resolve()
        if not py.exists():
            return ("3a_PDB2PDBQTbatch.py not found at " + str(py), 500)

        out_dir = (rec_dir.parent / "Receptors_PDBQT").resolve()
        log_path = (ws / "prep3a.log").resolve()

        cmd = [
            "python", str(py),
            "--folder", str(rec_dir.resolve()),
            "--headless",
            "--mode", "batch",
            "--backend", f.get("backend", "auto"),
            "--altloc", f.get("altloc", "collapse"),
            "--output-dir", str(out_dir),
            "--files", name
        ]

        rm = (f.get("remove_het_csv", "").strip() or "all")
        if rm.lower() == "all":
            cmd += ["--remove-het", "all"]
        elif rm:
            cmd += ["--remove-het", rm]

        ch = f.get("remove_chains_csv", "").strip()
        if ch:
            cmd += ["--remove-chains", ch]

        pdbqt_only = (f.get("pdbqt_only", "true").lower() != "false")
        cmd += ["--pdbqt-only"] if pdbqt_only else ["--no-pdbqt-only"]

        with open(log_path, "a") as logf:
            proc = subprocess.Popen(cmd, cwd=str(ws), stdout=logf, stderr=subprocess.STDOUT)

        st["prep_job"] = {"pid": proc.pid, "log": str(log_path), "out_dir": str(out_dir)}
        _save_state(ws, st)
        return jsonify({"pid": proc.pid})

    # ---------- 3a STATUS ----------
    @app.get("/api/prep/status")
    @login_required
    def api_prep_status():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        st = _load_state(ws)
        info = st.get("prep_job") or {}
        if not info:
            return jsonify({"running": False, "done": False, "log": ""})

        pid = info.get("pid")
        try:
            os.kill(pid, 0)
            running = True
        except Exception:
            running = False

        lp = Path(info.get("log", ""))
        log = lp.read_text()[-8000:] if lp.exists() else ""
        out_dir = Path(info.get("out_dir", "")).resolve()

        # Determine "done": each receptor must have at least one matching *.pdbqt
        st_now = _load_state(ws)
        receptor_names = [Path(r["rel"]).name for r in st_now.get("receptors", [])]
        files = list(out_dir.glob("*.pdbqt")) if out_dir.exists() else []
        matched = 0
        for rel in receptor_names:
            if any(_file_matches_receptor(rel, f) for f in files):
                matched += 1
        done = (matched == len(receptor_names) and matched > 0)

        # Mark prepped immediately when matching files appear
        if files:
            _mark_prepped_from_output(ws, st_now, out_dir)

        return jsonify({"running": running, "done": done, "log": log})

    # ---------- LIGANDS ----------
    @app.post("/api/ligands/upload")
    @login_required
    def api_lig_upload():
        jobname = request.form.get("jobname", "")
        mode = request.form.get("mode")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        lig_dir = ensure_subdir(ws, "Ligands")
        f = request.files.get("file")
        if not f:
            return ("no file", 400)

        if mode == "single":
            out = lig_dir / Path(f.filename).name
            out.parent.mkdir(parents=True, exist_ok=True)
            f.save(out)
        elif mode == "zip":
            save_uploaded_zip(f, lig_dir)
        else:
            return ("bad mode", 400)

        st = _load_state(ws)
        st["ligands_uploaded"] = True
        _save_state(ws, st)
        return jsonify({"ok": True})

    # ---------- BUILD ----------
    @app.post("/api/build")
    @login_required
    def api_build():
        f = request.form
        jobname = f.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)

        st = _load_state(ws)
        if not st["receptors"]:
            return ("add receptors first", 400)
        csv_map = _read_centers(ws, st)
        expected = [Path(r["rel"]).name for r in st["receptors"]]
        expected = [n.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt") for n in expected]
        if not all(n in csv_map for n in expected):
            return ("centers csv is incomplete", 400)
        if not st.get("ligands_uploaded"):
            return ("upload ligands first", 400)

        jobroot = assemble_job_tree(ws, ws / "Receptors", ws / "Ligands")

        # Ensure LSF directory exists BEFORE template writers touch it
        lsf_dir = (jobroot / "LSF")
        lsf_dir.mkdir(parents=True, exist_ok=True)

        poses_conf = int(f.get("poses_conf", 64))
        poses_vina = int(f.get("poses_vina", 9))

        build_confgen_lsfs(
            jobroot, lsf_dir, poses=poses_conf,
            workers=int(f.get("workers", 16)), queue=f.get("queue", "hihg"),
            project=f.get("project", "brd"), walltime=f.get("walltime", "96:00"),
            mem_per_core=int(f.get("mem_per_core", 2000)), email=current_user.email,
            env_line=f.get("env_line") or current_app.config["DEFAULT_ENV_LINE"],
            lig_mode=None, lig_filetype=None, csv_smiles_col=f.get("csv_smiles_col", ""),
            csv_id_col=f.get("csv_id_col", ""), single_sdf_rel=None
        )
        build_vina_lsfs(
            jobroot, lsf_dir, poses=poses_vina,
            workers=int(f.get("workers", 16)), queue=f.get("queue", "hihg"),
            project=f.get("project", "brd"), walltime=f.get("walltime", "96:00"),
            mem_per_core=int(f.get("mem_per_core", 2000)), email=current_user.email,
            env_line=f.get("env_line") or current_app.config["DEFAULT_ENV_LINE"],
            vina_path=f.get("vina_path") or None
        )

        rename_centers_with_tags(jobroot)
        z = zip_job_tree(jobroot)
        return jsonify({"zip": str(z)})

    # ---------- DOWNLOAD ----------
    @app.get("/download")
    @login_required
    def download():
        path = request.args.get("path", "")
        p = Path(path)
        if not p.exists():
            return ("not found", 404)
        return send_file(p, as_attachment=True, download_name=p.name)

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5050, debug=True)
