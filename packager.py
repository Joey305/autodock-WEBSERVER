# ==============================
# packager.py (final version)
# ==============================
from __future__ import annotations
from pathlib import Path
import os, io, zipfile, shutil, json, subprocess, sys
from typing import Dict, Any, Iterable
import csv

TOOLS_DIRNAME = "tools"

# ---------- workspace helpers ----------
def make_workspace(ws: Path) -> Path:
    ws.mkdir(parents=True, exist_ok=True)
    (ws / "Receptors").mkdir(exist_ok=True)
    (ws / "Ligands").mkdir(exist_ok=True)
    return ws

def ensure_subdir(ws: Path, name: str) -> Path:
    p = ws / name
    p.mkdir(exist_ok=True)
    return p

def save_uploaded_zip(file_storage, dest_dir: Path) -> str:
    buf = file_storage.read()
    with zipfile.ZipFile(io.BytesIO(buf)) as z:
        z.extractall(dest_dir)
    return str(dest_dir)

def move_dir(src: Path, dst: Path):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.move(str(src), str(dst))

# ---------- centers CSV utilities (canonical: PDB_ID,X,Y,Z,SIZE) ----------
def _read_any_centers(csv_path: Path) -> list[dict]:
    rows: list[dict] = []
    if not csv_path.exists():
        return rows
    with csv_path.open("r", newline="") as f:
        r = csv.DictReader(f)
        fields = { (h or "").strip().upper(): h for h in (r.fieldnames or []) }
        has_new = all(k in fields for k in ("PDB_ID","X","Y","Z"))
        has_old = all(k in fields for k in ("RECEPTOR_PDBQT","CENTER_X","CENTER_Y","CENTER_Z"))
        for row in r:
            try:
                if has_new:
                    key = (row.get(fields["PDB_ID"]) or "").strip()
                    x = float(row.get(fields["X"], "0")); y = float(row.get(fields["Y"], "0"))
                    z = float(row.get(fields["Z"], "0")); s = float(row.get(fields.get("SIZE","SIZE"), row.get("SIZE","20")))
                elif has_old:
                    key = (row.get(fields["RECEPTOR_PDBQT"]) or "").strip()
                    x = float(row.get(fields["CENTER_X"], "0")); y = float(row.get(fields["CENTER_Y"], "0"))
                    z = float(row.get(fields["CENTER_Z"], "0")); s = float(row.get(fields.get("SIZE","SIZE"), row.get("size","20")))
                else:
                    continue
                if key:
                    rows.append({"PDB_ID": key, "X": x, "Y": y, "Z": z, "SIZE": s})
            except Exception:
                continue
    return rows

def _write_canonical_centers(csv_path: Path, rows: list[dict]):
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["PDB_ID", "X", "Y", "Z", "SIZE"])
        w.writeheader()
        for r in rows:
            w.writerow({
                "PDB_ID": r["PDB_ID"],
                "X": float(r["X"]), "Y": float(r["Y"]), "Z": float(r["Z"]),
                "SIZE": float(r.get("SIZE", 20.0)),
            })

def write_centers_csv_row(csv_path: Path, pdbqt_name: str, center_xyz, size: float):
    csv_path = Path(csv_path)
    rows = _read_any_centers(csv_path)
    cx, cy, cz = map(float, center_xyz)

    found = False
    for r in rows:
        if r["PDB_ID"] == pdbqt_name:
            r.update({"X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})
            found = True
            break
    if not found:
        rows.append({"PDB_ID": pdbqt_name, "X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})

    _write_canonical_centers(csv_path, rows)

def normalize_centers_csv_to_canonical(src_csv: Path, dst_csv: Path | None = None) -> Path:
    src_csv = Path(src_csv)
    dst_csv = Path(dst_csv) if dst_csv else src_csv
    rows = _read_any_centers(src_csv)
    if not rows and not src_csv.exists():
        _write_canonical_centers(dst_csv, [])
    else:
        _write_canonical_centers(dst_csv, rows)
    return dst_csv

# ---------- runtime script staging ----------
RUNTIME_ROOT_FILES = [
    "1_ConformerGeneration.py",
    "3_Complete_batch_docking.py",
    "4A_PARSEvinascores.py",     # already matches
    "5C_BuildPymolSesh.py",
    "3B_ServerDocks.py",
    "1B_confgen_batch.py", 
    "5C_BuildPymolSesh.py",      # already matches
    # "4A_Parse_VinaResults.py",   # ADD this
    # "4B_Parse_VinaResults.py",
    "6_MDpymacs.py"   # ADD this
    "7_Graphs.py"
]

RUNTIME_TOOLS_FILES = [
    "3B_ServerDocks.py",
    "3a_PDB2PDBQTbatch.py",
    "lsf_templates.py",
    "runDOCKING-tmux.sh",
    "4A_PARSEvinascores.py",     # already matches
    "5C_BuildPymolSesh.py",      # already matches
    # "4A_Parse_VinaResults.py",   # ADD this
    # "4B_Parse_VinaResults.py",
    "6_MDpymacs.py"   # ADD this
    "7_Graphs.py"
]


RUNTIME_TOOL_DIRS = [
    "AutoDockTools_py3",
]

def _safe_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    if dst.suffix in {".sh", ".py"}:
        try:
            mode = os.stat(dst).st_mode
            os.chmod(dst, mode | 0o111)
        except Exception:
            pass

def _stage_runtime_scripts(jobroot: Path, repo_root: Path):
    tdir = jobroot / TOOLS_DIRNAME
    tdir.mkdir(exist_ok=True)

    for rel in RUNTIME_ROOT_FILES:
        src = repo_root / rel
        if src.exists():
            _safe_copy(src, jobroot / src.name)

    for rel in RUNTIME_TOOLS_FILES:
        src = repo_root / rel
        if src.exists():
            _safe_copy(src, tdir / src.name)

    for rel in RUNTIME_TOOL_DIRS:
        src_dir = repo_root / rel
        if src_dir.exists() and src_dir.is_dir():
            dst_dir = tdir / src_dir.name
            if dst_dir.exists():
                shutil.rmtree(dst_dir)
            shutil.copytree(src_dir, dst_dir)

def _chmod_scripts(root: Path):
    for p in root.rglob("*"):
        if p.is_file() and p.suffix in {".lsf", ".sh"}:
            try:
                mode = os.stat(p).st_mode
                os.chmod(p, mode | 0o111)
            except Exception:
                pass

# ---------- job bundle ----------
def assemble_job_tree(ws: Path, rec_dir: Path, lig_dir: Path) -> Path:
    jobroot = ws / "job"
    if jobroot.exists():
        shutil.rmtree(jobroot)
    jobroot.mkdir(parents=True, exist_ok=True)

    # 1) Ligands
    if lig_dir.exists():
        shutil.copytree(lig_dir, jobroot / "Ligands", dirs_exist_ok=True)

    # 2) Centers CSV
    centers_exact = ws / "vina_centers.csv"
    candidate = None
    if centers_exact.exists():
        candidate = centers_exact
    else:
        candidates = sorted(ws.glob("vina_centers*.csv"), key=lambda p: p.stat().st_mtime)
        if candidates:
            candidate = candidates[-1]
    if candidate:
        temp_norm = candidate.parent / "._tmp_canonical_vina_centers.csv"
        normalize_centers_csv_to_canonical(candidate, temp_norm)
        shutil.copy2(temp_norm, jobroot / "vina_centers.csv")
        temp_norm.unlink(missing_ok=True)
    else:
        _write_canonical_centers(jobroot / "vina_centers.csv", [])

    # 3) Rename raw receptors folder -> Receptors_PDB
    raw = ws / "Receptors"
    raw_pdb = ws / "Receptors_PDB"
    if raw.exists() and not raw_pdb.exists():
        raw.rename(raw_pdb)

    # 4) Prepared receptors: rename to Receptors
    prepared = None
    for name in ("Receptors_PDBQT", "Receptors_PDBQT_Converted"):
        p = ws / name
        if p.exists() and any(p.glob("*.pdbqt")):
            prepared = p
            break
    if prepared:
        target = ws / "Receptors"
        if target.exists():
            shutil.rmtree(target)
        prepared.rename(target)

        # also copy into jobroot/Receptors
        target_receptors = jobroot / "Receptors"
        target_receptors.mkdir(parents=True, exist_ok=True)
        for src in target.glob("*.pdbqt"):
            dst_name = src.name.replace(".converted.", ".")
            shutil.copy2(src, target_receptors / dst_name)

    # 5) Also include raw receptors in jobroot
    if raw_pdb.exists():
        shutil.copytree(raw_pdb, jobroot / "Receptors_PDB", dirs_exist_ok=True)

    # 6) Stage runtime scripts
    repo_root = Path(__file__).resolve().parent
    _stage_runtime_scripts(jobroot, repo_root)

    # 7) chmod job scripts
    _chmod_scripts(jobroot)

    return jobroot

def zip_job_tree(jobroot: Path) -> Path:
    zpath = jobroot.with_suffix(".zip")
    if zpath.exists():
        zpath.unlink()
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in jobroot.rglob("*"):
            z.write(p, p.relative_to(jobroot.parent))
    return zpath

# ---------- fetch PDB ----------
def fetch_pdb_and_prep(pdb_code: str, dest_dir: Path, chains: str = "") -> Dict[str, Any]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    pdb_code = pdb_code.strip().lower()
    raw_pdb = dest_dir / f"{pdb_code}.pdb"

    import urllib.request
    url = f"https://files.rcsb.org/view/{pdb_code}.pdb"
    with urllib.request.urlopen(url) as resp, open(raw_pdb, "wb") as fout:
        fout.write(resp.read())

    chains = (chains or "").replace(" ", "")
    if chains:
        wanted = set([c.upper() for c in chains.split(",") if c])
        filtered = dest_dir / f"{pdb_code}.pdb"
        with open(raw_pdb) as fin, open(filtered, "w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    ch = line[21].upper()
                    if ch in wanted:
                        fout.write(line)
                else:
                    fout.write(line)
        rec_pdb = filtered
    else:
        rec_pdb = raw_pdb

    return {"pdb_path": str(rec_pdb), "receptor_pdb": str(rec_pdb)}

# ---------- simple cleaner ----------
_STD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

def write_cleaned_pdb(src: Path, out: Path, keep_residues: set[str]) -> None:
    with open(src) as fin, open(out, "w") as fout:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                fout.write(line); continue
            rn = line[17:20].strip().upper()
            if rn in _STD or rn in keep_residues:
                fout.write(line)

# ---------- optional hook into 3a ----------
def prep_receptor_to_pdbqt(single_pdb: Path, out_dir: Path,
                            remove_het_codes: Iterable[str] | str = "all",
                            remove_chains: Iterable[str] | None = None,
                            keep_clean_pdb: bool = True) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("3a_PDB2PDBQTbatch.py")
    folder = single_pdb.parent
    files = [single_pdb.name]
    cfg = {
        files[0]: {
            "remove_het": ("all" if remove_het_codes == "all"
                           else list(remove_het_codes)),
            "remove_chains": list(remove_chains or [])
        }
    }
    cfg_path = out_dir / "_perfile_config.json"
    cfg_path.write_text(json.dumps(cfg))

    args = [
        sys.executable, str(script),
        "--folder", str(folder),
        "--headless", "--mode", "per-file",
        "--per-file-config", str(cfg_path),
        "--files", files[0],
        "--output-dir", str(out_dir),
        "--backend", "auto", "--altloc", "collapse",
    ]
    if keep_clean_pdb:
        args.append("--keep-clean-pdb")

    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out_dir / "_prep_log.txt").write_text(proc.stdout)
    return out_dir

# ---------- tag helpers ----------
def collect_receptor_tags(csv_path: Path) -> dict[str, str]:
    tags: dict[str,str] = {}
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return tags
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        fields = { (h or "").strip().upper(): h for h in (reader.fieldnames or []) }
        key_new = fields.get("PDB_ID")
        key_old = fields.get("RECEPTOR_PDBQT")
        tag_col = fields.get("TAG")
        for row in reader:
            if tag_col:
                if key_new and row.get(key_new):
                    base = Path(row[key_new]).stem
                    tags[base] = row.get(tag_col, "")
                elif key_old and row.get(key_old):
                    base = Path(row[key_old]).stem
                    tags[base] = row.get(tag_col, "")
    return tags

def rename_centers_with_tags(jobroot: Path):
    root_csv = jobroot / "vina_centers.csv"
    _ = collect_receptor_tags(root_csv)
    # no destructive rename in this version



# # ==============================
# # packager.py  (revamped)
# # ==============================
# from __future__ import annotations
# from pathlib import Path
# import os, io, zipfile, shutil, json, subprocess, sys
# from typing import Dict, Any, Iterable
# import csv

# TOOLS_DIRNAME = "tools"

# # ---------- workspace helpers ----------
# def make_workspace(ws: Path) -> Path:
#     ws.mkdir(parents=True, exist_ok=True)
#     (ws / "Receptors").mkdir(exist_ok=True)
#     (ws / "Ligands").mkdir(exist_ok=True)
#     return ws

# def ensure_subdir(ws: Path, name: str) -> Path:
#     p = ws / name
#     p.mkdir(exist_ok=True)
#     return p

# def save_uploaded_zip(file_storage, dest_dir: Path) -> str:
#     buf = file_storage.read()
#     with zipfile.ZipFile(io.BytesIO(buf)) as z:
#         z.extractall(dest_dir)
#     return str(dest_dir)

# def move_dir(src: Path, dst: Path):
#     if dst.exists():
#         shutil.rmtree(dst)
#     shutil.move(str(src), str(dst))

# # ---------- centers CSV utilities (canonical: PDB_ID,X,Y,Z,SIZE) ----------
# def _read_any_centers(csv_path: Path) -> list[dict]:
#     """Read centers CSV accepting either canonical or legacy headers; return rows in canonical dict shape."""
#     rows: list[dict] = []
#     if not csv_path.exists():
#         return rows
#     with csv_path.open("r", newline="") as f:
#         r = csv.DictReader(f)
#         fields = { (h or "").strip().upper(): h for h in (r.fieldnames or []) }
#         has_new = all(k in fields for k in ("PDB_ID","X","Y","Z"))
#         has_old = all(k in fields for k in ("RECEPTOR_PDBQT","CENTER_X","CENTER_Y","CENTER_Z"))
#         for row in r:
#             try:
#                 if has_new:
#                     key = (row.get(fields["PDB_ID"]) or "").strip()
#                     x = float(row.get(fields["X"], "0")); y = float(row.get(fields["Y"], "0"))
#                     z = float(row.get(fields["Z"], "0")); s = float(row.get(fields.get("SIZE","SIZE"), row.get("SIZE","20")))
#                 elif has_old:
#                     key = (row.get(fields["RECEPTOR_PDBQT"]) or "").strip()
#                     x = float(row.get(fields["CENTER_X"], "0")); y = float(row.get(fields["CENTER_Y"], "0"))
#                     z = float(row.get(fields["CENTER_Z"], "0")); s = float(row.get(fields.get("SIZE","SIZE"), row.get("size","20")))
#                 else:
#                     continue
#                 if key:
#                     rows.append({"PDB_ID": key, "X": x, "Y": y, "Z": z, "SIZE": s})
#             except Exception:
#                 continue
#     return rows

# def _write_canonical_centers(csv_path: Path, rows: list[dict]):
#     with csv_path.open("w", newline="") as f:
#         w = csv.DictWriter(f, fieldnames=["PDB_ID", "X", "Y", "Z", "SIZE"])
#         w.writeheader()
#         # ensure consistent order & types
#         for r in rows:
#             w.writerow({
#                 "PDB_ID": r["PDB_ID"],
#                 "X": float(r["X"]), "Y": float(r["Y"]), "Z": float(r["Z"]),
#                 "SIZE": float(r.get("SIZE", 20.0)),
#             })

# def write_centers_csv_row(csv_path: Path, pdbqt_name: str, center_xyz, size: float):
#     """
#     Upsert a row keyed by PDBQT filename using canonical headers:
#         PDB_ID, X, Y, Z, SIZE
#     Also tolerant to legacy files (rewrites as canonical on save).
#     """
#     csv_path = Path(csv_path)
#     rows = _read_any_centers(csv_path)
#     cx, cy, cz = map(float, center_xyz)

#     found = False
#     for r in rows:
#         if r["PDB_ID"] == pdbqt_name:
#             r.update({"X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})
#             found = True
#             break
#     if not found:
#         rows.append({"PDB_ID": pdbqt_name, "X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})

#     _write_canonical_centers(csv_path, rows)

# def normalize_centers_csv_to_canonical(src_csv: Path, dst_csv: Path | None = None) -> Path:
#     """
#     Read src_csv (any schema we support), then write the canonical schema to dst_csv (or overwrite src).
#     """
#     src_csv = Path(src_csv)
#     dst_csv = Path(dst_csv) if dst_csv else src_csv
#     rows = _read_any_centers(src_csv)
#     if not rows and not src_csv.exists():
#         # create a new empty canonical file
#         _write_canonical_centers(dst_csv, [])
#     else:
#         _write_canonical_centers(dst_csv, rows)
#     return dst_csv

# # ---------- runtime script staging ----------
# # Root-level scripts to land directly in jobroot/
# RUNTIME_ROOT_FILES = [
#     "1_ConformerGeneration.py",
#     "3_Complete_batch_docking.py",
# ]

# # Individual files to place under jobroot/tools/
# RUNTIME_TOOLS_FILES = [
#     "3B_ServerDocks.py",
#     "3a_PDB2PDBQTbatch.py",
#     "lsf_templates.py",
#     "runDOCKING-tmux.sh",
#     "4A_Parse_VinaResults.py",
#     "4B_Parse_VinaResults.py",
# ]

# # Whole directories to mirror under jobroot/tools/
# RUNTIME_TOOL_DIRS = [
#     "AutoDockTools_py3",
# ]

# def _safe_copy(src: Path, dst: Path):
#     dst.parent.mkdir(parents=True, exist_ok=True)
#     shutil.copy2(src, dst)
#     # Make shell/python scripts executable for convenience
#     if dst.suffix in {".sh", ".py"}:
#         try:
#             mode = os.stat(dst).st_mode
#             os.chmod(dst, mode | 0o111)
#         except Exception:
#             pass

# def _stage_runtime_scripts(jobroot: Path, repo_root: Path):
#     """
#     Copy selected scripts and tool directories from your repo into the job.
#     - Items in RUNTIME_ROOT_FILES go to jobroot/
#     - Items in RUNTIME_TOOLS_FILES go to jobroot/tools/
#     - Folders in RUNTIME_TOOL_DIRS are mirrored under jobroot/tools/
#     """
#     tdir = jobroot / TOOLS_DIRNAME
#     tdir.mkdir(exist_ok=True)

#     # Root-level entrypoints
#     for rel in RUNTIME_ROOT_FILES:
#         src = repo_root / rel
#         if src.exists():
#             _safe_copy(src, jobroot / src.name)

#     # Individual tool files
#     for rel in RUNTIME_TOOLS_FILES:
#         src = repo_root / rel
#         if src.exists():
#             _safe_copy(src, tdir / src.name)

#     # Tool directories
#     for rel in RUNTIME_TOOL_DIRS:
#         src_dir = repo_root / rel
#         if src_dir.exists() and src_dir.is_dir():
#             dst_dir = tdir / src_dir.name
#             if dst_dir.exists():
#                 shutil.rmtree(dst_dir)
#             shutil.copytree(src_dir, dst_dir)

# def _chmod_scripts(root: Path):
#     """Optional sweep to ensure .lsf/.sh in job get exec bit."""
#     for p in root.rglob("*"):
#         if p.is_file() and p.suffix in {".lsf", ".sh"}:
#             try:
#                 mode = os.stat(p).st_mode
#                 os.chmod(p, mode | 0o111)
#             except Exception:
#                 pass

# # ---------- job bundle ----------
# def _write_placeholders(jobroot: Path):
#     """
#     Fallback if you don't want to stage real scripts.
#     Not used by assemble_job_tree (kept for compatibility).
#     """
#     tdir = jobroot / TOOLS_DIRNAME
#     tdir.mkdir(exist_ok=True)
#     (jobroot / "1_ConformerGeneration.py").write_text("# placeholder\n")
#     (jobroot / "3_Complete_batch_docking.py").write_text("# placeholder\n")

# def assemble_job_tree(ws: Path, rec_dir: Path, lig_dir: Path) -> Path:
#     """
#     Build job/ with:
#       - Receptors/    <- prepared PDBQT set
#       - Ligands/      <- copied as-is
#       - vina_centers.csv (canonicalized to PDB_ID,X,Y,Z,SIZE)
#       - tools/        <- staged real scripts

#     Keep ws/Receptors for the frontend. At the very end, snapshot it to ws/Receptors_PDB.
#     """
#     jobroot = ws / "job"
#     if jobroot.exists():
#         shutil.rmtree(jobroot)
#     jobroot.mkdir(parents=True, exist_ok=True)

#     # 1) Ligands (as-is)
#     if lig_dir.exists():
#         shutil.copytree(lig_dir, jobroot / "Ligands", dirs_exist_ok=True)

#     # 2) Centers CSV: pick best candidate, normalize to canonical, and copy as job/vina_centers.csv
#     centers_exact = ws / "vina_centers.csv"
#     candidate = None
#     if centers_exact.exists():
#         candidate = centers_exact
#     else:
#         candidates = sorted(ws.glob("vina_centers*.csv"), key=lambda p: p.stat().st_mtime)
#         if candidates:
#             candidate = candidates[-1]
#     if candidate:
#         temp_norm = candidate.parent / ("._tmp_canonical_vina_centers.csv")
#         normalize_centers_csv_to_canonical(candidate, temp_norm)
#         shutil.copy2(temp_norm, jobroot / "vina_centers.csv")
#         try:
#             temp_norm.unlink()
#         except Exception:
#             pass
#     else:
#         # create an empty canonical file in jobroot
#         _write_canonical_centers(jobroot / "vina_centers.csv", [])

#     # 3) Prepared receptors -> job/Receptors (strip ".converted.")
#     target_receptors = jobroot / "Receptors"
#     target_receptors.mkdir(parents=True, exist_ok=True)
#     prepared = None
#     for name in ("Receptors_PDBQT", "Receptors_PDBQT_Converted"):
#         p = ws / name
#         if p.exists() and any(p.glob("*.pdbqt")):
#             prepared = p
#             break
#     if prepared:
#         for src in prepared.glob("*.pdbqt"):
#             dst_name = src.name.replace(".converted.", ".")
#             shutil.copy2(src, target_receptors / dst_name)

#     # 4) Stage runtime scripts (real files, not placeholders)
#     repo_root = Path(__file__).resolve().parent
#     _stage_runtime_scripts(jobroot, repo_root)

#     # 5) VERY END: snapshot raw Receptors -> Receptors_PDB (do NOT rename/delete original)
#     raw_receptors = ws / "Receptors"          # used by the frontend; do not touch
#     snapshot = ws / "Receptors_PDB"
#     if raw_receptors.exists():
#         if snapshot.exists():
#             shutil.rmtree(snapshot)
#         shutil.copytree(raw_receptors, snapshot)

#     # 6) Make sure job scripts are executable (belt-and-suspenders)
#     _chmod_scripts(jobroot)

#     return jobroot

# def zip_job_tree(jobroot: Path) -> Path:
#     zpath = jobroot.with_suffix(".zip")
#     if zpath.exists():
#         zpath.unlink()
#     with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
#         for p in jobroot.rglob("*"):
#             # include the 'job/' folder root inside the zip
#             z.write(p, p.relative_to(jobroot.parent))
#     return zpath

# # ---------- fetch PDB ----------
# def fetch_pdb_and_prep(pdb_code: str, dest_dir: Path, chains: str = "") -> Dict[str, Any]:
#     """
#     Download a PDB to dest_dir (ws/Receptors), optionally filter chains (comma list).
#     Returns dict with absolute path to downloaded/filtered PDB.
#     """
#     dest_dir.mkdir(parents=True, exist_ok=True)
#     pdb_code = pdb_code.strip().lower()
#     raw_pdb = dest_dir / f"{pdb_code}.pdb"

#     import urllib.request
#     url = f"https://files.rcsb.org/view/{pdb_code}.pdb"
#     with urllib.request.urlopen(url) as resp, open(raw_pdb, "wb") as fout:
#         fout.write(resp.read())

#     chains = (chains or "").replace(" ", "")
#     if chains:
#         wanted = set([c.upper() for c in chains.split(",") if c])
#         filtered = dest_dir / f"{pdb_code}.pdb"
#         with open(raw_pdb) as fin, open(filtered, "w") as fout:
#             for line in fin:
#                 if line.startswith(("ATOM", "HETATM")):
#                     ch = line[21].upper()
#                     if ch in wanted:
#                         fout.write(line)
#                 else:
#                     fout.write(line)
#         rec_pdb = filtered
#     else:
#         rec_pdb = raw_pdb

#     return {"pdb_path": str(rec_pdb), "receptor_pdb": str(rec_pdb)}

# # ---------- simple cleaner (unused by current front-end; kept for completeness) ----------
# _STD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
#         "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}

# def write_cleaned_pdb(src: Path, out: Path, keep_residues: set[str]) -> None:
#     with open(src) as fin, open(out, "w") as fout:
#         for line in fin:
#             if not line.startswith(("ATOM", "HETATM")):
#                 fout.write(line); continue
#             rn = line[17:20].strip().upper()
#             if rn in _STD or rn in keep_residues:
#                 fout.write(line)

# # ---------- optional hook into 3a for per-file ----------
# def prep_receptor_to_pdbqt(single_pdb: Path, out_dir: Path,
#                             remove_het_codes: Iterable[str] | str = "all",
#                             remove_chains: Iterable[str] | None = None,
#                             keep_clean_pdb: bool = True) -> Path:
#     """
#     Calls 3a_PDB2PDBQTbatch.py in headless mode for ONE file.
#     Returns the output directory path where PDBQT landed.
#     """
#     out_dir.mkdir(parents=True, exist_ok=True)
#     script = Path(__file__).with_name("3a_PDB2PDBQTbatch.py")
#     folder = single_pdb.parent
#     files = [single_pdb.name]
#     cfg = {
#         files[0]: {
#             "remove_het": ("all" if remove_het_codes == "all"
#                            else list(remove_het_codes)),
#             "remove_chains": list(remove_chains or [])
#         }
#     }
#     cfg_path = out_dir / "_perfile_config.json"
#     cfg_path.write_text(json.dumps(cfg))

#     args = [
#         sys.executable, str(script),
#         "--folder", str(folder),
#         "--headless", "--mode", "per-file",
#         "--per-file-config", str(cfg_path),
#         "--files", files[0],
#         "--output-dir", str(out_dir),
#         "--backend", "auto", "--altloc", "collapse",
#     ]
#     if keep_clean_pdb:
#         args.append("--keep-clean-pdb")

#     proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
#     (out_dir / "_prep_log.txt").write_text(proc.stdout)
#     return out_dir

# # ---------- tag helpers (no-op if missing) ----------
# def collect_receptor_tags(csv_path: Path) -> dict[str, str]:
#     """
#     Reads a TAG column if present. Accept PDB_ID or legacy receptor_pdbqt as the key.
#     """
#     tags: dict[str,str] = {}
#     csv_path = Path(csv_path)
#     if not csv_path.exists():
#         return tags
#     with csv_path.open() as f:
#         reader = csv.DictReader(f)
#         fields = { (h or "").strip().upper(): h for h in (reader.fieldnames or []) }
#         key_new = fields.get("PDB_ID")
#         key_old = fields.get("RECEPTOR_PDBQT")
#         tag_col = fields.get("TAG")
#         for row in reader:
#             if tag_col:
#                 if key_new and row.get(key_new):
#                     base = Path(row[key_new]).stem
#                     tags[base] = row.get(tag_col, "")
#                 elif key_old and row.get(key_old):
#                     base = Path(row[key_old]).stem
#                     tags[base] = row.get(tag_col, "")
#     return tags

# def rename_centers_with_tags(jobroot: Path):
#     # retained for compatibility; does nothing unless a TAG column exists
#     root_csv = jobroot / "vina_centers.csv"
#     _ = collect_receptor_tags(root_csv)
#     # no destructive rename in this version
