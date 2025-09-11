# ==============================
# packager.py  (RENOVATED)
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

# ---------- CSV / job bundle ----------
def write_centers_csv_row(csv_path: Path, pdbqt_name: str, center_xyz, size: float):
    csv_path = Path(csv_path)
    rows = []
    seen = False
    if csv_path.exists():
        with csv_path.open("r", newline="") as f:
            r = csv.DictReader(f)
            hdr = r.fieldnames or []
            for row in r:
                rows.append(row)
    else:
        hdr = ["PDBQT", "center_x", "center_y", "center_z", "size"]

    cx, cy, cz = map(float, center_xyz)
    for row in rows:
        key = row.get("PDBQT") or row.get("pdbqt")
        if key == pdbqt_name:
            row.update({"PDBQT": pdbqt_name, "center_x": cx, "center_y": cy, "center_z": cz, "size": float(size)})
            seen = True
            break
    if not seen:
        rows.append({"PDBQT": pdbqt_name, "center_x": cx, "center_y": cy, "center_z": cz, "size": float(size)})

    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["PDBQT", "center_x", "center_y", "center_z", "size"])
        w.writeheader()
        w.writerows(rows)

# --------- Helpers for PDBQT normalization ----------
def _strip_converted_in_tree(root: Path):
    """Rename files like X.converted.pdbqt -> X.pdbqt in a tree."""
    for p in root.rglob("*"):
        if p.is_file() and ".converted." in p.name and p.suffix.lower() == ".pdbqt":
            new_name = p.name.replace(".converted.", ".")
            p.rename(p.with_name(new_name))

def _copy_tree(src: Path, dst: Path):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)

def _find_pdbqt_dir(ws: Path) -> Path | None:
    # prefer Receptors_PDBQT, then Receptors_PDBQT_Converted
    cand1 = ws / "Receptors_PDBQT"
    cand2 = ws / "Receptors_PDBQT_Converted"
    if cand1.exists(): return cand1
    if cand2.exists(): return cand2
    return None

# ---------- assemble job tree ----------
def assemble_job_tree(ws: Path, rec_dir: Path, lig_dir: Path) -> Path:
    """
    Build job directory structure with normalized receptors:

    - If a PDBQT directory (Receptors_PDBQT or Receptors_PDBQT_Converted) exists,
      use it and rename files to drop '.converted.'; copy as job/Receptors.
    - Otherwise copy raw Receptors (best-effort fallback).
    - Copy Ligands and vina_centers*.csv.
    - Write tool placeholders (unchanged).
    """
    jobroot = ws / "job"
    jobroot.mkdir(exist_ok=True)

    # --- receptors ---
    pdbqt_src = _find_pdbqt_dir(ws)
    if pdbqt_src:
        tmp_norm = ws / "_tmp_norm_pdbqt"
        _copy_tree(pdbqt_src, tmp_norm)
        _strip_converted_in_tree(tmp_norm)
        _copy_tree(tmp_norm, jobroot / "Receptors")
        shutil.rmtree(tmp_norm, ignore_errors=True)
    else:
        # fallback: copy raw receptors
        if rec_dir.exists():
            _copy_tree(rec_dir, jobroot / "Receptors")

    # --- ligands ---
    if lig_dir.exists():
        _copy_tree(lig_dir, jobroot / "Ligands")

    # --- centers CSV(s) ---
    for c in ws.glob("vina_centers*.csv"):
        shutil.copy2(c, jobroot / c.name)

    # --- lightweight tools / placeholders ---
    tdir = jobroot / TOOLS_DIRNAME
    tdir.mkdir(exist_ok=True)
    (jobroot / "1_ConformerGeneration.py").write_text("# placeholder\n")
    (jobroot / "3_Complete_batch_docking.py").write_text("# placeholder\n")
    return jobroot

def zip_job_tree(jobroot: Path) -> Path:
    zpath = jobroot.with_suffix(".zip")
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

    return { "pdb_path": str(rec_pdb), "receptor_pdb": str(rec_pdb) }

# ---------- (optional) single-file prep hook ----------
def prep_receptor_to_pdbqt(single_pdb: Path, out_dir: Path,
                            remove_het_codes: Iterable[str] | str = "all",
                            remove_chains: Iterable[str] | None = None,
                            keep_clean_pdb: bool = True) -> Path:
    """
    Calls 3a_PDB2PDBQTbatch.py in headless mode for ONE file.
    Returns the output directory path where PDBQT landed.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("3a_PDB2PDBQTbatch.py")
    folder = single_pdb.parent
    files = [single_pdb.name]

    cfg = {
        files[0]: {
            "remove_het": ("all" if remove_het_codes == "all" else list(remove_het_codes)),
            "remove_chains": list(remove_chains or [])
        }
    }
    cfg_path = out_dir / "_perfile_config.json"
    cfg_path.write_text(json.dumps(cfg))

    # Prefer ADT first when auto
    args = [
        sys.executable, str(script),
        "--folder", str(folder),
        "--headless",
        "--mode", "per-file",
        "--per-file-config", str(cfg_path),
        "--files", files[0],
        "--output-dir", str(out_dir),
        "--backend", "auto",
        "--altloc", "collapse",
    ]
    if keep_clean_pdb:
        args.append("--keep-clean-pdb")

    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out_dir / "_prep_log.txt").write_text(proc.stdout)
    return out_dir

# ---------- CSV tag helpers (no-op unless you use them) ----------
def collect_receptor_tags(csv_path: Path) -> dict[str, str]:
    tags: dict[str,str] = {}
    if not csv_path.exists():
        return tags
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        headers = [h.strip() for h in reader.fieldnames or []]
        legacy = ("PDB_ID" in headers)
        modern = ("receptor_pdbqt" in headers)
        for row in reader:
            if legacy:
                base = Path(row["PDB_ID"]).stem
            elif modern:
                base = Path(row["receptor_pdbqt"]).stem
            else:
                continue
            tag = row.get("TAG") or ""
            tags[base] = tag
    return tags

def rename_centers_with_tags(jobroot: Path):
    # placeholder to keep compatibility with your previous import
    return
