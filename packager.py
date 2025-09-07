# ==============================
# packager.py  — clean, drop-in
# ==============================
from __future__ import annotations

import io
import os
import shutil
import time
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

TOOLS_DIRNAME = "tools"

# ---------- workspace helpers ----------

def make_workspace(ws: Path) -> Path:
    ws.mkdir(parents=True, exist_ok=True)
    (ws / "Receptors").mkdir(exist_ok=True)
    (ws / "Ligands").mkdir(exist_ok=True)
    return ws


def ensure_subdir(ws: Path, name: str) -> Path:
    p = ws / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def save_uploaded_zip(file_storage, dest_dir: Path) -> str:
    """Extract ZIP directly *into* dest_dir (no extra top-level folder)."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    buf = file_storage.read()
    with zipfile.ZipFile(io.BytesIO(buf)) as z:
        z.extractall(dest_dir)
    return str(dest_dir)


def move_dir(src: Path, dst: Path):
    if dst.exists():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))


def write_centers_csv_row(csv_path: Path, receptor_name: str, center, size: float):
    header_needed = not csv_path.exists()
    with open(csv_path, "a") as f:
        if header_needed:
            f.write("PDB_ID,X,Y,Z,SIZE\n")
        x, y, z = center
        f.write(f"{receptor_name},{x:.3f},{y:.3f},{z:.3f},{size}\n")


# ---------- job tree / packaging ----------

def assemble_job_tree(ws: Path, rec_dir: Path, lig_dir: Path) -> Path:
    """Create job/ with copied Receptors/, Ligands/, centers CSV and helper tools."""
    jobroot = ws / "job"
    jobroot.mkdir(exist_ok=True)

    if rec_dir.exists():
        shutil.copytree(rec_dir, jobroot / "Receptors", dirs_exist_ok=True)
    if lig_dir.exists():
        shutil.copytree(lig_dir, jobroot / "Ligands", dirs_exist_ok=True)

    # copy any centers CSV(s) from ws root
    for c in ws.glob("vina_centers*.csv"):
        shutil.copy2(c, jobroot / c.name)

    # placeholders (you’ll drop your real scripts in)
    (jobroot / "1_ConformerGeneration.py").write_text(
        "# PLACEHOLDER — drop your real 1_ConformerGeneration.py here\n"
    )
    (jobroot / "3_Complete_batch_docking.py").write_text(
        "# PLACEHOLDER — drop your real 3_Complete_batch_docking.py here\n"
    )

    # tiny tools folder
    tdir = jobroot / TOOLS_DIRNAME
    tdir.mkdir(exist_ok=True)
    (tdir / "pymol_center_helper.py").write_text(PYMOL_CENTER_HELPER)
    (tdir / "pdb2pdbqt_batch.py").write_text(PDB2PDBQT_BATCH)

    return jobroot


def zip_job_tree(jobroot: Path) -> Path:
    zpath = jobroot.with_suffix(".zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for p in jobroot.rglob("*"):
            z.write(p, p.relative_to(jobroot.parent))
    return zpath


# ---------- PDB fetch/prep (fixed) ----------

def fetch_pdb_and_prep(pdb_code: str, dest_dir: Path, chains: str = "") -> Dict[str, Any]:
    """
    Download PDB into *dest_dir* (already ws/Receptors), optionally filter chains,
    attempt a very-light ligand extraction, and return normalized paths.

    Returns:
      {
        "pdb_path": "/.../<job>/Receptors/3eky.pdb",
        "receptor_pdb": same as above,
        "ligand_sdf": None (or path if you wire in a real converter)
      }
    """
    import urllib.request

    dest_dir.mkdir(parents=True, exist_ok=True)
    code = pdb_code.strip().lower()
    raw_pdb = dest_dir / f"{code}.pdb"

    # fetch
    with urllib.request.urlopen(f"https://files.rcsb.org/view/{code}.pdb") as r, open(raw_pdb, "wb") as o:
        o.write(r.read())

    # optional chain filter
    chains = (chains or "").replace(" ", "")
    if chains:
        wanted = {c.upper() for c in chains.split(",") if c}
        filtered = dest_dir / f"{code}.pdb"
        with open(raw_pdb) as fin, open(filtered, "w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    if line[21].upper() in wanted:
                        fout.write(line)
                else:
                    fout.write(line)
        rec_pdb = filtered
    else:
        rec_pdb = raw_pdb

    # naive ligand extraction: collect non-standard HETATM lines
    std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
           "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
    nonstd = set()
    with open(rec_pdb) as f:
        for line in f:
            if line.startswith("HETATM"):
                rn = line[17:20].strip().upper()
                if rn not in std:
                    nonstd.add(rn)

    ligand_sdf = None
    if nonstd:
        # write a simple ligand .pdb for later conversion (keep path around if you integrate OpenBabel/RDKit)
        lig_pdb = dest_dir / f"{code}_ligand.pdb"
        with open(rec_pdb) as fin, open(lig_pdb, "w") as fout:
            for line in fin:
                if line.startswith("HETATM"):
                    rn = line[17:20].strip().upper()
                    if rn in nonstd:
                        fout.write(line)
        # keep ligand_sdf=None for now; your confgen step can consume PDB or do conversion

    return {
        "pdb_path": str(rec_pdb),
        "receptor_pdb": str(rec_pdb),
        "ligand_sdf": ligand_sdf,
    }


# ---------- centers filename helpers ----------

def collect_receptor_tags(csv_path: Path) -> list[str]:
    tags: list[str] = []
    if not csv_path.exists():
        return tags
    import csv
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            base = Path(row["PDB_ID"]).stem
            t = base.replace("Receptors_", "")
            if t not in tags:
                tags.append(t)
    return tags


def rename_centers_with_tags(jobroot: Path):
    """If vina_centers.csv exists, copy -> vina_centers_<TAG1>_<TAG2>...csv (keep original too)."""
    root_csv = jobroot / "vina_centers.csv"
    if not root_csv.exists():
        return
    tags = collect_receptor_tags(root_csv)
    if not tags:
        return
    (jobroot / f"vina_centers_{'_'.join(tags)}.csv").write_text(root_csv.read_text())


# ---------- tiny embedded helpers (placeholders) ----------

PYMOL_CENTER_HELPER = '''from pymol import cmd
import csv, os
receptor_dir = "./Receptors"
output_csv = "vina_centers.csv"
current_receptor = None
remaining_receptors = []
if not os.path.exists(output_csv) or os.stat(output_csv).st_size == 0:
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f); writer.writerow(["PDB_ID","X","Y","Z","SIZE"])

def center_this(size=20):
    global current_receptor
    if not current_receptor: print("No receptor loaded"); return
    model = cmd.get_model("pocket")
    if len(model.atom) == 0:
        print("No atoms in selection 'pocket'"); return
    x=y=z=0.0
    for a in model.atom:
        x+=a.coord[0]; y+=a.coord[1]; z+=a.coord[2]
    n=len(model.atom); center=(x/n, y/n, z/n)
    with open(output_csv, "a", newline="") as f:
        csv.writer(f).writerow([f"{current_receptor}.pdbqt", round(center[0],3), round(center[1],3), round(center[2],3), size])
    print(f"Center for {current_receptor}: {center}")
    load_next_receptor()

def load_next_receptor():
    global current_receptor, remaining_receptors
    if not remaining_receptors:
        print("All done"); return
    fname = remaining_receptors.pop(0)
    obj = fname.replace(".pdbqt","")
    cmd.reinitialize(); cmd.load(os.path.join(receptor_dir, fname), obj)
    current_receptor = obj
    print("Loaded", fname, "— define selection 'pocket' then run center_this()")

receptors = [f for f in os.listdir(receptor_dir) if f.endswith('.pdbqt')]
processed=set()
if os.path.exists(output_csv):
    with open(output_csv) as f:
        for row in csv.DictReader(f): processed.add(row['PDB_ID'])
remaining_receptors = [f for f in receptors if f not in processed]
load_next_receptor()
'''

PDB2PDBQT_BATCH = '''# Placeholder for 3a_PDB2PDBQTbatch.py
# Implement your receptor cleanup + AutoDockTools calls here if desired.
'''
