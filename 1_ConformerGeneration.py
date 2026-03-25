
#!/usr/bin/env python3
import os, sys, csv, shutil, subprocess
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool, cpu_count
import argparse
import re
import numpy as np
from rdkit import Chem
from pathlib import Path
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import Chem
from pathlib import Path
from rdkit.Chem.MolStandardize import rdMolStandardize

# ----------------- CLI -----------------
def build_args():
    p = argparse.ArgumentParser(
        description="Generate conformers → PDBQT with RDKit + OpenBabel. Interactive by default; flags enable headless."
    )
    # Mode selection
    p.add_argument("--mode", choices=["1","2","3"],
                   help="1=CSV (SMILES+ID), 2=Folder(s) of SDF/SMILES, 3=Single multi-record SDF")

    # Common knobs
    p.add_argument("--num-confs", type=int, help="Poses per molecule (default 64)")
    p.add_argument("--workers", type=int, help="CPU cores to use (default=min(8, cpu_count))")
    p.add_argument("--obabel-bin", help="Path to obabel (overrides PATH detection)")
    p.add_argument("--remove-tmp", action="store_true",
               help="Delete the TMP SDF folder at the end (default: keep TMP SDFs)")


    # Mode 1 (CSV)
    p.add_argument("--csv", help="Path to CSV (mode=1)")
    p.add_argument("--smiles-col", dest="smiles_col", help="CSV column name for SMILES")
    p.add_argument("--id-col", dest="id_col", help="CSV column name for ligand ID")

    # Mode 2 (Folder) — now supports multiple
    p.add_argument("--folder", nargs="+",
                   help="One or more folders containing .sdf or .smiles files (mode=2). "
                        "Space-separated or a single comma-separated token is accepted.")
    p.add_argument("--filetype", choices=["sdf","smiles"], help="Folder filetype for mode=2")

    # Mode 3 (SDF)
    p.add_argument("--sdf", help="Path to a single multi-record .sdf (mode=3)")

    return p.parse_args()

# ----------------- utils -----------------
def obabel_path(cli_path: str | None):
    p = cli_path or os.environ.get("OBABEL_BIN") or shutil.which("obabel")
    if not p:
        print("❌ Cannot find 'obabel'. Install via conda-forge:\n"
              "   conda install -c conda-forge openbabel\n"
              "or set OBABEL_BIN or use --obabel-bin.", file=sys.stderr)
        sys.exit(2)
    return p



def sanitize_id(name: str) -> str:
    safe = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in (name or "ligand"))
    return safe.strip("_") or "ligand"

def create_output_dirs(base_name: str, tag: str) -> tuple[Path, Path]:
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    pdbqt_dir = Path(f"{base_name}_Ligands_PDBQT_{tag}_{ts}")
    tmp_dir   = Path(f"{base_name}_TMP_SDF_{tag}_{ts}")
    pdbqt_dir.mkdir(exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)
    return pdbqt_dir, tmp_dir





def embed_and_optimize(mol: Chem.Mol, num_confs: int) -> tuple[Chem.Mol, list[int]]:
    """
    Robust conformer generation:
    - ETKDG embedding (1-by-1 for stability)
    - UFF/MMFF optimization (best-effort)
    - NO conformer object reuse (vector-safe)
    - Multiprocessing safe
    - Returns conformer IDs on the FINAL mol
    """

    # --- Keep largest connected fragment only (THIS CREATES A NEW MOL) ---
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    mol = max(frags, key=lambda m: m.GetNumAtoms())

    # --- Ensure hydrogens ---
    mol = Chem.AddHs(mol, addCoords=True)

    # --- ETKDG parameters ---
    try:
        params = AllChem.ETKDGv3()
    except Exception:
        params = AllChem.ETKDG()

    params.useRandomCoords = True
    params.numThreads = 0  # multiprocessing-safe

    accepted_conf_ids = []
    attempts = 0
    max_attempts = num_confs * 5

    while len(accepted_conf_ids) < num_confs and attempts < max_attempts:
        attempts += 1

        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=1, params=params)

        for cid in conf_ids:
            # Optimize best-effort
            try:
                try:
                    AllChem.UFFOptimizeMolecule(mol, confId=cid)
                except Exception:
                    AllChem.MMFFOptimizeMolecule(mol, confId=cid)
            except Exception:
                pass

            # Geometry sanity check
            conf = mol.GetConformer(cid)
            coords = np.array([[conf.GetAtomPosition(i).x,
                                conf.GetAtomPosition(i).y,
                                conf.GetAtomPosition(i).z]
                               for i in range(mol.GetNumAtoms())])
            centroid = coords.mean(axis=0)
            dists = np.linalg.norm(coords - centroid, axis=1)

            if np.all(dists < 50.0):
                accepted_conf_ids.append(cid)
            else:
                mol.RemoveConformer(cid)

            if len(accepted_conf_ids) >= num_confs:
                break

    # Repack conformers to clean sequential IDs
    clean = Chem.Mol(mol)
    clean.RemoveAllConformers()

    for cid in accepted_conf_ids[:num_confs]:
        conf = mol.GetConformer(cid)
        clean.AddConformer(Chem.Conformer(conf), assignId=True)

    return clean, [c.GetId() for c in clean.GetConformers()]






def pick_canonical_tautomer(mol: Chem.Mol,
                            max_tautomers: int = 32,
                            max_transforms: int = 200) -> Chem.Mol:
    """
    Return ONE tautomer (canonical) to avoid:
      - exploding enumeration (1000 tautomers)
      - multiplying conformers per state
      - breaking downstream file naming

    Uses RDKit MolStandardize TautomerEnumerator with strict caps.
    """
    if mol is None:
        return None

    # Keep only largest fragment (consistent with your conformer generator)
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    except Exception:
        pass

    # Cleanup (best effort)
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.Reionizer().reionize(mol)
    except Exception:
        pass

    try:
        te = rdMolStandardize.TautomerEnumerator()
        # These two lines are the key to preventing the "1000 tautomers" blow-up
        te.SetMaxTautomers(int(max_tautomers))
        te.SetMaxTransforms(int(max_transforms))

        # Canonicalize picks one deterministically
        canon = te.Canonicalize(mol)
        return canon if canon is not None else mol
    except Exception:
        return mol







def is_connected(mol: Chem.Mol) -> bool:
    """Return True if the molecule is a single fragment (no disconnected atoms)."""
    if mol is None:
        return False
    frags = Chem.GetMolFrags(mol, asMols=False)
    return len(frags) == 1


def crest_conformers(ligand_id: str, sdf_path: Path, tmp_dir: Path, num_confs: int) -> list[Path]:
    """
    Run CREST on a molecule if RDKit conformer generation fails.
    Input: single conformer SDF
    Output: list of conformer SDF files
    """
    workdir = tmp_dir / f"{ligand_id}_crest"
    workdir.mkdir(exist_ok=True)

    # --- Ensure we have a valid RDKit mol ---
    mol = Chem.MolFromMolFile(str(sdf_path), sanitize=False, removeHs=False)
    if mol is None:
        print(f"❌ Could not load SDF for {ligand_id}")
        return []

    if not mol.GetNumConformers():
        # Generate 3D if missing
        mol = Chem.AddHs(mol, addCoords=True)
        try:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        except Exception:
            print(f"⚠️ {ligand_id}: failed to embed before CREST")

    # --- Write RDKit XYZ ---
    xyz_in = workdir / f"{ligand_id}.xyz"
    conf = mol.GetConformer()
    with open(xyz_in, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():<3} {pos.x: .6f} {pos.y: .6f} {pos.z: .6f}\n")

    if not xyz_in.exists() or xyz_in.stat().st_size < 10:
        print(f"❌ Failed to write XYZ for {ligand_id}")
        return []

    # --- Run CREST ---
    cmd = ["crest", str(xyz_in), "--quick", "--nconf", str(num_confs)]
    try:
        subprocess.run(cmd, cwd=workdir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ CREST crashed for {ligand_id}: {e}")
        return []

    # --- Collect CREST output ---
    xyz_out = workdir / "crest_conformers.xyz"
    if not xyz_out.exists() or xyz_out.stat().st_size < 10:
        print(f"❌ CREST produced no conformers for {ligand_id}")
        return []

    # --- Split XYZ into SDFs via OpenBabel ---
    try:
        subprocess.run(
            ["obabel", "-ixyz", str(xyz_out), "-osdf",
             "-O", str(workdir / f"{ligand_id}_pose.sdf"), "-m"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError:
        print(f"❌ Failed to convert CREST XYZ → SDF for {ligand_id}")
        return []

    return list(workdir.glob(f"{ligand_id}_pose*.sdf"))




def rebuild_from_smiles(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Rebuild a molecule via SMILES → 3D ETKDG.
    Useful if direct parsing or force field fails.
    """
    try:
        smiles = Chem.MolToSmiles(mol)
        m2 = Chem.MolFromSmiles(smiles)
        if not m2:
            return None
        m2 = Chem.AddHs(m2)
        AllChem.EmbedMolecule(m2, AllChem.ETKDG())
        return m2
    except Exception:
        return None



# ----------------- workers -----------------

def generate_poses(task):
    """
    task: (ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL)

    Output contract (matches your older script):
      TMP:  tmp_dir/<ligand_id>/<ligand_id>_pose{i}.sdf   (kept unless --remove-tmp)
      OUT:  pdbqt_dir/<ligand_id>_pose{i}.pdbqt           (top-level)

    Adds tautomer generation WITHOUT exploding state counts:
      - chooses ONE canonical tautomer
      - generates num_confs conformers for that tautomer
      - converts EVERY SDF → MOL2 → PDBQT
    """
    ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL = task

    lig_tmp_dir = tmp_dir / ligand_id
    lig_tmp_dir.mkdir(exist_ok=True, parents=True)

    # 1) Load molblock
    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if mol is None:
        print(f"❌ {ligand_id}: could not parse molblock")
        return 0

    # 2) Rebuild via SMILES for consistency (same as your old script)
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"❌ {ligand_id}: could not rebuild from SMILES")
            return 0
    except Exception:
        print(f"❌ {ligand_id}: SMILES rebuild failed")
        return 0

    # 3) Pick ONE tautomer (canonical) — prevents 1000-tautomer slowdown
    mol = pick_canonical_tautomer(mol, max_tautomers=32, max_transforms=200)
    if mol is None:
        print(f"❌ {ligand_id}: tautomer selection failed")
        return 0

    # 4) Add Hs then generate conformers
    try:
        mol, conf_ids = embed_and_optimize(mol, num_confs)
    except Exception as e:
        print(f"❌ {ligand_id}: embed/opt failed: {e}")
        return 0


    if not conf_ids:
        print(f"❌ {ligand_id}: no conformers generated")
        return 0

    written = 0
    total = len(conf_ids)

    for i, cid in enumerate(conf_ids, start=1):
        sdf_path  = lig_tmp_dir / f"{ligand_id}_pose{i}.sdf"
        mol2_path = lig_tmp_dir / f"{ligand_id}_pose{i}.mol2"
        pdbqt_path = pdbqt_dir / f"{ligand_id}_pose{i}.pdbqt"

        # Write SDF pose
        try:
            Chem.MolToMolFile(mol, str(sdf_path), confId=int(cid))
        except Exception as e:
            print(f"❌ {ligand_id}: failed writing SDF pose {i} (confId={cid}): {e}")
            continue


        # Convert SDF → MOL2 (with gasteiger)
        try:
            subprocess.run(
                [OBABEL, "-isdf", str(sdf_path),
                 "-omol2", "-O", str(mol2_path),
                 "--partialcharge", "gasteiger"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError:
            print(f"❌ OpenBabel failed (SDF→MOL2) for {ligand_id} pose {i}")
            continue

        # Convert MOL2 → PDBQT
        try:
            subprocess.run(
                [OBABEL, "-imol2", str(mol2_path),
                 "-opdbqt", "-O", str(pdbqt_path)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            written += 1
        except subprocess.CalledProcessError:
            print(f"❌ OpenBabel failed (MOL2→PDBQT) for {ligand_id} pose {i}")
            continue
        finally:
            # Keep SDF in TMP (by design). MOL2 is intermediate: delete it.
            try:
                if mol2_path.exists():
                    mol2_path.unlink()
            except Exception:
                pass

    print(f"[OK] {ligand_id}: {written}/{total} poses written → {pdbqt_dir}")
    return written




def has_exotic_atoms(mol: Chem.Mol) -> bool:
    """
    Return True if the molecule contains atoms that UFF/MMFF
    cannot handle well (e.g. hypervalent S, P, metals).
    """
    exotic = {"S", "P", "Se", "Te", "Si", "B", "Al", "Fe", "Zn", "Mg", "Na", "K", "Ca"}
    return any(atom.GetSymbol() in exotic for atom in mol.GetAtoms())



def needs_embedding(mol: Chem.Mol) -> bool:
    """Return True if mol has no conformers or all coords are (0,0,0)."""
    if mol is None or mol.GetNumAtoms() == 0:
        return True
    if mol.GetNumConformers() == 0:
        return True
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    if np.allclose(coords, 0.0):
        return True
    return False





def try_rdkit_load(sdf_path, sanitize=True):
    """Try loading with RDKit. Returns mol or None."""
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=sanitize, removeHs=False)
        mols = [m for m in suppl if m is not None]
        return mols if mols else None
    except Exception as e:
        return None



def mol_from_sdf_as_smiles(sdf_path: Path):
    """
    Yield (ligand_id, molblock) for each record in an SDF,
    but always rebuild via SMILES → 3D for consistency.
    """
    supp = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=True)
    if not supp:
        return

    base = sdf_path.stem
    seen = {}

    for idx, mol in enumerate(supp, start=1):
        if mol is None:
            continue

        # Canonical isomeric SMILES
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

        # Rebuild clean mol from SMILES
        m2 = Chem.MolFromSmiles(smiles)
        if not m2:
            continue
        m2 = Chem.AddHs(m2)

        # Decide ligand_id
        raw_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        if not raw_name:
            raw_name = base if len(supp) == 1 else f"{base}_rec{idx}"

        ligand_id = "".join(c if c.isalnum() or c in ("-", "_", ".") else "_" for c in raw_name)
        if ligand_id in seen:
            seen[ligand_id] += 1
            ligand_id = f"{ligand_id}_{seen[ligand_id]:02d}"
        else:
            seen[ligand_id] = 1

        yield ligand_id, Chem.MolToMolBlock(m2)


def mol_from_smiles_text(smiles: str, ligand_id: str):
    """
    Convert SMILES text directly into a hydrogenated MolBlock.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mol = Chem.AddHs(mol)
    return ligand_id, Chem.MolToMolBlock(mol)


# ----------------- input collectors -----------------
def iter_mols_from_sdf_file(sdf_path, obabel_bin="obabel"):
    """
    Robust iterator over molecules from a .sdf file with multi-step fallback:
    1. RDKit sanitize=True
    2. RDKit sanitize=False
    3. OpenBabel repair
    4. SMILES roundtrip → regenerate SDF → retry RDKit

    Each record tries to use the molecule title (mol.GetProp("_Name")),
    falling back to file stem + index if missing.
    Yields (ligand_id, molblock).
    """
    sdf_path = Path(sdf_path)
    base = sdf_path.stem

    def _yield_mols(mols):
        for idx, m in enumerate(mols, start=1):
            if m is None:
                continue
            # use title if present, otherwise fallback
            if m.HasProp("_Name") and m.GetProp("_Name").strip():
                ligand_id = sanitize_id(m.GetProp("_Name").strip())
            else:
                ligand_id = sanitize_id(f"{base}_{idx}")
            yield ligand_id, Chem.MolToMolBlock(m)

    # --- Step 1: RDKit sanitize=True ---
    mols = try_rdkit_load(sdf_path, sanitize=True)
    if mols:
        yield from _yield_mols(mols)
        return

    print(f"⚠️ {sdf_path.name}: RDKit failed (sanitize=True) → retry sanitize=False")
    mols = try_rdkit_load(sdf_path, sanitize=False)
    if mols:
        yield from _yield_mols(mols)
        return

    # --- Step 2: OpenBabel repair ---
    repaired = str(sdf_path.with_suffix(".repaired.sdf"))
    cmd = [obabel_bin, "-isdf", str(sdf_path), "-osdf", "-O", repaired, "--gen3d"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        mols = try_rdkit_load(repaired, sanitize=False)
        if mols:
            print(f"🔧 {sdf_path.name} repaired via OpenBabel")
            yield from _yield_mols(mols)
            return
    except Exception:
        pass

    # --- Step 3: SMILES roundtrip ---
    smiles_tmp = str(sdf_path.with_suffix(".smi"))
    sdf_regen  = str(sdf_path.with_suffix(".regen.sdf"))
    try:
        subprocess.run([obabel_bin, "-isdf", str(sdf_path), "-osmi", "-O", smiles_tmp],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run([obabel_bin, "-ismi", smiles_tmp, "-osdf", "-O", sdf_regen, "--gen3d"],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        mols = try_rdkit_load(sdf_regen, sanitize=False)
        if mols:
            print(f"♻️ {sdf_path.name} salvaged via SMILES roundtrip")
            yield from _yield_mols(mols)
            return
    except Exception:
        pass

    # --- Step 4: give up ---
    print(f"❌ {sdf_path.name}: could not parse after all attempts.")
    return







def iter_mols_from_sdf_folder(folder: Path):
    for sdf in sorted(folder.glob("*.sdf")):
        yield from iter_mols_from_sdf_file(sdf)

def iter_mols_from_smiles_folder(folder: Path):
    for smi in sorted(folder.glob("*.smiles")):
        try:
            smiles = smi.read_text().strip()
        except Exception:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            ligand_id = sanitize_id(smi.stem)
            mol = Chem.AddHs(mol)
            yield ligand_id, Chem.MolToMolBlock(mol)

# ----------------- interactive helpers -----------------
def choose(prompt, options):
    print(f"\n{prompt}")
    for i, opt in enumerate(options):
        print(f" [{i}] {opt}")
    while True:
        try:
            idx = int(input("Enter index: ").strip())
            if 0 <= idx < len(options):
                return options[idx]
        except Exception:
            pass
        print("Invalid choice—try again.")

def choose_multi(prompt, options):
    print(f"\n{prompt}")
    for i, opt in enumerate(options):
        print(f" [{i}] {opt}")
    while True:
        raw = input("Enter indices (comma-separated; ranges OK e.g. 0,2,5-7): ").strip()
        picks = parse_index_list(raw, len(options))
        if picks:
            return [options[i] for i in picks]
        print("Invalid selection—try again.")

def parse_index_list(s: str, n: int) -> list[int]:
    s = s.strip()
    if not s: return []
    picks = set()
    for part in s.split(","):
        part = part.strip()
        if not part: continue
        if "-" in part:
            a, b = part.split("-", 1)
            if a.strip().isdigit() and b.strip().isdigit():
                lo, hi = int(a), int(b)
                if lo > hi: lo, hi = hi, lo
                for k in range(lo, hi+1):
                    if 0 <= k < n: picks.add(k)
        else:
            if part.isdigit():
                k = int(part)
                if 0 <= k < n: picks.add(k)
    return sorted(picks)

def prompt_int(question, default):
    while True:
        s = input(f"{question} [{default}]: ").strip()
        if s == "":
            return default
        try:
            v = int(s)
            if v > 0:
                return v
        except Exception:
            pass
        print("Please enter a positive integer.")

# ----------------- per-folder runner -----------------
def run_one_folder(folder: Path, ft: str, num_confs: int, num_workers: int, OBABEL: str, remove_tmp: bool):
    print(f"\n=== Processing folder: {folder} (type={ft}) ===")

    tasks: list[tuple[str, str]] = []
    base_name = folder.name

    # Collect ligands from the folder
    if ft == "sdf":
        for ligand_id, molblock in iter_mols_from_sdf_folder(folder):
            tasks.append((ligand_id, molblock))
    else:
        for ligand_id, molblock in iter_mols_from_smiles_folder(folder):
            tasks.append((ligand_id, molblock))

    if not tasks:
        print(f"⚠️  No valid molecules found in {folder}. Skipping.")
        return 0, None, None

    # Create output dirs
    tag = f"{num_confs}Poses"
    pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

    print(f"🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
    print(f"📦 Output PDBQT: {pdbqt_dir}")
    print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if remove_tmp else 'be KEPT'})")

    # Build work items
    work_items = [
        (ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL)
        for (ligand_id, molblock) in tasks
    ]

    written_total = 0
    try:
        with Pool(processes=num_workers) as pool:
            for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                written_total += int(n or 0)
    finally:
        # 🧹 Only remove TMP if explicitly requested
        if remove_tmp:
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
                print(f"🧹 TMP removed: {tmp_dir}")
            except Exception as e:
                print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

    print(f"✅ Done! {written_total} poses written to: {pdbqt_dir}")
    return written_total, pdbqt_dir, tmp_dir


# ----------------- main -----------------
if __name__ == "__main__":
    args = build_args()
    mode = args.mode

    def _or_prompt(val, fn):
        return val if val is not None else fn()

    # Determine workers / poses
    total_cores = cpu_count() or 1
    default_workers = min(8, total_cores)

    num_confs = _or_prompt(
        args.num_confs,
        lambda: prompt_int("How many poses per molecule?", 64)
    )

    num_workers = _or_prompt(
        args.workers,
        lambda: prompt_int(f"How many CPU cores to use? (system has {total_cores})", default_workers)
    )

    # Resolve obabel path now (fail early if missing)
    OBABEL = obabel_path(args.obabel_bin)

    # Interactive fallback if mode not provided
    if mode is None:
        print("📘 Choose input type:")
        print(" [1] CSV file (SMILES + ID columns)")
        print(" [2] Folder of ligands (.sdf or .smiles)  ← supports MULTIPLE folders")
        print(" [3] Single .sdf file (multi-record)")
        mode = input("Enter 1, 2, or 3: ").strip()

    # --- Mode 1: CSV ---
    if mode == "1":
        # CSV path
        if args.csv:
            csv_path = Path(args.csv)
            if not csv_path.is_file():
                print(f"❌ CSV not found: {csv_path}"); sys.exit(1)
        else:
            csv_files = sorted(Path(".").glob("*.csv"), key=lambda p: p.name.lower())
            if not csv_files:
                print("❌ No CSV files found here."); sys.exit(1)
            pick = choose("📄 Select a CSV:", [p.name for p in csv_files])
            csv_path = next(p for p in csv_files if p.name == pick)

        base_name = csv_path.stem

        tasks: list[tuple[str, str]] = []
        with open(csv_path, "r", newline="") as f:
            r = csv.DictReader(f)
            headers = r.fieldnames or []
            if not headers:
                print("❌ CSV has no header."); sys.exit(1)

            if args.smiles_col and args.id_col:
                smi_col = args.smiles_col
                id_col  = args.id_col
            else:
                print("\n📊 CSV Columns:")
                for i, h in enumerate(headers):
                    print(f" [{i}] {h}")
                smi_col = headers[prompt_int("Index of SMILES column?", 0)]
                id_col  = headers[prompt_int("Index of ID column?", 1)]

            for row in r:
                smiles = (row.get(smi_col) or "").strip()
                ligand_id = sanitize_id((row.get(id_col) or "").strip())
                if not smiles or not ligand_id:
                    continue
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mol = Chem.AddHs(mol)
                    tasks.append((ligand_id, Chem.MolToMolBlock(mol)))

        if not tasks:
            print("⚠️ No valid molecules found."); sys.exit(2)

        tag = f"{num_confs}Poses"
        pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

        print(f"\n🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
        print(f"📦 Output PDBQT: {pdbqt_dir}")
        print(f"🗂  TMP SDF:     {tmp_dir} (will {'NOT ' if args.remove_tmp else ''}be removed at end)")

        work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]

        written_total = 0
        try:
            with Pool(processes=num_workers) as pool:
                for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                    written_total += int(n or 0)
        finally:
            print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if args.remove_tmp else 'be KEPT'})")
            if args.remove_tmp:
                shutil.rmtree(tmp_dir, ignore_errors=True)

                try:
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    print(f"🧹 TMP removed: {tmp_dir}")
                except Exception as e:
                    print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

        print(f"\n✅ Done! {written_total} poses written to: {pdbqt_dir}")

    # --- Mode 2: MULTIPLE folders (sdf OR smiles) ---
    elif mode == "2":
        # Resolve folder list (headless or interactive)
        folders: list[Path] = []
        if args.folder:
            # Support both space-separated and a single comma-separated token
            raw = []
            for tok in args.folder:
                raw.extend([t for t in re.split(r"\s*,\s*", tok) if t])
            for r in raw:
                p = Path(r)
                if not p.is_dir():
                    print(f"❌ Folder not found: {p}"); sys.exit(1)
                folders.append(p.resolve())
            # Sort headless selection deterministically
            folders = sorted(folders, key=lambda p: p.name.lower())
            ft = (args.filetype or "sdf").lower()
        else:
            subdirs = sorted([d for d in Path(".").iterdir() if d.is_dir()],
                            key=lambda d: d.name.lower())
            if not subdirs:
                print("❌ No subdirectories here."); sys.exit(1)
            multi = choose_multi("📁 Select one or more folders containing ligand files:",
                                [d.name for d in subdirs])
            folders = [next(d for d in subdirs if d.name == name) for name in multi]
            # Ensure deterministic order in case user typed ranges in mixed order
            folders = sorted(folders, key=lambda p: p.name.lower())
            ft = input("🔹 Selected folder(s) contain [sdf] or [smiles] files? ").strip().lower() or "sdf"

        if ft not in ("sdf","smiles"):
            print("❌ Invalid filetype for mode=2 (use sdf/smiles)."); sys.exit(1)

        # Process each folder sequentially with same settings
        total_written = 0
        outputs = []
        for folder in folders:
            written, pdbqt_dir, tmp_dir = run_one_folder(folder, ft, num_confs, num_workers, OBABEL, args.remove_tmp)
            total_written += written
            if pdbqt_dir is not None:
                outputs.append(pdbqt_dir)

        # Sort outputs for a clean summary print
        outputs = sorted(outputs, key=lambda p: str(p).lower())

        print("\n====== SUMMARY ======")
        print(f"Folders processed : {len(folders)}")
        print(f"Total poses written: {total_written}")
        if outputs:
            print("PDBQT outputs:")
            for o in outputs:
                print(f" - {o}")
        print("=====================\n")


    # --- Mode 3: Single multi-record SDF ---
    elif mode == "3":
        if args.sdf:
            sdf_path = Path(args.sdf)
            if not sdf_path.is_file():
                print(f"❌ SDF not found: {sdf_path}"); sys.exit(1)
        else:
            sdf_files = sorted(Path(".").glob("*.sdf"), key=lambda p: p.name.lower())
            if not sdf_files:
                print("❌ No .sdf files here."); sys.exit(1)
            pick = choose("🧪 Select the SDF file:", [p.name for p in sdf_files])
            sdf_path = next(p for p in sdf_files if p.name == pick)
        base_name = sdf_path.stem

        tasks: list[tuple[str, str]] = []
        for ligand_id, molblock in iter_mols_from_sdf_file(sdf_path):
            tasks.append((ligand_id, molblock))

        if not tasks:
            print("⚠️ No valid molecules found."); sys.exit(2)

        tag = f"{num_confs}Poses"
        pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

        print(f"\n🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
        print(f"📦 Output PDBQT: {pdbqt_dir}")
        print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if args.remove_tmp else 'be KEPT'})")


        work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]

        written_total = 0
        try:
            with Pool(processes=num_workers) as pool:
                for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                    written_total += int(n or 0)
        finally:
            print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if args.remove_tmp else 'be KEPT'})")
            if args.remove_tmp:
                shutil.rmtree(tmp_dir, ignore_errors=True)

                try:
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    print(f"🧹 TMP removed: {tmp_dir}")
                except Exception as e:
                    print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

        print(f"\n✅ Done! {written_total} poses written to: {pdbqt_dir}")

    else:
        print("❌ Invalid selection."); sys.exit(1)

