#!/usr/bin/env python3
import os, sys, csv, shutil, subprocess
from pathlib import Path
from datetime import datetime
from multiprocessing import Pool, cpu_count
import argparse
import re

from rdkit import Chem
from rdkit.Chem import AllChem

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
    p.add_argument("--keep-tmp", action="store_true",
                   help="Do NOT delete the TMP SDF folder at the end (default: delete to save space)")

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

def embed_and_optimize(mol, num_confs: int) -> list[int]:
    try:
        params = AllChem.ETKDGv3()
    except Exception:
        params = AllChem.ETKDG()
    params.numThreads = 0
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    for cid in conf_ids:
        try:
            AllChem.UFFOptimizeMolecule(mol, confId=cid)
        except Exception:
            pass
    return list(conf_ids)

# ----------------- workers -----------------
def generate_poses(task):
    """
    task: (ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL)
    """
    ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL = task

    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if mol is None:
        return 0

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass

    mol = Chem.AddHs(mol, addCoords=True)
    conf_ids = embed_and_optimize(mol, num_confs)

    written = 0
    for i, cid in enumerate(conf_ids, start=1):
        sdf_path   = tmp_dir  / f"{ligand_id}_pose{i}.sdf"
        pdbqt_path = pdbqt_dir / f"{ligand_id}__pose{i}.pdbqt"

        # write SDF for that conformer
        Chem.MolToMolFile(mol, str(sdf_path), confId=cid)

        # SDF -> PDBQT (Gasteiger charges)
        subprocess.run(
            [OBABEL, "-isdf", str(sdf_path), "-O", str(pdbqt_path), "--partialcharge", "gasteiger"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        print(f"✅ {ligand_id} pose {i} → {pdbqt_path}")
        written += 1

    return written

# ----------------- input collectors -----------------
def iter_mols_from_sdf_file(sdf_path: Path):
    """Yield (ligand_id, molblock) for every record in a multi-record SDF."""
    supp = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
    if not supp:
        return
    base = sanitize_id(sdf_path.stem)
    for idx, mol in enumerate(supp, start=1):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"{base}_rec{idx}"
        ligand_id = sanitize_id(name)
        yield ligand_id, Chem.MolToMolBlock(mol)

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
def run_one_folder(folder: Path, ft: str, num_confs: int, num_workers: int, OBABEL: str, keep_tmp: bool):
    print(f"\n=== Processing folder: {folder} (type={ft}) ===")

    tasks: list[tuple[str, str]] = []
    base_name = folder.name

    if ft == "sdf":
        for ligand_id, molblock in iter_mols_from_sdf_folder(folder):
            tasks.append((ligand_id, molblock))
    else:
        for ligand_id, molblock in iter_mols_from_smiles_folder(folder):
            tasks.append((ligand_id, molblock))

    if not tasks:
        print(f"⚠️  No valid molecules found in {folder}. Skipping.")
        return 0, None, None

    tag = f"{num_confs}Poses"
    pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

    print(f"🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
    print(f"📦 Output PDBQT: {pdbqt_dir}")
    print(f"🗂  TMP SDF:     {tmp_dir} (will {'NOT ' if keep_tmp else ''}be removed at end)")

    work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]

    written_total = 0
    try:
        with Pool(processes=num_workers) as pool:
            for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                written_total += int(n or 0)
    finally:
        if not keep_tmp:
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
        print(f"🗂  TMP SDF:     {tmp_dir} (will {'NOT ' if args.keep_tmp else ''}be removed at end)")

        work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]

        written_total = 0
        try:
            with Pool(processes=num_workers) as pool:
                for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                    written_total += int(n or 0)
        finally:
            if not args.keep_tmp:
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
            written, pdbqt_dir, tmp_dir = run_one_folder(folder, ft, num_confs, num_workers, OBABEL, args.keep_tmp)
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
        print(f"🗂  TMP SDF:     {tmp_dir} (will {'NOT ' if args.keep_tmp else ''}be removed at end)")

        work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]

        written_total = 0
        try:
            with Pool(processes=num_workers) as pool:
                for n in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                    written_total += int(n or 0)
        finally:
            if not args.keep_tmp:
                try:
                    shutil.rmtree(tmp_dir, ignore_errors=True)
                    print(f"🧹 TMP removed: {tmp_dir}")
                except Exception as e:
                    print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

        print(f"\n✅ Done! {written_total} poses written to: {pdbqt_dir}")

    else:
        print("❌ Invalid selection."); sys.exit(1)

