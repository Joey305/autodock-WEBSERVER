#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
6_PyMacsPREP.py

Prepare MD-ready folder structure from PyMOL docking outputs.

INPUTS
------
- vina_centers.csv
- 11_AA_PDBOUTPUT/*.pdb

OUTPUT
------
MD_pymacs/
├── <RECEPTOR>/
│   ├── <COMPOUND_NAME>/
│   │   ├── complex.pdb
│   │   └── UNK.mol2   (hydrogenated, molecule name = UNK)
"""

import csv
import re
import shutil
import subprocess
from pathlib import Path

# ================= CONFIG =================

PDB_INPUT_DIR = Path("11_AA_PDBOUTPUT")
VINA_CSV      = Path("vina_centers.csv")
OUT_ROOT      = Path("MD_pymacs")

# Extract compound name between Top##_ and _state
COMPOUND_RE = re.compile(r"_Top\d{2}_(.+?)_state", re.IGNORECASE)

# ================= HELPERS =================

def load_receptors(vina_csv: Path) -> set[str]:
    receptors = set()
    with open(vina_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            receptors.add(row["PDB_ID"].replace(".pdbqt", "").strip())
    return receptors


def extract_compound_name(stem: str) -> str | None:
    """
    Extract compound name between Top##_ and _state
    """
    m = COMPOUND_RE.search(stem)
    return m.group(1) if m else None


def extract_unk_mol2(pdb_path: Path, out_mol2: Path):
    """
    Extract UNK residue from PDB and write hydrogenated MOL2
    with molecule name forced to 'UNK'.
    """
    tmp_pdb = out_mol2.with_suffix(".tmp.pdb")

    # --- Extract UNK atoms only ---
    with open(pdb_path) as fin, open(tmp_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "UNK":
                fout.write(line)
        fout.write("END\n")

    # --- Convert to MOL2 with hydrogens and correct title ---
    subprocess.run(
        [
            "obabel",
            str(tmp_pdb),
            "-O", str(out_mol2),
            "-h",
            "--title", "UNK"
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    tmp_pdb.unlink(missing_ok=True)


# ================= MAIN =================

def main():
    if not VINA_CSV.exists():
        raise SystemExit("❌ vina_centers.csv not found")

    if not PDB_INPUT_DIR.exists():
        raise SystemExit("❌ 11_AA_PDBOUTPUT directory not found")

    OUT_ROOT.mkdir(exist_ok=True)

    receptors = load_receptors(VINA_CSV)

    pdb_files = list(PDB_INPUT_DIR.glob("*.pdb"))
    if not pdb_files:
        raise SystemExit("❌ No PDB files found in 11_AA_PDBOUTPUT")

    for pdb in pdb_files:
        stem = pdb.stem

        # Identify receptor
        receptor = next((r for r in receptors if stem.startswith(r + "_")), None)
        if receptor is None:
            continue

        # Extract compound name
        compound = extract_compound_name(stem)
        if compound is None:
            continue

        # Create output directories
        rec_dir = OUT_ROOT / receptor
        cmp_dir = rec_dir / compound
        cmp_dir.mkdir(parents=True, exist_ok=True)

        # Copy and rename complex
        complex_pdb = cmp_dir / "complex.pdb"
        shutil.copy2(pdb, complex_pdb)

        # Generate UNK.mol2
        mol2_path = cmp_dir / "UNK.mol2"
        try:
            extract_unk_mol2(complex_pdb, mol2_path)
        except Exception as e:
            print(f"⚠️ Failed UNK extraction for {pdb.name}: {e}")

    print("\n✅ MD_pymacs preparation complete.")
    print(f"📁 Output directory: {OUT_ROOT.resolve()}")


if __name__ == "__main__":
    main()
