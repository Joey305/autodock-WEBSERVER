#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

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
import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from ligand_naming import parse_ligand_variant

# ================= CONFIG =================

PDB_INPUT_DIR = Path("11_AA_PDBOUTPUT")
VINA_CSV      = Path("vina_centers.csv")
OUT_ROOT      = Path("MD_pymacs")
PYMACS_REPO_URL = "https://github.com/schurerlab/pymacs.git"
PYMACS_SKIP_PATTERNS = (
    ".git",
    "Example_Choices",
)
PYMACS_HEAVY_SUFFIXES = (
    ".xtc",
    ".trr",
    ".tpr",
    ".cpt",
    ".edr",
    ".xvg",
    ".log",
    ".pdf",
    ".zip",
    ".tar",
    ".tar.gz",
)

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare MD-ready PyMACS folders from PyMOL docking outputs."
    )
    parser.add_argument(
        "--skip-pymacs-install",
        action="store_true",
        help="Skip copying the lightweight PyMACS scaffold into each generated folder.",
    )
    parser.add_argument(
        "--overwrite-pymacs",
        action="store_true",
        help="Overwrite existing PyMACS files/directories inside generated folders.",
    )
    return parser.parse_args()


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


def ensure_git_available():
    if shutil.which("git") is None:
        raise SystemExit(
            "❌ git is required to install the lightweight PyMACS scaffold. "
            "Install git or rerun with --skip-pymacs-install."
        )


def prune_heavy_pymacs_assets(repo_root: Path):
    example_dir = repo_root / "Example_Choices"
    if example_dir.exists():
        shutil.rmtree(example_dir)

    for path in repo_root.rglob("*"):
        if not path.is_file():
            continue
        name = path.name
        if any(name.endswith(suffix) for suffix in PYMACS_HEAVY_SUFFIXES):
            path.unlink()


def clone_lightweight_pymacs() -> tuple[tempfile.TemporaryDirectory, Path]:
    ensure_git_available()
    tmp_dir = tempfile.TemporaryDirectory()
    repo_dir = Path(tmp_dir.name) / "pymacs"
    subprocess.run(
        ["git", "clone", "--depth", "1", PYMACS_REPO_URL, str(repo_dir)],
        check=True,
        env={**os.environ, "GIT_LFS_SKIP_SMUDGE": "1"},
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    prune_heavy_pymacs_assets(repo_dir)
    return tmp_dir, repo_dir


def copy_lightweight_pymacs(source_dir: Path, target_dir: Path, overwrite: bool = False) -> list[str]:
    copied = []
    for item in source_dir.iterdir():
        if item.name in PYMACS_SKIP_PATTERNS:
            continue
        destination = target_dir / item.name
        if destination.exists():
            if not overwrite:
                continue
            if destination.is_dir() and not destination.is_symlink():
                shutil.rmtree(destination)
            else:
                destination.unlink()
        if item.is_dir():
            shutil.copytree(item, destination)
        else:
            shutil.copy2(item, destination)
        copied.append(item.name)
    return copied


# ================= MAIN =================

def main():
    args = parse_args()

    if not VINA_CSV.exists():
        raise SystemExit("❌ vina_centers.csv not found")

    if not PDB_INPUT_DIR.exists():
        raise SystemExit("❌ 11_AA_PDBOUTPUT directory not found")

    OUT_ROOT.mkdir(exist_ok=True)

    receptors = load_receptors(VINA_CSV)
    pymacs_tmp = None
    pymacs_source = None
    if not args.skip_pymacs_install:
        try:
            pymacs_tmp, pymacs_source = clone_lightweight_pymacs()
            print("📦 Prepared lightweight PyMACS scaffold for generated MD folders.")
        except subprocess.CalledProcessError as exc:
            stderr = (exc.stderr or "").strip()
            raise SystemExit(
                "❌ Failed to clone the PyMACS repository for folder bootstrap.\n"
                f"{stderr}\n"
                "If you want to generate MD folders without installing PyMACS, rerun with --skip-pymacs-install."
            )

    try:
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
            compound_variant = extract_compound_name(stem)
            if compound_variant is None:
                continue
            meta = parse_ligand_variant(compound_variant)
            compound_base = meta["LigandBase"] or compound_variant

            # Create output directories
            rec_dir = OUT_ROOT / receptor
            cmp_dir = rec_dir / compound_base
            if compound_variant != compound_base:
                cmp_dir = cmp_dir / compound_variant
            cmp_dir.mkdir(parents=True, exist_ok=True)

            # Copy and rename complex
            complex_pdb = cmp_dir / "complex.pdb"
            shutil.copy2(pdb, complex_pdb)

            if pymacs_source is not None:
                copied = copy_lightweight_pymacs(
                    pymacs_source,
                    cmp_dir,
                    overwrite=args.overwrite_pymacs,
                )
                if copied:
                    print(f"📥 Installed lightweight PyMACS scaffold into {cmp_dir}")
                else:
                    print(f"ℹ️ PyMACS scaffold already present in {cmp_dir}; nothing overwritten.")

            # Generate UNK.mol2
            mol2_path = cmp_dir / "UNK.mol2"
            try:
                extract_unk_mol2(complex_pdb, mol2_path)
            except Exception as e:
                print(f"⚠️ Failed UNK extraction for {pdb.name}: {e}")

        print("\n✅ MD_pymacs preparation complete.")
        print(f"📁 Output directory: {OUT_ROOT.resolve()}")
    finally:
        if pymacs_tmp is not None:
            pymacs_tmp.cleanup()


if __name__ == "__main__":
    main()
