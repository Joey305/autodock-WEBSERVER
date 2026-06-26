
#!/usr/bin/env python3
from __future__ import annotations
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
from ligand_manifest import build_ligand_state_metadata, write_ligand_state_manifest
# =============================
# 1) ADD THESE IMPORTS
# =============================

try:
    from dimorphite_dl import protonate_smiles
    HAVE_DIMORPHITE = True
except Exception:
    protonate_smiles = None
    HAVE_DIMORPHITE = False


# ----------------- CLI -----------------
def build_args():
    p = argparse.ArgumentParser(
        description="Generate conformers → PDBQT with RDKit + OpenBabel. Interactive by default; flags enable headless."
    )
    # Mode selection
    p.add_argument("--mode", choices=["1","2","3","4"],
                   help="1=CSV (SMILES+ID), 2=Folder(s) of SDF/SMILES, 3=Single multi-record SDF, 4=Folder containing a CSV")

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
    # =============================
    # 2) ADD THESE CLI ARGS
    # inside build_args()
    # =============================

    p.add_argument("--enumerate-protomers", action="store_true",
                   help="Enumerate protonation/ionization states with Dimorphite-DL")
    p.add_argument("--ph-min", type=float, default=6.8,
                   help="Minimum pH for protonation-state enumeration")
    p.add_argument("--ph-max", type=float, default=7.4,
                   help="Maximum pH for protonation-state enumeration")
    p.add_argument("--ph-precision", type=float, default=0.5,
                   help="Dimorphite precision parameter")
    p.add_argument("--max-protomers", type=int, default=4,
                   help="Maximum protomer states to keep per ligand")
    p.add_argument("--max-tautomers", type=int, default=4,
                   help="Maximum tautomer states to keep per protomer")
    p.add_argument("--max-transforms", type=int, default=200,
                   help="RDKit tautomer transform cap")

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


def ligand_id_from_source_file(path: Path, record_index: int | None = None, total_records: int | None = None) -> str:
    base = sanitize_id(path.stem)
    if record_index is not None:
        width = max(3, len(str(total_records or record_index)))
        return f"{base}_rec{int(record_index):0{width}d}"
    return base


def ensure_unique_ligand_id(ligand_id: str, seen: dict[str, int], context: str = "") -> str:
    ligand_id = sanitize_id(ligand_id)
    count = seen.get(ligand_id, 0)
    if count == 0:
        seen[ligand_id] = 1
        return ligand_id
    seen[ligand_id] = count + 1
    unique = f"{ligand_id}_{seen[ligand_id]:03d}"
    where = f" from {context}" if context else ""
    print(f"⚠️ Ligand ID collision for {ligand_id}{where}; using {unique}")
    return unique


def _ignore_ligand_entry(name: str) -> bool:
    return (
        not name
        or name in {"__MACOSX", ".DS_Store"}
        or name.startswith(".")
        or name.startswith("._")
    )


def _ignore_ligand_dir(name: str) -> bool:
    return (
        _ignore_ligand_entry(name)
        or name.startswith("Ligands_TMP_SDF_")
        or (name.startswith("Ligands_") and "_PDBQT_" in name)
    )


def discover_ligand_input_files(folder: Path, suffixes: tuple[str, ...]) -> list[Path]:
    folder = Path(folder)
    found: list[Path] = []
    seen: set[Path] = set()
    nested_hits = 0

    for root, dirnames, filenames in os.walk(folder):
        root_path = Path(root)
        dirnames[:] = [name for name in dirnames if not _ignore_ligand_dir(name)]
        for filename in sorted(filenames):
            if _ignore_ligand_entry(filename):
                continue
            path = root_path / filename
            if path.suffix.lower() not in suffixes:
                continue
            try:
                resolved = path.resolve()
            except Exception:
                resolved = path
            if resolved in seen:
                continue
            seen.add(resolved)
            if root_path != folder:
                nested_hits += 1
            found.append(path)

    found.sort()
    if nested_hits:
        print(f"🔎 Detected nested ligand files under {folder}; using recursive discovery.")
    print(f"🔎 Discovered {len(found)} ligand input file(s) under {folder}")
    return found

def create_output_dirs(base_name: str, tag: str) -> tuple[Path, Path]:
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    pdbqt_dir = Path(f"{base_name}_Ligands_PDBQT_{tag}_{ts}")
    tmp_dir   = Path(f"{base_name}_TMP_SDF_{tag}_{ts}")
    pdbqt_dir.mkdir(exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)
    return pdbqt_dir, tmp_dir


def write_manifest_for_outputs(pdbqt_dir: Path, rows: list[dict]) -> Path | None:
    if not rows:
        return None
    manifest_path = Path(pdbqt_dir) / "ligand_state_manifest.csv"
    write_ligand_state_manifest(manifest_path, rows)
    print(f"🗒️  Wrote ligand state manifest: {manifest_path}")
    return manifest_path




def embed_and_optimize(mol: Chem.Mol, num_confs: int) -> tuple[Chem.Mol, list[int]]:
    """
    Generate genuinely distinct RDKit conformers for one tautomer/protomer state.

    Fixes the duplicate-conformer bug caused by repeatedly calling
    EmbedMultipleConfs(..., numConfs=1) and then repacking stale conformer IDs.
    """

    # Keep largest connected fragment only
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    mol = max(frags, key=lambda m: m.GetNumAtoms())

    # Ensure hydrogens
    mol = Chem.AddHs(mol, addCoords=True)

    # ETKDG parameters
    try:
        params = AllChem.ETKDGv3()
    except Exception:
        params = AllChem.ETKDG()

    params.useRandomCoords = True
    params.numThreads = 0

    # Optional but useful: avoid near-identical conformers.
    # Lower value = keep more conformers; higher value = stricter pruning.
    try:
        params.pruneRmsThresh = 0.25
    except Exception:
        pass

    # Generate all conformers in ONE call so conformer IDs remain valid.
    conf_ids = list(AllChem.EmbedMultipleConfs(
        mol,
        numConfs=int(num_confs),
        params=params
    ))

    accepted_conf_ids = []

    for cid in conf_ids:
        # Optimize best-effort
        try:
            try:
                AllChem.UFFOptimizeMolecule(mol, confId=int(cid))
            except Exception:
                AllChem.MMFFOptimizeMolecule(mol, confId=int(cid))
        except Exception:
            pass

        # Geometry sanity check
        try:
            conf = mol.GetConformer(int(cid))
            coords = np.array([
                [
                    conf.GetAtomPosition(i).x,
                    conf.GetAtomPosition(i).y,
                    conf.GetAtomPosition(i).z,
                ]
                for i in range(mol.GetNumAtoms())
            ])

            centroid = coords.mean(axis=0)
            dists = np.linalg.norm(coords - centroid, axis=1)

            if np.all(dists < 50.0):
                accepted_conf_ids.append(int(cid))
        except Exception:
            continue

    # Repack accepted conformers into a clean molecule with clean IDs
    clean = Chem.Mol(mol)
    clean.RemoveAllConformers()

    for cid in accepted_conf_ids[:num_confs]:
        conf = mol.GetConformer(int(cid))
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


# =============================
# 3) ADD THESE HELPERS
# place them near your other utils
# =============================

def keep_largest_fragment(mol: Chem.Mol) -> Chem.Mol | None:
    if mol is None:
        return None
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        return max(frags, key=lambda m: m.GetNumAtoms()) if frags else mol
    except Exception:
        return mol


def mol_key(mol: Chem.Mol) -> str:
    """Canonical key for deduplicating states."""
    try:
        m = Chem.Mol(mol)
        m = Chem.RemoveHs(m)
        return Chem.MolToSmiles(m, isomericSmiles=True)
    except Exception:
        try:
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        except Exception:
            return f"mol_{id(mol)}"


def standardize_preserve_state(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Standardize gently, but DO NOT reionize.
    Reionizer would tend to collapse your protomer diversity.
    """
    if mol is None:
        return None

    mol = keep_largest_fragment(mol)
    if mol is None:
        return None

    try:
        mol = rdMolStandardize.Cleanup(mol)
    except Exception:
        pass

    return mol


def enumerate_protomers(
    mol: Chem.Mol,
    enumerate_protomers_flag: bool = False,
    ph_min: float = 6.8,
    ph_max: float = 7.4,
    ph_precision: float = 0.5,
    max_protomers: int = 4,
) -> list[Chem.Mol]:
    """
    Return one or more protonation/ionization states.
    If enumeration is off, returns a single standardized state.
    """
    mol = standardize_preserve_state(mol)
    if mol is None:
        return []

    if not enumerate_protomers_flag:
        return [mol]

    if not HAVE_DIMORPHITE:
        raise RuntimeError(
            "--enumerate-protomers requested, but dimorphite_dl is not installed. "
            "Install with: pip install dimorphite_dl"
        )

    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)

    try:
        variants = protonate_smiles(
            smiles,
            ph_min=ph_min,
            ph_max=ph_max,
            precision=ph_precision,
            max_variants=max_protomers
        )
    except Exception as e:
        print(f"⚠️ Protomer enumeration failed; falling back to single state: {e}")
        return [mol]

    out = []
    seen = set()

    for item in variants:
        # Some versions may return strings; if labels are enabled in future,
        # handle tuples/lists defensively.
        psmi = item[0] if isinstance(item, (tuple, list)) else item
        pmol = Chem.MolFromSmiles(psmi)
        if pmol is None:
            continue

        pmol = standardize_preserve_state(pmol)
        if pmol is None:
            continue

        key = mol_key(pmol)
        if key not in seen:
            seen.add(key)
            out.append(pmol)

        if len(out) >= max_protomers:
            break

    return out or [mol]


def enumerate_tautomers(
    mol: Chem.Mol,
    max_tautomers: int = 4,
    max_transforms: int = 200,
) -> list[Chem.Mol]:
    """
    Enumerate tautomers for ONE protomer.
    """
    mol = standardize_preserve_state(mol)
    if mol is None:
        return []

    try:
        te = rdMolStandardize.TautomerEnumerator()
        te.SetMaxTautomers(int(max_tautomers))
        te.SetMaxTransforms(int(max_transforms))

        res = te.Enumerate(mol)
        candidates = list(res.tautomers)
    except Exception as e:
        print(f"⚠️ Tautomer enumeration failed; falling back to single tautomer: {e}")
        candidates = [mol]

    out = []
    seen = set()

    for tmol in candidates:
        tmol = standardize_preserve_state(tmol)
        if tmol is None:
            continue

        key = mol_key(tmol)
        if key not in seen:
            seen.add(key)
            out.append(tmol)

        if len(out) >= max_tautomers:
            break

    return out or [mol]




# ----------------- workers -----------------
# =============================
# 4) REPLACE generate_poses()
# =============================

def generate_poses(task):
    """
    task:
      (
        ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL, state_opts
      )

    Enumerates:
      protomer -> tautomer -> conformers

    Output names:
      ligand_p01_t01_c001.pdbqt
      ligand_p01_t01_c002.pdbqt
      ligand_p01_t02_c001.pdbqt
      ...
    """
    (
        ligand_id,
        molblock,
        tmp_dir,
        pdbqt_dir,
        num_confs,
        OBABEL,
        state_opts,
        source_input,
        source_input_type,
        source_ligand_id,
        source_record_index,
        source_record_name,
        original_mol_name,
    ) = task

    lig_tmp_dir = tmp_dir / ligand_id
    lig_tmp_dir.mkdir(exist_ok=True, parents=True)
    manifest_rows = []

    # 1) Load molblock
    mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
    if mol is None:
        print(f"❌ {ligand_id}: could not parse molblock")
        return {"written": 0, "rows": manifest_rows}

    # 2) Rebuild via SMILES for consistency
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"❌ {ligand_id}: could not rebuild from SMILES")
            return {"written": 0, "rows": manifest_rows}
    except Exception:
        print(f"❌ {ligand_id}: SMILES rebuild failed")
        return {"written": 0, "rows": manifest_rows}

    canonical_mol = Chem.Mol(mol)

    # 3) Enumerate protomers
    try:
        protomers = enumerate_protomers(
            mol,
            enumerate_protomers_flag=state_opts["enumerate_protomers"],
            ph_min=state_opts["ph_min"],
            ph_max=state_opts["ph_max"],
            ph_precision=state_opts["ph_precision"],
            max_protomers=state_opts["max_protomers"],
        )
    except Exception as e:
        print(f"❌ {ligand_id}: protomer enumeration failed: {e}")
        return {"written": 0, "rows": manifest_rows}

    if not protomers:
        print(f"❌ {ligand_id}: no protomers generated")
        return {"written": 0, "rows": manifest_rows}

    written = 0

    for p_idx, pmol in enumerate(protomers, start=1):
        # 4) Enumerate tautomers for each protomer
        tautomers = enumerate_tautomers(
            pmol,
            max_tautomers=state_opts["max_tautomers"],
            max_transforms=state_opts["max_transforms"],
        )

        if not tautomers:
            continue

        for t_idx, tmol in enumerate(tautomers, start=1):
            # 5) Conformer generation per tautomer
            try:
                conf_mol, conf_ids = embed_and_optimize(Chem.Mol(tmol), num_confs)
                print(f"[STATE] {ligand_id} p{p_idx:02d} t{t_idx:02d}: {len(conf_ids)}/{num_confs} conformers generated")
            except Exception as e:
                print(f"⚠️ {ligand_id} p{p_idx:02d} t{t_idx:02d}: embed/opt failed: {e}")
                continue

            if not conf_ids:
                continue

            state_prefix = f"{ligand_id}_p{p_idx:02d}_t{t_idx:02d}"
            state_tmp_dir = lig_tmp_dir / state_prefix
            state_tmp_dir.mkdir(exist_ok=True, parents=True)

            for c_idx, cid in enumerate(conf_ids, start=1):
                ligand_variant = f"{state_prefix}_c{c_idx:03d}"
                sdf_path = state_tmp_dir / f"{ligand_variant}.sdf"
                mol2_path = state_tmp_dir / f"{ligand_variant}.mol2"
                pdbqt_path = pdbqt_dir / f"{ligand_variant}.pdbqt"
                manifest_kwargs = {
                    "ligand_base": ligand_id,
                    "ligand_variant": ligand_variant,
                    "source_ligand_id": source_ligand_id or ligand_id,
                    "source_input": source_input,
                    "source_input_type": source_input_type,
                    "source_record_index": source_record_index,
                    "source_record_name": source_record_name,
                    "sdf_title": source_record_name,
                    "original_mol_name": original_mol_name,
                    "canonical_mol": canonical_mol,
                    "state_mol": tmol,
                    "pdbqt_file": pdbqt_path.name,
                    "sdf_file": sdf_path.name,
                    "tmp_sdf_file": str(sdf_path.relative_to(tmp_dir)),
                }

                # Write SDF conformer
                try:
                    Chem.MolToMolFile(conf_mol, str(sdf_path), confId=int(cid))
                except Exception as e:
                    print(f"❌ {state_prefix}: failed writing conformer {c_idx}: {e}")
                    manifest_rows.append(
                        build_ligand_state_metadata(
                            **manifest_kwargs,
                            generation_status="failed_sdf",
                            generation_warning=str(e),
                        )
                    )
                    continue

                # Convert SDF -> MOL2 with Gasteiger charges
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
                    print(f"❌ OpenBabel failed (SDF→MOL2) for {state_prefix} c{c_idx:03d}")
                    manifest_rows.append(
                        build_ligand_state_metadata(
                            **manifest_kwargs,
                            generation_status="failed_mol2",
                            generation_warning="OpenBabel failed during SDF to MOL2 conversion",
                        )
                    )
                    continue

                # Convert MOL2 -> PDBQT
                try:
                    subprocess.run(
                        [OBABEL, "-imol2", str(mol2_path),
                         "-opdbqt", "-O", str(pdbqt_path)],
                        check=True,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL
                    )
                    written += 1
                    manifest_rows.append(
                        build_ligand_state_metadata(
                            **manifest_kwargs,
                            generation_status="ok",
                            generation_warning="",
                        )
                    )
                except subprocess.CalledProcessError:
                    print(f"❌ OpenBabel failed (MOL2→PDBQT) for {state_prefix} c{c_idx:03d}")
                    manifest_rows.append(
                        build_ligand_state_metadata(
                            **manifest_kwargs,
                            generation_status="failed_pdbqt",
                            generation_warning="OpenBabel failed during MOL2 to PDBQT conversion",
                        )
                    )
                    continue
                finally:
                    try:
                        if mol2_path.exists():
                            mol2_path.unlink()
                    except Exception:
                        pass

    print(f"[OK] {ligand_id}: {written} total structures written → {pdbqt_dir}")
    return {"written": written, "rows": manifest_rows}


# def generate_poses(task):
#     """
#     task: (ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL)

#     Output contract (matches your older script):
#       TMP:  tmp_dir/<ligand_id>/<ligand_id>_pose{i}.sdf   (kept unless --remove-tmp)
#       OUT:  pdbqt_dir/<ligand_id>_pose{i}.pdbqt           (top-level)

#     Adds tautomer generation WITHOUT exploding state counts:
#       - chooses ONE canonical tautomer
#       - generates num_confs conformers for that tautomer
#       - converts EVERY SDF → MOL2 → PDBQT
#     """
#     ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL = task

#     lig_tmp_dir = tmp_dir / ligand_id
#     lig_tmp_dir.mkdir(exist_ok=True, parents=True)

#     # 1) Load molblock
#     mol = Chem.MolFromMolBlock(molblock, sanitize=True, removeHs=False)
#     if mol is None:
#         print(f"❌ {ligand_id}: could not parse molblock")
#         return 0

#     # 2) Rebuild via SMILES for consistency (same as your old script)
#     try:
#         smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
#         mol = Chem.MolFromSmiles(smiles)
#         if mol is None:
#             print(f"❌ {ligand_id}: could not rebuild from SMILES")
#             return 0
#     except Exception:
#         print(f"❌ {ligand_id}: SMILES rebuild failed")
#         return 0

#     # 3) Pick ONE tautomer (canonical) — prevents 1000-tautomer slowdown
#     mol = pick_canonical_tautomer(mol, max_tautomers=32, max_transforms=200)
#     if mol is None:
#         print(f"❌ {ligand_id}: tautomer selection failed")
#         return 0

#     # 4) Add Hs then generate conformers
#     try:
#         mol, conf_ids = embed_and_optimize(mol, num_confs)
#     except Exception as e:
#         print(f"❌ {ligand_id}: embed/opt failed: {e}")
#         return 0


#     if not conf_ids:
#         print(f"❌ {ligand_id}: no conformers generated")
#         return 0

#     written = 0
#     total = len(conf_ids)

#     for i, cid in enumerate(conf_ids, start=1):
#         sdf_path  = lig_tmp_dir / f"{ligand_id}_pose{i}.sdf"
#         mol2_path = lig_tmp_dir / f"{ligand_id}_pose{i}.mol2"
#         pdbqt_path = pdbqt_dir / f"{ligand_id}_pose{i}.pdbqt"

#         # Write SDF pose
#         try:
#             Chem.MolToMolFile(mol, str(sdf_path), confId=int(cid))
#         except Exception as e:
#             print(f"❌ {ligand_id}: failed writing SDF pose {i} (confId={cid}): {e}")
#             continue


#         # Convert SDF → MOL2 (with gasteiger)
#         try:
#             subprocess.run(
#                 [OBABEL, "-isdf", str(sdf_path),
#                  "-omol2", "-O", str(mol2_path),
#                  "--partialcharge", "gasteiger"],
#                 check=True,
#                 stdout=subprocess.DEVNULL,
#                 stderr=subprocess.DEVNULL
#             )
#         except subprocess.CalledProcessError:
#             print(f"❌ OpenBabel failed (SDF→MOL2) for {ligand_id} pose {i}")
#             continue

#         # Convert MOL2 → PDBQT
#         try:
#             subprocess.run(
#                 [OBABEL, "-imol2", str(mol2_path),
#                  "-opdbqt", "-O", str(pdbqt_path)],
#                 check=True,
#                 stdout=subprocess.DEVNULL,
#                 stderr=subprocess.DEVNULL
#             )
#             written += 1
#         except subprocess.CalledProcessError:
#             print(f"❌ OpenBabel failed (MOL2→PDBQT) for {ligand_id} pose {i}")
#             continue
#         finally:
#             # Keep SDF in TMP (by design). MOL2 is intermediate: delete it.
#             try:
#                 if mol2_path.exists():
#                     mol2_path.unlink()
#             except Exception:
#                 pass

#     print(f"[OK] {ligand_id}: {written}/{total} poses written → {pdbqt_dir}")
#     return written




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
    total_records = sum(1 for mol in supp if mol is not None)
    supp = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=True)

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

        record_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
        raw_id = ligand_id_from_source_file(sdf_path, idx if total_records > 1 else None, total_records)
        ligand_id = ensure_unique_ligand_id(raw_id, seen, context=f"{sdf_path.name} record {idx}")

        yield {
            "ligand_id": ligand_id,
            "molblock": Chem.MolToMolBlock(m2),
            "source_input": str(sdf_path),
            "source_input_type": "sdf",
            "source_ligand_id": ligand_id,
            "source_record_index": idx,
            "source_record_name": record_name,
            "original_mol_name": record_name,
        }


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

    Each record uses the source filename as the primary ligand identity.
    Internal SDF titles are preserved only as metadata.
    """
    sdf_path = Path(sdf_path)
    total_records_hint = 0
    try:
        raw_supp = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        total_records_hint = sum(1 for m in raw_supp if m is not None)
    except Exception:
        total_records_hint = 0

    def _yield_mols(mols):
        seen = {}
        for idx, m in enumerate(mols, start=1):
            if m is None:
                continue
            record_name = m.GetProp("_Name").strip() if m.HasProp("_Name") else ""
            use_record_index = idx if total_records_hint > 1 else None
            raw_id = ligand_id_from_source_file(sdf_path, use_record_index, total_records_hint or None)
            ligand_id = ensure_unique_ligand_id(raw_id, seen, context=f"{sdf_path.name} record {idx}")
            yield {
                "ligand_id": ligand_id,
                "molblock": Chem.MolToMolBlock(m),
                "source_input": str(sdf_path),
                "source_input_type": "sdf",
                "source_ligand_id": ligand_id,
                "source_record_index": idx,
                "source_record_name": record_name,
                "original_mol_name": record_name,
            }

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
    for sdf in discover_ligand_input_files(folder, (".sdf",)):
        yield from iter_mols_from_sdf_file(sdf)

def iter_mols_from_smiles_folder(folder: Path):
    seen = {}
    for smi in discover_ligand_input_files(folder, (".smiles", ".smi")):
        try:
            smiles = smi.read_text().strip()
        except Exception:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            ligand_id = ensure_unique_ligand_id(ligand_id_from_source_file(smi), seen, context=smi.name)
            mol = Chem.AddHs(mol)
            yield {
                "ligand_id": ligand_id,
                "molblock": Chem.MolToMolBlock(mol),
                "source_input": str(smi),
                "source_input_type": "smiles",
                "source_ligand_id": ligand_id,
                "source_record_index": 1,
                "source_record_name": "",
                "original_mol_name": "",
            }

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

def prompt_int(question, default, *, min_value=1):
    while True:
        s = input(f"{question} [{default}]: ").strip()
        if s == "":
            return default
        try:
            v = int(s)
            if v >= min_value:
                return v
        except Exception:
            pass
        print(f"Please enter an integer greater than or equal to {min_value}.")


def choose_csv_column(headers: list[str], question: str, default_index: int) -> str:
    return headers[prompt_int(question, default_index, min_value=0)]


def resolve_csv_path(args) -> Path:
    if args.csv:
        csv_path = Path(args.csv)
        if not csv_path.is_file():
            print(f"❌ CSV not found: {csv_path}")
            sys.exit(1)
        return csv_path

    csv_files = sorted(Path(".").glob("*.csv"), key=lambda p: p.name.lower())
    if not csv_files:
        print("❌ No CSV files found here.")
        sys.exit(1)
    pick = choose("📄 Select a CSV:", [p.name for p in csv_files])
    return next(p for p in csv_files if p.name == pick)


def resolve_csv_path_from_folder(args) -> Path:
    if args.folder:
        raw = []
        for tok in args.folder:
            raw.extend([t for t in re.split(r"\s*,\s*", tok) if t])
        folder = Path(raw[0]) if raw else Path(".")
        if not folder.is_dir():
            print(f"❌ Folder not found: {folder}")
            sys.exit(1)
    else:
        subdirs = sorted([p for p in Path(".").iterdir() if p.is_dir()],
                         key=lambda p: p.name.lower())
        if not subdirs:
            print("❌ No subdirectories here.")
            sys.exit(1)
        pick = choose("📁 Select the folder containing your CSV:", [p.name for p in subdirs])
        folder = next(p for p in subdirs if p.name == pick)

    csv_files = sorted(folder.glob("*.csv"), key=lambda p: p.name.lower())
    if not csv_files:
        print(f"❌ No CSV files found in {folder}.")
        sys.exit(1)
    if len(csv_files) == 1:
        return csv_files[0]

    pick = choose(f"📄 Select a CSV from {folder}:", [p.name for p in csv_files])
    return next(p for p in csv_files if p.name == pick)


def run_csv_mode(csv_path: Path, args, num_confs: int, num_workers: int, OBABEL: str, state_opts: dict):
    base_name = csv_path.stem

    tasks: list[dict[str, str]] = []
    seen_csv_ids: dict[str, int] = {}
    with open(csv_path, "r", newline="") as f:
        r = csv.DictReader(f)
        headers = r.fieldnames or []
        if not headers:
            print("❌ CSV has no header.")
            sys.exit(1)

        if args.smiles_col and args.id_col:
            smi_col = args.smiles_col
            id_col = args.id_col
        else:
            print("\n📊 CSV Columns:")
            for i, h in enumerate(headers):
                print(f" [{i}] {h}")
            smi_col = choose_csv_column(headers, "Index of SMILES column?", 0)
            id_col = choose_csv_column(headers, "Index of ID column?", 1 if len(headers) > 1 else 0)

        for row_index, row in enumerate(r, start=1):
            smiles = (row.get(smi_col) or "").strip()
            ligand_id = ensure_unique_ligand_id(
                sanitize_id((row.get(id_col) or "").strip()),
                seen_csv_ids,
                context=f"{csv_path.name} row {row_index}",
            )
            if not smiles or not ligand_id:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                tasks.append(
                    {
                        "ligand_id": ligand_id,
                        "molblock": Chem.MolToMolBlock(mol),
                        "source_input": str(csv_path),
                        "source_input_type": "csv",
                        "source_ligand_id": ligand_id,
                        "source_record_index": row_index,
                        "source_record_name": "",
                        "original_mol_name": "",
                    }
                )

    if not tasks:
        print("⚠️ No valid molecules found.")
        sys.exit(2)

    tag = f"{num_confs}Poses"
    pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

    print(f"\n🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
    print(f"📦 Output PDBQT: {pdbqt_dir}")
    print(f"🗂  TMP SDF:     {tmp_dir} (will {'NOT ' if args.remove_tmp else ''}be removed at end)")

    work_items = [
        (
            task["ligand_id"],
            task["molblock"],
            tmp_dir,
            pdbqt_dir,
            num_confs,
            OBABEL,
            state_opts,
            task["source_input"],
            task["source_input_type"],
            task["source_ligand_id"],
            task.get("source_record_index", ""),
            task.get("source_record_name", ""),
            task.get("original_mol_name", ""),
        )
        for task in tasks
    ]

    written_total = 0
    manifest_rows: list[dict] = []
    try:
        with Pool(processes=num_workers) as pool:
            for result in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                written_total += int((result or {}).get("written", 0))
                manifest_rows.extend((result or {}).get("rows", []))
    finally:
        if args.remove_tmp:
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
                print(f"🧹 TMP removed: {tmp_dir}")
            except Exception as e:
                print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

    print(f"\n✅ Done! {written_total} poses written to: {pdbqt_dir}")
    write_manifest_for_outputs(pdbqt_dir, manifest_rows)

# ----------------- per-folder runner -----------------
def run_one_folder(folder: Path, ft: str, num_confs: int, num_workers: int, OBABEL: str, remove_tmp: bool):
    print(f"\n=== Processing folder: {folder} (type={ft}) ===")

    tasks: list[dict[str, str]] = []
    base_name = folder.name

    # Collect ligands from the folder
    if ft == "sdf":
        for task in iter_mols_from_sdf_folder(folder):
            tasks.append(task)
    else:
        for task in iter_mols_from_smiles_folder(folder):
            tasks.append(task)

    if not tasks:
        print(f"⚠️  No valid molecules found in {folder}. Skipping.")
        return 0, None, None

    # Create output dirs
    tag = f"{num_confs}Poses"
    pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

    print(f"🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
    print(f"📦 Output PDBQT: {pdbqt_dir}")
    print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if remove_tmp else 'be KEPT'})")

    # # Build work items
    # work_items = [
    #     (ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL)
    #     for (ligand_id, molblock) in tasks
    # ]
    work_items = [
        (
            task["ligand_id"],
            task["molblock"],
            tmp_dir,
            pdbqt_dir,
            num_confs,
            OBABEL,
            state_opts,
            task["source_input"],
            task["source_input_type"],
            task["source_ligand_id"],
            task.get("source_record_index", ""),
            task.get("source_record_name", ""),
            task.get("original_mol_name", ""),
        )
        for task in tasks
    ]

    written_total = 0
    manifest_rows: list[dict] = []
    try:
        with Pool(processes=num_workers) as pool:
            for result in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                written_total += int((result or {}).get("written", 0))
                manifest_rows.extend((result or {}).get("rows", []))
    finally:
        # 🧹 Only remove TMP if explicitly requested
        if remove_tmp:
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
                print(f"🧹 TMP removed: {tmp_dir}")
            except Exception as e:
                print(f"⚠️ Could not remove TMP dir {tmp_dir}: {e}", file=sys.stderr)

    print(f"✅ Done! {written_total} poses written to: {pdbqt_dir}")
    write_manifest_for_outputs(pdbqt_dir, manifest_rows)
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
    # =============================
    # 5) ADD THIS IN main
    # after OBABEL = obabel_path(...)
    # =============================

    if args.enumerate_protomers and not HAVE_DIMORPHITE:
        print("❌ --enumerate-protomers requires dimorphite_dl.")
        print("   Install with: pip install dimorphite_dl")
        sys.exit(2)

    state_opts = {
        "enumerate_protomers": args.enumerate_protomers,
        "ph_min": args.ph_min,
        "ph_max": args.ph_max,
        "ph_precision": args.ph_precision,
        "max_protomers": args.max_protomers,
        "max_tautomers": args.max_tautomers,
        "max_transforms": args.max_transforms,
    }

    # Interactive fallback if mode not provided
    if mode is None:
        print("📘 Choose input type:")
        print(" [1] CSV file (SMILES + ID columns)")
        print(" [2] Folder of ligands (.sdf or .smiles)  ← supports MULTIPLE folders")
        print(" [3] Single .sdf file (multi-record)")
        print(" [4] Folder containing a CSV file")
        mode = input("Enter 1, 2, 3, or 4: ").strip()

    # --- Mode 1: CSV ---
    if mode == "1":
        run_csv_mode(resolve_csv_path(args), args, num_confs, num_workers, OBABEL, state_opts)

    # --- Mode 4: Folder containing a CSV ---
    elif mode == "4":
        run_csv_mode(resolve_csv_path_from_folder(args), args, num_confs, num_workers, OBABEL, state_opts)

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

        tasks: list[dict[str, str]] = []
        for task in iter_mols_from_sdf_file(sdf_path):
            tasks.append(task)

        if not tasks:
            print("⚠️ No valid molecules found."); sys.exit(2)

        tag = f"{num_confs}Poses"
        pdbqt_dir, tmp_dir = create_output_dirs(base_name, tag)

        print(f"\n🚀 Generating {num_confs} poses for {len(tasks)} molecules using {num_workers} cores…")
        print(f"📦 Output PDBQT: {pdbqt_dir}")
        print(f"🗂  TMP SDF:     {tmp_dir} (will {'be removed' if args.remove_tmp else 'be KEPT'})")


        # work_items = [(ligand_id, molblock, tmp_dir, pdbqt_dir, num_confs, OBABEL) for (ligand_id, molblock) in tasks]
        work_items = [
            (
                task["ligand_id"],
                task["molblock"],
                tmp_dir,
                pdbqt_dir,
                num_confs,
                OBABEL,
                state_opts,
                task["source_input"],
                task["source_input_type"],
                task["source_ligand_id"],
                task.get("source_record_index", ""),
                task.get("source_record_name", ""),
                task.get("original_mol_name", ""),
            )
            for task in tasks
        ]

        written_total = 0
        manifest_rows: list[dict] = []
        try:
            with Pool(processes=num_workers) as pool:
                for result in pool.imap_unordered(generate_poses, work_items, chunksize=1):
                    written_total += int((result or {}).get("written", 0))
                    manifest_rows.extend((result or {}).get("rows", []))
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
        write_manifest_for_outputs(pdbqt_dir, manifest_rows)

    else:
        print("❌ Invalid selection."); sys.exit(1)

