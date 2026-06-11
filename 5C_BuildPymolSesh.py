#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ligand_manifest import (
    CHEMICAL_METADATA_COLUMNS,
    find_ligand_state_manifests,
    load_ligand_state_manifest,
    merge_ligand_metadata,
)
from ligand_naming import parse_ligand_variant

try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS, rdMolAlign
    from rdkit.Chem.rdmolfiles import SDWriter
    HAVE_RDKIT = True
except Exception:  # pragma: no cover - optional runtime dependency
    Chem = Any  # type: ignore[assignment]
    rdFMCS = None  # type: ignore[assignment]
    rdMolAlign = None  # type: ignore[assignment]
    SDWriter = None  # type: ignore[assignment]
    HAVE_RDKIT = False

try:  # pragma: no cover - depends on runtime
    import pymol as _pymol
    HAVE_PYMOL = True
except Exception:  # pragma: no cover - depends on runtime
    _pymol = None
    HAVE_PYMOL = False

cmd = None


AUDIT_COLUMNS = [
    "timestamp",
    "csv_style",
    "selection_mode",
    "receptor",
    "rank_for_receptor",
    "base_ligand",
    "ligand_variant",
    "state_tag",
    "protomer_tag",
    "tautomer_tag",
    "conformer_tag",
    "vina_pose",
    "binding_affinity",
    "dock_pdbqt",
    "selected_pose_pdbqt",
    "reference_sdf",
    "reference_lookup_method",
    "reference_lookup_warning",
    "coordinate_transfer_method",
    "coordinate_transfer_warning",
    "hydrogen_mode",
    "hydrogen_policy_result",
    "heavy_atom_count_before_h_policy",
    "heavy_atom_count_after_h_policy",
    "total_atom_count_before_h_policy",
    "total_atom_count_after_h_policy",
    "hydrogen_count_before_h_policy",
    "hydrogen_count_after_h_policy",
    "receptor_file",
    "saved_complex_pdb",
    "saved_corrected_sdf",
    "saved_reference_sdf",
    "saved_raw_pdbqt",
    "pymol_object",
    "pymol_group",
    "state_index",
] + CHEMICAL_METADATA_COLUMNS


def build_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build PyMOL sessions and top-hit exports using exact generated SDF chemistry plus docked PDBQT coordinates."
    )
    parser.add_argument("--csv", help="Docking score CSV to use")
    parser.add_argument("--mode", choices=["per_ligand", "per_receptor"], help="Selection mode")
    parser.add_argument("--top", type=int, help="Top N rows to keep for the selected mode")
    parser.add_argument("--receptor-roots", nargs="+", help="One or more receptor root folders")
    parser.add_argument("--outdir", help="Output directory root (defaults to current directory)")
    parser.add_argument("--obabel-bin", help="Path to obabel executable")
    parser.add_argument("--hydrogen-mode", choices=["none", "nonpolar", "all"], default="nonpolar", help="Hydrogen cleanup policy for corrected docked ligand exports")
    parser.add_argument("--no-hydrogens", action="store_true", help="Shortcut for --hydrogen-mode none")
    parser.add_argument("--non-interactive", action="store_true", help="Fail instead of prompting if required args are missing")
    return parser.parse_args()


def safe_float(x: Any, default: float = 1e9) -> float:
    try:
        return float(x)
    except Exception:
        return default


def sanitize_filename(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")
    return cleaned or "item"


def sanitize_pymol_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    if not cleaned or not cleaned[0].isalpha():
        cleaned = f"obj_{cleaned or 'item'}"
    return cleaned[:120]


def obabel_path(cli_path: Optional[str]) -> str:
    found = cli_path or os.environ.get("OBABEL_BIN") or shutil.which("obabel")
    if not found:
        raise SystemExit(
            "❌ Cannot find Open Babel executable. Install `obabel` or set OBABEL_BIN/--obabel-bin."
        )
    return found


def launch_pymol() -> bool:  # pragma: no cover - depends on runtime
    global cmd
    if not HAVE_PYMOL:
        print("⚠️ PyMOL is not available; session and complex PDB exports will be skipped.")
        return False
    try:
        _pymol.pymol_argv = ["pymol", "-cq"]
        _pymol.finish_launching()
    except Exception:
        pass
    from pymol import cmd as pymol_cmd  # type: ignore

    cmd = pymol_cmd
    try:
        cmd.reinitialize()
        cmd.feedback("disable", "all", "warnings")
        cmd.set("pdb_conect_all", 0)
        cmd.set("pdb_conect_nodup", 1)
    except Exception:
        pass
    return True


def list_and_select(prompt: str, options: Sequence[str]) -> str:
    print(f"\n📁 {prompt}")
    for idx, opt in enumerate(options, 1):
        print(f"{idx}. {opt}")
    chosen = int(input("🔢 Select: ").strip()) - 1
    return options[chosen]


def list_and_select_multiple(prompt: str, options: Sequence[str]) -> List[str]:
    print(f"\n📁 {prompt}")
    for idx, opt in enumerate(options, 1):
        print(f"{idx}. {opt}")
    print("\n🔢 Select one or more by number, comma-separated, or type 'all'.")
    raw = input("🔢 Select: ").strip()
    if raw.lower() == "all":
        return list(options)
    picks: List[str] = []
    seen = set()
    for piece in raw.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if not piece.isdigit():
            raise SystemExit(f"❌ Invalid selection: {piece}")
        idx = int(piece) - 1
        if idx < 0 or idx >= len(options):
            raise SystemExit(f"❌ Selection out of range: {piece}")
        value = options[idx]
        if value not in seen:
            picks.append(value)
            seen.add(value)
    if not picks:
        raise SystemExit("❌ No receptors selected.")
    return picks


def detect_csv_style(rows: Sequence[Dict[str, str]]) -> str:
    if not rows:
        return "simple"
    if any((row.get("OutFile") or "").strip() for row in rows):
        return "provenance"
    if any((row.get("ResultsRoot") or "").strip() for row in rows):
        return "provenance"
    return "simple"


def normalize_score_row(row: Dict[str, str]) -> Dict[str, Any]:
    ligand_variant = (
        (row.get("LigandVariant") or "").strip()
        or Path((row.get("OutFile") or "").strip()).parent.name
        or (row.get("LigandDir") or "").strip()
        or (row.get("Ligand") or "").strip()
    )
    merged = merge_ligand_metadata(ligand_variant, base_row=dict(row))
    merged["LigandVariant"] = ligand_variant or merged.get("LigandVariant") or merged.get("Ligand") or ""
    merged["LigandBase"] = merged.get("LigandBase") or merged.get("Ligand") or parse_ligand_variant(merged["LigandVariant"])["LigandBase"]
    merged["Ligand"] = merged.get("Ligand") or merged["LigandBase"]
    return merged


def load_csv_rows(csv_path: Path) -> Tuple[List[Dict[str, Any]], str]:
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = [normalize_score_row(row) for row in csv.DictReader(handle)]
    if not rows:
        raise SystemExit("❌ CSV is empty.")
    for key in ("Receptor", "Ligand", "Pose"):
        if key not in rows[0]:
            raise SystemExit(f"❌ CSV missing required column: {key}")
    bind_col = "Binding_Affinity" if "Binding_Affinity" in rows[0] else "Binding_Affinity_kcal_per_mol"
    rows.sort(key=lambda row: safe_float(row.get(bind_col, "1e9")))
    return rows, detect_csv_style(rows)


def select_rows(rows: Sequence[Dict[str, Any]], top_n: int, mode: str) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    if mode == "per_ligand":
        counts: Dict[Tuple[str, str], int] = defaultdict(int)
        for row in rows:
            key = (row["Receptor"], row["LigandBase"])
            if counts[key] >= top_n:
                continue
            grouped[row["Receptor"]].append(row)
            counts[key] += 1
        return grouped

    if mode == "per_receptor":
        seen_bases: Dict[str, set[str]] = defaultdict(set)
        for row in rows:
            receptor = row["Receptor"]
            base = row["LigandBase"]
            if base in seen_bases[receptor] or len(seen_bases[receptor]) >= top_n:
                continue
            grouped[receptor].append(row)
            seen_bases[receptor].add(base)
        return grouped

    raise SystemExit(f"❌ Unknown mode: {mode}")


def resolve_outfile_unified(row: Dict[str, Any], cwd: Path) -> Optional[Path]:
    out = (row.get("OutFile") or "").strip()
    if out:
        p = Path(out)
        if p.is_file():
            return p
        if p.is_absolute():
            tail = Path(*p.parts[-4:]) if len(p.parts) >= 4 else p
            candidate = cwd / tail
            if candidate.is_file():
                return candidate

    rr = (row.get("ResultsRoot") or "").strip()
    rd = (row.get("ReceptorDir") or "").strip()
    ld = (row.get("LigandDir") or "").strip()
    if rr and rd and ld:
        guess = cwd / rr / rd / ld / "out.pdbqt"
        if guess.is_file():
            return guess

    ligand_variant = row.get("LigandVariant") or row.get("Ligand") or ""
    receptor = row.get("Receptor") or ""
    for root in cwd.glob("Docking_Results*"):
        candidate = root / receptor / ligand_variant / "out.pdbqt"
        if candidate.is_file():
            return candidate
    return None


def infer_tmp_sdf_relative_path(ligand_variant: str) -> Optional[Path]:
    parsed = parse_ligand_variant(ligand_variant)
    base = parsed["LigandBase"]
    state = parsed["StateTag"]
    if not base or not state:
        return None
    return Path(base) / f"{base}_{state}" / f"{ligand_variant}.sdf"


def _suffix_after(marker: str, text: str) -> str:
    if marker not in text:
        return ""
    return text.split(marker, 1)[1]


class ReferenceResolver:
    def __init__(self, cwd: Path):
        self.cwd = Path(cwd).resolve()
        self.tmp_roots = sorted([p.resolve() for p in self.cwd.glob("Ligands_TMP_SDF_*") if p.is_dir()])
        self.variant_sdf_index: Dict[str, List[Path]] = defaultdict(list)
        self.base_sdf_index: Dict[str, List[Path]] = defaultdict(list)
        self.manifest_index: Dict[str, List[Tuple[Path, Dict[str, str]]]] = defaultdict(list)
        self._index_files()

    def _index_files(self):
        search_roots: List[Path] = self.tmp_roots[:]
        search_roots.extend([p.resolve() for p in self.cwd.glob("Ligands*") if p.is_dir()])
        out_ligands = self.cwd / "11_AA_Ligands"
        if out_ligands.exists():
            search_roots.append(out_ligands.resolve())
        seen = set()
        for root in search_roots:
            root_key = str(root)
            if root_key in seen:
                continue
            seen.add(root_key)
            for path in root.rglob("*.sdf"):
                if any(part.startswith(".") or part == "__MACOSX" for part in path.parts):
                    continue
                resolved = path.resolve()
                self.variant_sdf_index[resolved.stem].append(resolved)
                parsed = parse_ligand_variant(resolved.stem)
                self.base_sdf_index[parsed["LigandBase"]].append(resolved)

        for manifest_path in find_ligand_state_manifests([self.cwd]):
            for variant, row in load_ligand_state_manifest(manifest_path).items():
                self.manifest_index[variant].append((manifest_path.resolve(), row))

    def _manifest_candidate_paths(self, manifest_path: Path, row: Dict[str, str]) -> Iterable[Tuple[Path, str]]:
        manifest_dir = manifest_path.parent
        suffix = _suffix_after("_PDBQT_", manifest_dir.name)
        aligned_tmp_roots = [root for root in self.tmp_roots if suffix and root.name.endswith(suffix)]
        if not aligned_tmp_roots:
            aligned_tmp_roots = self.tmp_roots

        for field in ("TmpSDFFile", "SDFFile", "SourceSDF"):
            value = (row.get(field) or "").strip()
            if not value:
                continue
            raw = Path(value)
            if raw.is_absolute():
                yield raw, f"manifest:{field}"
                continue
            for root in aligned_tmp_roots:
                yield root / raw, f"manifest:{field}"
            yield manifest_dir / raw, f"manifest:{field}"
            yield self.cwd / raw, f"manifest:{field}"
            if field == "SDFFile":
                for root in aligned_tmp_roots:
                    yield root / raw.name, f"manifest:{field}"

    def resolve(self, row: Dict[str, Any]) -> Tuple[Optional[Path], str, str, Dict[str, str]]:
        ligand_variant = (row.get("LigandVariant") or "").strip()
        base_ligand = (row.get("LigandBase") or "").strip() or parse_ligand_variant(ligand_variant)["LigandBase"]

        for manifest_path, manifest_row in self.manifest_index.get(ligand_variant, []):
            for candidate, method in self._manifest_candidate_paths(manifest_path, manifest_row):
                if candidate.exists():
                    return candidate.resolve(), method, "", manifest_row

        exact_paths = self.variant_sdf_index.get(ligand_variant, [])
        if exact_paths:
            preferred = sorted(exact_paths, key=lambda p: ("Ligands_TMP_SDF_" not in str(p), len(str(p))))
            return preferred[0], "exact_filename_search", "", {}

        inferred = infer_tmp_sdf_relative_path(ligand_variant)
        if inferred is not None:
            for root in self.tmp_roots:
                candidate = root / inferred
                if candidate.exists():
                    return candidate.resolve(), "structured_tmp_inference", "", {}

        source_input = (row.get("SourceInput") or "").strip()
        if source_input:
            direct = self.cwd / source_input
            if direct.exists():
                warning = "Fell back to base ligand source SDF because exact variant SDF was not found."
                return direct.resolve(), "base_source_fallback", warning, {}
            named = self.cwd / "Ligands" / Path(source_input).name
            if named.exists():
                warning = "Fell back to base ligand source SDF because exact variant SDF was not found."
                return named.resolve(), "base_source_fallback", warning, {}

        base_paths = self.base_sdf_index.get(base_ligand, [])
        if base_paths:
            warning = "Fell back to base ligand SDF because exact variant SDF was not found."
            preferred = sorted(base_paths, key=lambda p: ("Ligands_TMP_SDF_" in str(p), len(str(p))), reverse=True)
            return preferred[0], "base_source_fallback", warning, {}

        return None, "not_found", "No exact or fallback ligand reference SDF was found.", {}


def extract_selected_pose_pdbqt(source_pdbqt: Path, pose_index: int, dest_pdbqt: Path) -> Tuple[Path, str]:
    text = source_pdbqt.read_text()
    lines = text.splitlines()
    models: List[Tuple[int, List[str]]] = []
    current_number: Optional[int] = None
    current_lines: List[str] = []

    for line in lines:
        if line.startswith("MODEL"):
            if current_lines:
                models.append((current_number or len(models) + 1, current_lines[:]))
            current_lines = [line]
            match = re.search(r"MODEL\s+(\d+)", line)
            current_number = int(match.group(1)) if match else len(models) + 1
            continue
        if current_lines:
            current_lines.append(line)
            if line.startswith("ENDMDL"):
                models.append((current_number or len(models) + 1, current_lines[:]))
                current_lines = []
                current_number = None
    if current_lines:
        models.append((current_number or len(models) + 1, current_lines[:]))

    dest_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if not models:
        shutil.copy2(source_pdbqt, dest_pdbqt)
        return dest_pdbqt, "No MODEL records found; copied whole docked PDBQT."

    wanted = pose_index if pose_index > 0 else 1
    selected = next((block for number, block in models if number == wanted), None)
    if selected is None:
        selected = models[0][1]
        warning = f"Requested pose {wanted} not found; used MODEL {models[0][0]} instead."
    else:
        warning = ""
    dest_pdbqt.write_text("\n".join(selected) + "\n")
    return dest_pdbqt, warning


def load_reference_mol(reference_sdf: Path) -> Optional[Chem.Mol]:
    if not HAVE_RDKIT:
        return None
    try:
        supplier = Chem.SDMolSupplier(str(reference_sdf), sanitize=False, removeHs=False)
        return next((mol for mol in supplier if mol is not None), None)
    except Exception:
        return None


def pdbqt_pose_to_mol(pdbqt_file: Path, obabel_bin: str) -> Chem.Mol:
    if not HAVE_RDKIT:
        raise RuntimeError("RDKit is required for coordinate transfer.")
    with tempfile.TemporaryDirectory(prefix="dock_pose_") as tmpdir:
        mol_path = Path(tmpdir) / "dock_pose.mol"
        subprocess.run(
            [obabel_bin, str(pdbqt_file), "-O", str(mol_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        mol = Chem.MolFromMolFile(str(mol_path), sanitize=False, removeHs=False)
        if mol is None:
            raise RuntimeError(f"Open Babel conversion did not yield a molecule for {pdbqt_file}")
        return mol


def atom_count_summary(mol: Chem.Mol) -> Dict[str, int]:
    total = mol.GetNumAtoms()
    hydrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    return {
        "total": total,
        "hydrogen": hydrogen,
        "heavy": total - hydrogen,
    }


def _heavy_atom_signature(mol: Chem.Mol) -> Tuple[List[int], List[Tuple[float, float, float]]]:
    conf = mol.GetConformer()
    atomic_numbers: List[int] = []
    coords: List[Tuple[float, float, float]] = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        atomic_numbers.append(atom.GetAtomicNum())
        coords.append((float(pos.x), float(pos.y), float(pos.z)))
    return atomic_numbers, coords


def heavy_atoms_preserved(before: Chem.Mol, after: Chem.Mol, tol: float = 1e-4) -> bool:
    before_nums, before_coords = _heavy_atom_signature(before)
    after_nums, after_coords = _heavy_atom_signature(after)
    if before_nums != after_nums or len(before_coords) != len(after_coords):
        return False
    for (bx, by, bz), (ax, ay, az) in zip(before_coords, after_coords):
        if abs(bx - ax) > tol or abs(by - ay) > tol or abs(bz - az) > tol:
            return False
    return True


def strip_all_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    return Chem.RemoveHs(Chem.Mol(mol), sanitize=False)


def add_nonpolar_hydrogens_only(mol: Chem.Mol) -> Chem.Mol:
    no_h = strip_all_hydrogens(mol)
    with_h = Chem.AddHs(Chem.Mol(no_h), addCoords=True)
    editable = Chem.RWMol(with_h)
    remove_idxs: List[int] = []
    for atom in editable.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        neighbors = atom.GetNeighbors()
        neighbor = neighbors[0] if neighbors else None
        if neighbor is None or neighbor.GetAtomicNum() != 6:
            remove_idxs.append(atom.GetIdx())
    for idx in sorted(remove_idxs, reverse=True):
        editable.RemoveAtom(idx)
    return editable.GetMol()


def apply_hydrogen_policy(mol: Chem.Mol, mode: str) -> Tuple[Chem.Mol, str]:
    requested = mode or "nonpolar"
    if requested == "none":
        processed = strip_all_hydrogens(mol)
        if not heavy_atoms_preserved(mol, processed):
            fallback = strip_all_hydrogens(mol)
            return fallback, "Requested no hydrogens; heavy atoms preserved check failed unexpectedly, kept hydrogen-stripped fallback."
        return processed, "Removed all hydrogens."

    if requested == "all":
        processed = Chem.Mol(mol)
        if not heavy_atoms_preserved(mol, processed):
            fallback = strip_all_hydrogens(mol)
            return fallback, "Preserve-all-hydrogen mode failed heavy-atom preservation check; fell back to no hydrogens."
        return processed, "Preserved all hydrogens."

    processed = add_nonpolar_hydrogens_only(mol)
    if not heavy_atoms_preserved(mol, processed):
        fallback = strip_all_hydrogens(mol)
        return fallback, "Nonpolar-hydrogen mode failed heavy-atom preservation check; fell back to no hydrogens."
    return processed, "Removed all hydrogens and added back only carbon-bound hydrogens."


def _copy_coords_by_atom_map(reference_mol: Chem.Mol, dock_mol: Chem.Mol, atom_map: Sequence[Tuple[int, int]]) -> Chem.Mol:
    copied = Chem.Mol(reference_mol)
    ref_conf = copied.GetConformer()
    dock_conf = dock_mol.GetConformer()
    for ref_idx, dock_idx in atom_map:
        ref_conf.SetAtomPosition(ref_idx, dock_conf.GetAtomPosition(dock_idx))
    return copied


def _full_atom_order_map(reference_mol: Chem.Mol, dock_mol: Chem.Mol) -> Optional[List[Tuple[int, int]]]:
    if reference_mol.GetNumAtoms() != dock_mol.GetNumAtoms():
        return None
    ref_atoms = [atom.GetAtomicNum() for atom in reference_mol.GetAtoms()]
    dock_atoms = [atom.GetAtomicNum() for atom in dock_mol.GetAtoms()]
    if ref_atoms != dock_atoms:
        return None
    return [(idx, idx) for idx in range(reference_mol.GetNumAtoms())]


def _no_h_with_map(mol: Chem.Mol) -> Tuple[Chem.Mol, List[int]]:
    full_idx = [idx for idx, atom in enumerate(mol.GetAtoms()) if atom.GetAtomicNum() > 1]
    return Chem.RemoveHs(Chem.Mol(mol), sanitize=False), full_idx


def _heavy_atom_mcs_map(reference_mol: Chem.Mol, dock_mol: Chem.Mol) -> List[Tuple[int, int]]:
    if not HAVE_RDKIT or rdFMCS is None:
        return []
    ref_no_h, ref_map = _no_h_with_map(reference_mol)
    dock_no_h, dock_map = _no_h_with_map(dock_mol)
    mcs = rdFMCS.FindMCS([ref_no_h, dock_no_h], timeout=10)
    if mcs.canceled or not mcs.smartsString:
        return []
    patt = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_no_h.GetSubstructMatch(patt)
    dock_match = dock_no_h.GetSubstructMatch(patt)
    if not ref_match or len(ref_match) != len(dock_match):
        return []
    return [(ref_map[ref_idx], dock_map[dock_idx]) for ref_idx, dock_idx in zip(ref_match, dock_match)]


def _rigid_align_with_map(reference_mol: Chem.Mol, dock_mol: Chem.Mol, atom_map: Sequence[Tuple[int, int]]) -> Optional[Chem.Mol]:
    if not atom_map or not HAVE_RDKIT or rdMolAlign is None:
        return None
    aligned = Chem.Mol(reference_mol)
    try:
        rdMolAlign.AlignMol(aligned, dock_mol, atomMap=list(atom_map))
        return aligned
    except Exception:
        return None


def transfer_docked_coordinates(reference_mol: Optional[Chem.Mol], dock_mol: Chem.Mol) -> Tuple[Chem.Mol, str, str]:
    if reference_mol is None:
        return dock_mol, "coords_only", "No reference chemistry available; using docked coordinates only."

    atom_order_map = _full_atom_order_map(reference_mol, dock_mol)
    if atom_order_map:
        return _copy_coords_by_atom_map(reference_mol, dock_mol, atom_order_map), "atom_order", ""

    mcs_map = _heavy_atom_mcs_map(reference_mol, dock_mol)
    if mcs_map:
        aligned = _rigid_align_with_map(reference_mol, dock_mol, mcs_map)
        if aligned is not None:
            corrected = _copy_coords_by_atom_map(aligned, dock_mol, mcs_map)
            return corrected, "heavy_atom_mcs", "Used heavy-atom MCS coordinate transfer because atom order did not match."
        rigid_only = _rigid_align_with_map(reference_mol, dock_mol, mcs_map)
        if rigid_only is not None:
            return rigid_only, "rigid_mcs_fit", "Used rigid MCS fit because direct coordinate transfer was not possible."

    return dock_mol, "coords_only", "Failed to map docked coordinates onto reference chemistry; using coordinates-only fallback."


def write_sdf(path: Path, mol: Chem.Mol):
    if not HAVE_RDKIT or SDWriter is None:
        raise RuntimeError("RDKit is required to write SDF output.")
    writer = SDWriter(str(path))
    writer.SetKekulize(False)
    writer.write(mol)
    writer.close()


def find_receptor_file_in_root(receptors_root: Path, rec_name: str) -> Optional[Path]:
    stems = {rec_name, rec_name.replace(".pdbqt", "").replace(".converted", "")}
    for root, _, files in os.walk(receptors_root):
        for filename in files:
            candidate = Path(root) / filename
            if candidate.stem in stems or candidate.name.startswith(rec_name):
                return candidate
    return None


def find_receptor_file(receptor_roots: Sequence[Path], rec_name: str) -> Optional[Path]:
    for root in receptor_roots:
        hit = find_receptor_file_in_root(root, rec_name)
        if hit is not None:
            return hit
    return None


def build_pymol_names(receptor: str, ligand_base: str, ligand_variant: str, pose_index: int) -> Dict[str, str]:
    receptor_obj = sanitize_pymol_name(f"obj_{receptor}")
    receptor_group = sanitize_pymol_name(f"grp_{receptor}")
    ligand_group = sanitize_pymol_name(f"grp_{receptor}_{ligand_base}")
    ligand_obj = sanitize_pymol_name(f"{receptor}_{ligand_variant}_pose{pose_index}")
    return {
        "receptor_obj": receptor_obj,
        "receptor_group": receptor_group,
        "ligand_group": ligand_group,
        "ligand_obj": ligand_obj,
    }


def ensure_pymol_receptor_loaded(pymol_enabled: bool, receptor_file: Path, receptor: str):
    if not pymol_enabled:
        return
    names = build_pymol_names(receptor, receptor, receptor, 1)
    try:
        if names["receptor_obj"] in set(cmd.get_names("objects")):
            return
    except Exception:
        pass
    cmd.load(str(receptor_file), names["receptor_obj"])
    try:
        cmd.hide("everything", f"{names['receptor_obj']} and hydro")
    except Exception:
        pass
    cmd.group(names["receptor_group"], names["receptor_obj"])


def write_audit_csv(path: Path, rows: Sequence[Dict[str, Any]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_COLUMNS)
        writer.writeheader()
        for row in rows:
            cooked = {key: row.get(key, "") for key in AUDIT_COLUMNS}
            writer.writerow(cooked)


def prompt_for_inputs(args: argparse.Namespace, cwd: Path) -> Tuple[Path, str, int, List[Path], Path]:
    csv_path = Path(args.csv).expanduser() if args.csv else None
    mode = args.mode
    top_n = args.top
    receptor_roots = [Path(p).expanduser() for p in (args.receptor_roots or [])]
    out_root = Path(args.outdir).expanduser() if args.outdir else cwd

    if args.non_interactive:
        if not csv_path or not mode or not top_n or not receptor_roots:
            raise SystemExit("❌ Non-interactive mode requires --csv, --mode, --top, and --receptor-roots.")
        return csv_path, mode, int(top_n), receptor_roots, out_root

    if not csv_path:
        csv_files = sorted([f.name for f in cwd.glob("*.csv")])
        if not csv_files:
            raise SystemExit("❌ No CSV files found in current directory.")
        csv_path = cwd / list_and_select("Available CSV files:", csv_files)

    if not mode:
        print("\n⚙️ Selection mode:")
        print("1. Top N poses per ligand")
        print("2. Top N ligands total")
        mode_choice = input("🔢 Choose [1/2]: ").strip()
        mode = "per_ligand" if mode_choice == "1" else ("per_receptor" if mode_choice == "2" else None)
        if not mode:
            raise SystemExit("❌ Invalid choice.")

    if not top_n:
        top_n = int(input("🔝 How many top entries? (e.g., 5): ").strip())

    if not receptor_roots:
        rec_dirs = sorted([d.name for d in cwd.iterdir() if d.is_dir() and "Receptor" in d.name])
        if not rec_dirs:
            raise SystemExit("❌ No receptor folders found in current directory.")
        selected = list_and_select_multiple("Receptor folders (choose one or more):", rec_dirs)
        receptor_roots = [cwd / item for item in selected]

    return csv_path, mode, int(top_n), receptor_roots, out_root


def main():
    args = build_args()
    cwd = Path(".").resolve()
    hydrogen_mode = "none" if args.no_hydrogens else args.hydrogen_mode
    csv_path, mode, top_n, receptor_roots, out_root = prompt_for_inputs(args, cwd)
    out_root = out_root.resolve()

    if not HAVE_RDKIT:
        raise SystemExit("❌ RDKit is required for corrected ligand reconstruction.")

    obabel_bin = obabel_path(args.obabel_bin)
    out_dir_pdb = out_root / "11_AA_PDBOUTPUT"
    out_dir_sdf = out_root / "11_AA_Ligands"
    out_dir_pdbqt = out_root / "11_AA_PDBQT"
    out_dir_pdb.mkdir(parents=True, exist_ok=True)
    out_dir_sdf.mkdir(parents=True, exist_ok=True)
    out_dir_pdbqt.mkdir(parents=True, exist_ok=True)

    rows, csv_style = load_csv_rows(csv_path)
    grouped_rows = select_rows(rows, top_n, mode)
    resolver = ReferenceResolver(cwd)
    pymol_enabled = launch_pymol()

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    audit_rows: List[Dict[str, Any]] = []

    for receptor, selected_rows in grouped_rows.items():
        receptor_file = find_receptor_file(receptor_roots, receptor)
        if receptor_file is None:
            print(f"⚠️ Receptor file not found for {receptor}; skipping receptor.")
            continue

        ensure_pymol_receptor_loaded(pymol_enabled, receptor_file, receptor)
        rank_for_receptor = 0

        for row in selected_rows:
            rank_for_receptor += 1
            ligand_variant = row["LigandVariant"]
            ligand_base = row["LigandBase"]
            state_tag = row.get("StateTag", "")
            pose_index = int(float(row.get("Pose") or 1))

            dock_pdbqt = resolve_outfile_unified(row, cwd)
            if dock_pdbqt is None or not dock_pdbqt.exists():
                print(f"⚠️ Missing docked PDBQT for {receptor}/{ligand_variant}; skipping.")
                continue

            prefix = sanitize_filename(
                f"{receptor}_Top{rank_for_receptor:02d}_{ligand_variant}_vinaPose{pose_index}"
            )
            selected_pose_pdbqt = out_dir_pdbqt / f"{prefix}_raw_docked.pdbqt"
            _, pose_warning = extract_selected_pose_pdbqt(dock_pdbqt, pose_index, selected_pose_pdbqt)

            reference_sdf, ref_method, ref_warning, manifest_row = resolver.resolve(row)
            merged_row = merge_ligand_metadata(ligand_variant, base_row=row, manifest_row=manifest_row)

            reference_mol = load_reference_mol(reference_sdf) if reference_sdf else None
            if reference_sdf and reference_mol is None:
                ref_warning = (ref_warning + " " if ref_warning else "") + "Reference SDF could not be parsed."

            try:
                dock_mol = pdbqt_pose_to_mol(selected_pose_pdbqt, obabel_bin)
            except Exception as exc:
                print(f"⚠️ Open Babel/RDKit failed for {selected_pose_pdbqt.name}: {exc}")
                continue

            corrected_mol, coord_method, coord_warning = transfer_docked_coordinates(reference_mol, dock_mol)
            if pose_warning:
                coord_warning = f"{coord_warning} {pose_warning}".strip()
            counts_before = atom_count_summary(corrected_mol)
            corrected_mol, hydrogen_policy_result = apply_hydrogen_policy(corrected_mol, hydrogen_mode)
            counts_after = atom_count_summary(corrected_mol)

            corrected_sdf = out_dir_sdf / f"{prefix}_corrected_docked.sdf"
            reference_copy: Optional[Path] = out_dir_sdf / f"{prefix}_reference.sdf"
            complex_pdb = out_dir_pdb / f"{prefix}_complex.pdb"

            write_sdf(corrected_sdf, corrected_mol)
            if reference_sdf and reference_sdf.exists():
                shutil.copy2(reference_sdf, reference_copy)
            else:
                reference_copy = None

            pymol_names = build_pymol_names(receptor, ligand_base, ligand_variant, pose_index)
            pymol_group = pymol_names["ligand_group"]
            ligand_object = pymol_names["ligand_obj"]
            saved_complex = ""

            if pymol_enabled:
                try:
                    cmd.load(str(corrected_sdf), ligand_object)
                    cmd.hide("everything", ligand_object)
                    cmd.show("sticks", ligand_object)
                    cmd.group(pymol_group, ligand_object)
                    cmd.group(pymol_names["receptor_group"], pymol_group)
                    selection = f"({pymol_names['receptor_obj']}) or ({ligand_object})"
                    cmd.save(str(complex_pdb), selection=selection)
                    saved_complex = str(complex_pdb)
                except Exception as exc:
                    print(f"⚠️ PyMOL export failed for {ligand_variant}: {exc}")

            audit_row: Dict[str, Any] = {
                "timestamp": timestamp,
                "csv_style": csv_style,
                "selection_mode": mode,
                "receptor": receptor,
                "rank_for_receptor": rank_for_receptor,
                "base_ligand": ligand_base,
                "ligand_variant": ligand_variant,
                "state_tag": state_tag,
                "protomer_tag": merged_row.get("ProtomerTag", ""),
                "tautomer_tag": merged_row.get("TautomerTag", ""),
                "conformer_tag": merged_row.get("ConformerTag", ""),
                "vina_pose": pose_index,
                "binding_affinity": merged_row.get("Binding_Affinity", merged_row.get("Binding_Affinity_kcal_per_mol", "")),
                "dock_pdbqt": str(dock_pdbqt),
                "selected_pose_pdbqt": str(selected_pose_pdbqt),
                "reference_sdf": str(reference_sdf) if reference_sdf else "",
                "reference_lookup_method": ref_method,
                "reference_lookup_warning": ref_warning,
                "coordinate_transfer_method": coord_method,
                "coordinate_transfer_warning": coord_warning,
                "hydrogen_mode": hydrogen_mode,
                "hydrogen_policy_result": hydrogen_policy_result,
                "heavy_atom_count_before_h_policy": counts_before["heavy"],
                "heavy_atom_count_after_h_policy": counts_after["heavy"],
                "total_atom_count_before_h_policy": counts_before["total"],
                "total_atom_count_after_h_policy": counts_after["total"],
                "hydrogen_count_before_h_policy": counts_before["hydrogen"],
                "hydrogen_count_after_h_policy": counts_after["hydrogen"],
                "receptor_file": str(receptor_file),
                "saved_complex_pdb": saved_complex,
                "saved_corrected_sdf": str(corrected_sdf),
                "saved_reference_sdf": str(reference_copy) if reference_copy is not None else "",
                "saved_raw_pdbqt": str(selected_pose_pdbqt),
                "pymol_object": ligand_object,
                "pymol_group": pymol_group,
                "state_index": 1,
            }
            for key in CHEMICAL_METADATA_COLUMNS:
                audit_row[key] = merged_row.get(key, "")
            audit_rows.append(audit_row)

    session_path = out_dir_pdb / f"Top{top_n}_{mode}_Docking_Results_{timestamp}.pse"
    if pymol_enabled:
        try:
            cmd.save(str(session_path))
        except Exception as exc:
            print(f"⚠️ Failed to save PyMOL session: {exc}")

    audit_csv = out_dir_pdb / f"top_hits_audit_{timestamp}.csv"
    write_audit_csv(audit_csv, audit_rows)

    print("\n🎉 Done!")
    if pymol_enabled:
        print(f"🧪 PyMOL session: {session_path}")
    print(f"📦 PDB complexes: {out_dir_pdb}")
    print(f"📦 Corrected SDFs: {out_dir_sdf}")
    print(f"📦 Raw selected PDBQT: {out_dir_pdbqt}")
    print(f"🧾 Audit CSV: {audit_csv}")


if __name__ == "__main__":
    main()
