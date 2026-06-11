from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.rdmolops import GetFormalCharge
    HAVE_RDKIT = True
except Exception:  # pragma: no cover - depends on optional runtime
    Chem = Any  # type: ignore[assignment]
    Descriptors = None  # type: ignore[assignment]
    rdMolDescriptors = None  # type: ignore[assignment]
    GetFormalCharge = None  # type: ignore[assignment]
    HAVE_RDKIT = False

from ligand_naming import parse_ligand_variant


CHEMICAL_METADATA_COLUMNS = [
    "CanonicalSMILES",
    "IsomericSMILES",
    "StateCanonicalSMILES",
    "StateIsomericSMILES",
    "Formula",
    "FormalCharge",
    "NetCharge",
    "NumAtoms",
    "NumHeavyAtoms",
    "NumRotatableBonds",
    "ExactMolWt",
    "SourceInput",
    "SourceInputType",
    "SourceLigandID",
    "SourceRecordIndex",
    "SourceRecordName",
    "SDFTitle",
    "OriginalMolName",
]

MANIFEST_COLUMNS = [
    "LigandBase",
    "LigandVariant",
    "SourceLigandID",
    "SourceInput",
    "SourceInputType",
    "SourceRecordIndex",
    "SourceRecordName",
    "SDFTitle",
    "OriginalMolName",
    "ProtomerTag",
    "ProtomerIndex",
    "TautomerTag",
    "TautomerIndex",
    "ConformerTag",
    "ConformerIndex",
    "StateTag",
    "CanonicalSMILES",
    "IsomericSMILES",
    "StateCanonicalSMILES",
    "StateIsomericSMILES",
    "Formula",
    "FormalCharge",
    "NetCharge",
    "NumAtoms",
    "NumHeavyAtoms",
    "NumRotatableBonds",
    "ExactMolWt",
    "PDBQTFile",
    "SDFFile",
    "TmpSDFFile",
    "GenerationStatus",
    "GenerationWarning",
]

SCORE_COLUMNS = [
    "LigandBase",
    "LigandVariant",
    "ProtomerTag",
    "ProtomerIndex",
    "TautomerTag",
    "TautomerIndex",
    "ConformerTag",
    "ConformerIndex",
    "StateTag",
    "LegacyPoseTag",
] + CHEMICAL_METADATA_COLUMNS


def _safe_smiles(mol: Optional[Chem.Mol], isomeric: bool = False) -> str:
    if mol is None:
        return ""
    if not HAVE_RDKIT:
        return ""
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=isomeric)
    except Exception:
        try:
            return Chem.MolToSmiles(Chem.Mol(mol), isomericSmiles=isomeric)
        except Exception:
            return ""


def _safe_descriptor(func, mol: Optional[Chem.Mol]) -> Any:
    if mol is None:
        return ""
    if func is None:
        return ""
    if not HAVE_RDKIT:
        return ""
    try:
        return func(mol)
    except Exception:
        return ""


def _stringify_manifest_row(row: Dict[str, Any]) -> Dict[str, Any]:
    cooked = {}
    for key in MANIFEST_COLUMNS:
        value = row.get(key, "")
        if value is None:
            cooked[key] = ""
        else:
            cooked[key] = value
    return cooked


def build_ligand_state_metadata(
    ligand_base: str,
    ligand_variant: str,
    source_ligand_id: str,
    source_input: str,
    source_input_type: str,
    source_record_index: Any = "",
    source_record_name: str = "",
    sdf_title: str = "",
    original_mol_name: str = "",
    canonical_mol: Optional[Chem.Mol] = None,
    state_mol: Optional[Chem.Mol] = None,
    pdbqt_file: str = "",
    sdf_file: str = "",
    tmp_sdf_file: str = "",
    generation_status: str = "ok",
    generation_warning: str = "",
) -> Dict[str, Any]:
    parsed = parse_ligand_variant(ligand_variant)
    canonical_smiles = _safe_smiles(canonical_mol, isomeric=False)
    isomeric_smiles = _safe_smiles(canonical_mol, isomeric=True)
    state_canonical_smiles = _safe_smiles(state_mol, isomeric=False)
    state_isomeric_smiles = _safe_smiles(state_mol, isomeric=True)
    formula = _safe_descriptor(rdMolDescriptors.CalcMolFormula if rdMolDescriptors else None, state_mol)
    formal_charge = _safe_descriptor(GetFormalCharge, state_mol)
    exact_mw = _safe_descriptor(Descriptors.ExactMolWt if Descriptors else None, state_mol)
    num_atoms = _safe_descriptor(lambda m: m.GetNumAtoms(), state_mol)
    num_heavy_atoms = _safe_descriptor(lambda m: m.GetNumHeavyAtoms(), state_mol)
    num_rotatable_bonds = _safe_descriptor(rdMolDescriptors.CalcNumRotatableBonds if rdMolDescriptors else None, state_mol)
    row = {
        "LigandBase": ligand_base or parsed["LigandBase"],
        "LigandVariant": ligand_variant,
        "SourceLigandID": source_ligand_id or ligand_base,
        "SourceInput": source_input,
        "SourceInputType": source_input_type,
        "SourceRecordIndex": source_record_index,
        "SourceRecordName": source_record_name,
        "SDFTitle": sdf_title or source_record_name,
        "OriginalMolName": original_mol_name or source_record_name or sdf_title,
        "ProtomerTag": parsed["ProtomerTag"],
        "ProtomerIndex": parsed["ProtomerIndex"],
        "TautomerTag": parsed["TautomerTag"],
        "TautomerIndex": parsed["TautomerIndex"],
        "ConformerTag": parsed["ConformerTag"],
        "ConformerIndex": parsed["ConformerIndex"],
        "StateTag": parsed["StateTag"],
        "CanonicalSMILES": canonical_smiles,
        "IsomericSMILES": isomeric_smiles,
        "StateCanonicalSMILES": state_canonical_smiles,
        "StateIsomericSMILES": state_isomeric_smiles,
        "Formula": formula,
        "FormalCharge": formal_charge,
        "NetCharge": formal_charge,
        "NumAtoms": num_atoms,
        "NumHeavyAtoms": num_heavy_atoms,
        "NumRotatableBonds": num_rotatable_bonds,
        "ExactMolWt": exact_mw,
        "PDBQTFile": pdbqt_file,
        "SDFFile": sdf_file,
        "TmpSDFFile": tmp_sdf_file,
        "GenerationStatus": generation_status,
        "GenerationWarning": generation_warning,
    }
    return _stringify_manifest_row(row)


def write_ligand_state_manifest(path: Path, rows: Iterable[Dict[str, Any]]) -> Path:
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(_stringify_manifest_row(row))
    return path


def load_ligand_state_manifest(path: Path) -> Dict[str, Dict[str, str]]:
    manifest = {}
    path = Path(path)
    if not path.exists():
        return manifest
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                variant = (row.get("LigandVariant") or "").strip()
                if not variant:
                    continue
                manifest[variant] = row
    except Exception:
        return {}
    return manifest


def find_ligand_state_manifests(search_roots: Iterable[Path]) -> List[Path]:
    found: List[Path] = []
    seen = set()
    for root in search_roots:
        if root is None:
            continue
        root = Path(root)
        if not root.exists():
            continue
        if root.is_file():
            candidate = root if root.name == "ligand_state_manifest.csv" else None
            if candidate and str(candidate.resolve()) not in seen:
                seen.add(str(candidate.resolve()))
                found.append(candidate.resolve())
            continue
        direct = root / "ligand_state_manifest.csv"
        if direct.exists():
            resolved = direct.resolve()
            if str(resolved) not in seen:
                seen.add(str(resolved))
                found.append(resolved)
        for candidate in root.rglob("ligand_state_manifest.csv"):
            resolved = candidate.resolve()
            if str(resolved) in seen:
                continue
            seen.add(str(resolved))
            found.append(resolved)
    return found


def merge_ligand_metadata(ligand_variant: str, base_row: Optional[Dict[str, Any]] = None, manifest_row: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    parsed = parse_ligand_variant(ligand_variant)
    row = dict(base_row or {})
    manifest_row = manifest_row or {}
    row["LigandVariant"] = row.get("LigandVariant") or manifest_row.get("LigandVariant") or parsed["LigandVariant"]
    row["LigandBase"] = row.get("LigandBase") or manifest_row.get("LigandBase") or parsed["LigandBase"]
    row["Ligand"] = row.get("Ligand") or row["LigandBase"]
    for key in ("ProtomerTag", "TautomerTag", "ConformerTag", "StateTag", "LegacyPoseTag"):
        row[key] = row.get(key) or manifest_row.get(key) or parsed.get(key, "")
    for key in ("ProtomerIndex", "TautomerIndex", "ConformerIndex"):
        row[key] = row.get(key) if row.get(key) not in ("", None) else manifest_row.get(key, parsed.get(key))
    for key in CHEMICAL_METADATA_COLUMNS:
        row[key] = row.get(key) or manifest_row.get(key, "")
    return row
