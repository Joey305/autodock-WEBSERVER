from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shlex
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple


class CenterResolutionError(Exception):
    def __init__(
        self,
        error: str,
        message: str,
        details: Optional[Dict[str, Any]] = None,
        status_code: int = 400,
    ):
        super().__init__(message)
        self.error = error
        self.message = message
        self.details = details or {}
        self.status_code = status_code

    def to_dict(self) -> Dict[str, Any]:
        return {"error": self.error, "message": self.message, "details": self.details}


@dataclass(frozen=True)
class StructureAtom:
    record: str
    atom_name: str
    resname: str
    chain: str
    resi: str
    insertion_code: str
    x: float
    y: float
    z: float
    element: str

    @property
    def instance_key(self) -> Tuple[str, str, str, str, str]:
        return (self.record, self.resname, self.chain, self.resi, self.insertion_code)


SUPPORTED_STRUCTURE_SUFFIXES = {".pdb", ".ent", ".pdbqt", ".cif", ".mmcif"}


def parse_pdb_atoms(path: Path) -> List[StructureAtom]:
    suffix = path.suffix.lower()
    if suffix not in SUPPORTED_STRUCTURE_SUFFIXES:
        raise CenterResolutionError(
            "unsupported_structure_format",
            f"Headless center resolution currently supports PDB, ENT, PDBQT, and mmCIF-like coordinate files; got {path.name}.",
            {"supported_suffixes": sorted(SUPPORTED_STRUCTURE_SUFFIXES), "path": str(path)},
            status_code=415,
        )
    if suffix in {".cif", ".mmcif"}:
        return parse_mmcif_atoms(path)
    atoms: List[StructureAtom] = []
    with path.open("r", errors="replace") as handle:
        for line in handle:
            atom = parse_pdb_atom_line(line)
            if atom is not None:
                atoms.append(atom)
    return atoms


def parse_mmcif_atoms(path: Path) -> List[StructureAtom]:
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    atoms: List[StructureAtom] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        headers: List[str] = []
        j = i + 1
        while j < len(lines) and lines[j].strip().startswith("_atom_site."):
            headers.append(lines[j].strip())
            j += 1
        if not headers or "_atom_site.group_PDB" not in headers:
            i += 1
            continue
        atoms.extend(_parse_mmcif_atom_loop(lines, j, headers))
        if atoms:
            return atoms
        i = j
    return atoms


def _parse_mmcif_atom_loop(lines: Sequence[str], start_idx: int, headers: Sequence[str]) -> List[StructureAtom]:
    idx = {header: pos for pos, header in enumerate(headers)}
    required = [
        "_atom_site.group_PDB",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]
    if any(name not in idx for name in required):
        return []
    atoms: List[StructureAtom] = []
    row_tokens: List[str] = []
    expected_cols = len(headers)
    for raw_line in lines[start_idx:]:
        stripped = raw_line.strip()
        if not stripped or stripped == "#":
            if stripped == "#":
                break
            continue
        if stripped == "loop_" or stripped.startswith("_"):
            break
        if stripped.startswith(";"):
            continue
        row_tokens.extend(_mmcif_split_tokens(stripped))
        while len(row_tokens) >= expected_cols:
            row = row_tokens[:expected_cols]
            row_tokens = row_tokens[expected_cols:]
            atom = _structure_atom_from_mmcif_row(row, idx)
            if atom is not None:
                atoms.append(atom)
    return atoms


def _mmcif_split_tokens(line: str) -> List[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def _mmcif_value(row: Sequence[str], idx: Dict[str, int], key: str) -> str:
    pos = idx.get(key)
    if pos is None or pos >= len(row):
        return ""
    value = row[pos].strip()
    if value in {"?", "."}:
        return ""
    return value


def _structure_atom_from_mmcif_row(row: Sequence[str], idx: Dict[str, int]) -> Optional[StructureAtom]:
    record = _mmcif_value(row, idx, "_atom_site.group_PDB").upper()
    if record not in {"ATOM", "HETATM"}:
        return None
    atom_name = _mmcif_value(row, idx, "_atom_site.auth_atom_id") or _mmcif_value(row, idx, "_atom_site.label_atom_id")
    resname = (_mmcif_value(row, idx, "_atom_site.auth_comp_id") or _mmcif_value(row, idx, "_atom_site.label_comp_id")).upper()
    chain = (_mmcif_value(row, idx, "_atom_site.auth_asym_id") or _mmcif_value(row, idx, "_atom_site.label_asym_id")).upper()
    resi = _mmcif_value(row, idx, "_atom_site.auth_seq_id") or _mmcif_value(row, idx, "_atom_site.label_seq_id")
    insertion_code = _mmcif_value(row, idx, "_atom_site.pdbx_PDB_ins_code").upper()
    element = _mmcif_value(row, idx, "_atom_site.type_symbol").upper()
    try:
        x = float(_mmcif_value(row, idx, "_atom_site.Cartn_x"))
        y = float(_mmcif_value(row, idx, "_atom_site.Cartn_y"))
        z = float(_mmcif_value(row, idx, "_atom_site.Cartn_z"))
    except Exception:
        return None
    return StructureAtom(
        record=record,
        atom_name=atom_name.upper(),
        resname=resname,
        chain=chain,
        resi=resi,
        insertion_code=insertion_code,
        x=x,
        y=y,
        z=z,
        element=element,
    )


def parse_pdb_atom_line(line: str) -> Optional[StructureAtom]:
    record = line[:6].strip().upper()
    if record not in {"ATOM", "HETATM"}:
        return None
    try:
        return StructureAtom(
            record=record,
            atom_name=line[12:16].strip(),
            resname=line[17:20].strip().upper(),
            chain=line[21].strip().upper(),
            resi=line[22:26].strip(),
            insertion_code=line[26].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            element=line[76:78].strip().upper() if len(line) >= 78 else "",
        )
    except Exception:
        return None


def resolve_center_from_file(path: Path, payload: Dict[str, Any]) -> Dict[str, Any]:
    method = normalize_method(payload.get("method"))
    if method == "xyz":
        return resolve_xyz(payload)
    atoms = parse_pdb_atoms(path)
    if not atoms:
        raise CenterResolutionError(
            "no_atoms_found",
            f"No ATOM or HETATM coordinate records were found in receptor {path.name}.",
            {"receptor": path.name},
        )
    return resolve_center_from_atoms(atoms, payload, receptor_name=path.name)


def resolve_center_from_atoms(
    atoms: Sequence[StructureAtom],
    payload: Dict[str, Any],
    receptor_name: str = "",
) -> Dict[str, Any]:
    method = normalize_method(payload.get("method"))
    if method == "xyz":
        return resolve_xyz(payload)
    if method == "residue":
        return resolve_residue(atoms, payload, receptor_name)
    if method == "hetatm":
        return resolve_hetatm(atoms, payload, receptor_name)
    if method == "atom":
        return resolve_atom(atoms, payload, receptor_name)
    if method == "selection":
        selection = payload.get("selection")
        if not isinstance(selection, dict):
            raise CenterResolutionError("bad_request", "Selection mode requires a selection object.")
        return resolve_selection(atoms, selection, payload, receptor_name)
    raise CenterResolutionError(
        "unsupported_center_method",
        f"Unsupported center resolution method: {payload.get('method')!r}.",
        {"supported_methods": ["xyz", "residue", "hetatm", "atom", "selection"]},
    )


def resolve_xyz(payload: Dict[str, Any]) -> Dict[str, Any]:
    center = payload.get("center")
    if not isinstance(center, (list, tuple)) or len(center) != 3:
        raise CenterResolutionError("invalid_xyz", "XYZ center must be an array of exactly three numeric values.")
    try:
        coords = [float(center[0]), float(center[1]), float(center[2])]
    except Exception as exc:
        raise CenterResolutionError("invalid_xyz", "XYZ center values must be numeric.") from exc
    return build_result("xyz", coords, payload, matched={"atom_count": 0})


def resolve_residue(atoms: Sequence[StructureAtom], payload: Dict[str, Any], receptor_name: str) -> Dict[str, Any]:
    chain = norm_chain(payload.get("chain"))
    resi = norm_resi(payload.get("resi"))
    resname = norm_resname(payload.get("resname"))
    if not chain or not resi:
        raise CenterResolutionError("bad_request", "Residue center requires chain and resi.")
    selected = [
        atom for atom in atoms
        if atom.record == "ATOM"
        and atom.chain == chain
        and atom.resi == resi
        and (not resname or atom.resname == resname)
    ]
    if not selected:
        return no_match("residue", receptor_name, {"chain": chain, "resi": resi, "resname": resname})
    groups = group_instances(selected)
    if len(groups) > 1:
        raise_ambiguous("Multiple residue instances matched. Provide resname or insertion code.", groups)
    return centroid_result("residue", selected, payload, receptor_name)


def resolve_hetatm(atoms: Sequence[StructureAtom], payload: Dict[str, Any], receptor_name: str) -> Dict[str, Any]:
    het = norm_resname(payload.get("het") or payload.get("hetatm") or payload.get("ligand") or payload.get("resname"))
    chain = norm_chain(payload.get("chain"))
    resi = norm_resi(payload.get("resi"))
    if not het:
        raise CenterResolutionError("bad_request", "HETATM center requires het, hetatm, ligand, or resname.")
    selected = [
        atom for atom in atoms
        if atom.record == "HETATM"
        and atom.resname == het
        and (not chain or atom.chain == chain)
        and (not resi or atom.resi == resi)
    ]
    if not selected:
        return no_match("hetatm", receptor_name, {"het": het, "chain": chain, "resi": resi})
    groups = group_instances(selected)
    if len(groups) > 1 and (not chain or not resi):
        raise_ambiguous(f"Multiple {het} instances matched. Provide chain and resi.", groups)
    return centroid_result("hetatm", selected, payload, receptor_name)


def resolve_atom(atoms: Sequence[StructureAtom], payload: Dict[str, Any], receptor_name: str) -> Dict[str, Any]:
    atom_name = (payload.get("atom_name") or payload.get("atom") or "").strip().upper()
    if not atom_name:
        raise CenterResolutionError("bad_request", "Atom center requires atom_name.")
    selection = {
        "record": payload.get("record"),
        "resname": payload.get("resname") or payload.get("het") or payload.get("ligand"),
        "chain": payload.get("chain"),
        "resi": payload.get("resi"),
        "atom_name": atom_name,
    }
    selected = filter_atoms(atoms, selection)
    if not selected:
        return no_match("atom", receptor_name, selection)
    if len(selected) > 1:
        raise CenterResolutionError(
            "ambiguous_selection",
            "Multiple atoms matched. Provide record, resname, chain, and resi to disambiguate.",
            {"candidates": [atom_metadata(atom) for atom in selected[:25]], "match_count": len(selected)},
        )
    atom = selected[0]
    return build_result(
        "atom",
        [atom.x, atom.y, atom.z],
        payload,
        receptor_name,
        matched={**instance_metadata(selected), "atom_count": 1, "atom_names": [atom.atom_name]},
    )


def resolve_selection(
    atoms: Sequence[StructureAtom],
    selection: Dict[str, Any],
    payload: Dict[str, Any],
    receptor_name: str,
) -> Dict[str, Any]:
    selected = filter_atoms(atoms, selection)
    if not selected:
        return no_match("selection", receptor_name, selection)
    groups = group_instances(selected)
    if len(groups) > 1 and not payload.get("allow_multiple"):
        raise_ambiguous("Selection matched multiple residue or ligand instances. Set allow_multiple=true to centroid all matches.", groups)
    return centroid_result("selection", selected, payload, receptor_name)


def filter_atoms(atoms: Sequence[StructureAtom], selection: Dict[str, Any]) -> List[StructureAtom]:
    record = (selection.get("record") or "").strip().upper()
    resname = norm_resname(selection.get("resname") or selection.get("het") or selection.get("ligand"))
    chain = norm_chain(selection.get("chain"))
    resi = norm_resi(selection.get("resi"))
    atom_name = (selection.get("atom_name") or selection.get("atom") or "").strip().upper()
    insertion_code = (selection.get("insertion_code") or selection.get("icode") or "").strip().upper()
    return [
        atom for atom in atoms
        if (not record or atom.record == record)
        and (not resname or atom.resname == resname)
        and (not chain or atom.chain == chain)
        and (not resi or atom.resi == resi)
        and (not atom_name or atom.atom_name.upper() == atom_name)
        and (not insertion_code or atom.insertion_code == insertion_code)
    ]


def centroid_result(method: str, atoms: Sequence[StructureAtom], payload: Dict[str, Any], receptor_name: str) -> Dict[str, Any]:
    coords = [
        sum(atom.x for atom in atoms) / len(atoms),
        sum(atom.y for atom in atoms) / len(atoms),
        sum(atom.z for atom in atoms) / len(atoms),
    ]
    return build_result(
        method,
        coords,
        payload,
        receptor_name,
        matched={**instance_metadata(atoms), "atom_count": len(atoms), "atom_names": sorted({atom.atom_name for atom in atoms})},
    )


def build_result(
    method: str,
    center: Sequence[float],
    payload: Dict[str, Any],
    receptor_name: str = "",
    matched: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    coords = [round(float(v), 6) for v in center]
    size = float(payload.get("size") or 20)
    result = {
        "method": method,
        "receptor": receptor_name,
        "center": coords,
        "size": size,
        "matched": matched or {},
        "save_payload": {"method": "xyz", "center": coords, "size": size},
    }
    return result


def no_match(method: str, receptor_name: str, filters: Dict[str, Any]) -> Dict[str, Any]:
    readable = " ".join(f"{key}={value}" for key, value in filters.items() if value)
    raise CenterResolutionError(
        "no_atoms_matched",
        f"No atoms matched {readable} in receptor {receptor_name}.".strip(),
        {"method": method, "receptor": receptor_name, "filters": filters},
        status_code=404,
    )


def raise_ambiguous(message: str, groups: Dict[Tuple[str, str, str, str, str], List[StructureAtom]]) -> None:
    candidates = []
    for group_atoms in groups.values():
        candidates.append({**instance_metadata(group_atoms), "atom_count": len(group_atoms)})
    raise CenterResolutionError(
        "ambiguous_selection",
        message,
        {"candidates": sorted(candidates, key=lambda c: (c.get("resname", ""), c.get("chain", ""), c.get("resi", "")))},
    )


def group_instances(atoms: Iterable[StructureAtom]) -> Dict[Tuple[str, str, str, str, str], List[StructureAtom]]:
    groups: Dict[Tuple[str, str, str, str, str], List[StructureAtom]] = {}
    for atom in atoms:
        groups.setdefault(atom.instance_key, []).append(atom)
    return groups


def instance_metadata(atoms: Sequence[StructureAtom]) -> Dict[str, Any]:
    first = atoms[0]
    records = sorted({atom.record for atom in atoms})
    return {
        "record": records[0] if len(records) == 1 else records,
        "resname": first.resname,
        "chain": first.chain,
        "resi": first.resi,
        "insertion_code": first.insertion_code,
    }


def atom_metadata(atom: StructureAtom) -> Dict[str, Any]:
    return {
        "record": atom.record,
        "atom_name": atom.atom_name,
        "resname": atom.resname,
        "chain": atom.chain,
        "resi": atom.resi,
        "insertion_code": atom.insertion_code,
        "center": [atom.x, atom.y, atom.z],
    }


def normalize_method(method: Any) -> str:
    value = (method or "xyz").strip().lower()
    aliases = {"ligand": "hetatm", "het": "hetatm", "hetatm": "hetatm"}
    return aliases.get(value, value)


def norm_resname(value: Any) -> str:
    return ("" if value is None else str(value)).strip().upper()


def norm_chain(value: Any) -> str:
    return ("" if value is None else str(value)).strip().upper()


def norm_resi(value: Any) -> str:
    return ("" if value is None else str(value)).strip()
