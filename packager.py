# ==============================
# packager.py
# ==============================
from __future__ import annotations

import csv
import io
import json
import os
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

RUNTIME_ROOT_FILES = [
    "0_LIGSPLIT.py",
    "1_ConformerGeneration.py",
    "3_Complete_batch_docking.py",
    "4_ParseScores.py",
    "4C_ConcatenateScores.py",
    "5_BuidlHTMLViz.py",
    "5_CompactedHTMLViz.py",
    "5_COMPACTED_SDF_HTML.py",
    "5C_BuildPymolSesh.py",
    "6_MDpymacs.py",
    "7_Graphs.py",
    "create_vina_env.sh",
    "docking.yaml",
    "ligand_naming.py",
    "ligand_manifest.py",
]

RUNTIME_LSF_FILES = [
    "1B_confgen_batch.py",
    "3B_ServerDocks.py",
    "3a_PDB2PDBQTbatch.py",
    "4B_LSFbatch.py",
    "hpc_profiles.py",
    "lsf_templates.py",
    "runDOCKING-tmux.sh",
]

RUNTIME_TOOL_DIRS = [
    "AutoDockTools_py3",
    "viewer_templates",
]

RUNTIME_SOURCE_ALIASES = {
    "7_Graphs.py": ["7_Graph.py"],
    "3a_PDB2PDBQTbatch.py": ["2a_PDB2PDBQTbatch.py"],
}

PORTABLE_PACKAGE_README = "README_RUN_LOCAL.md"
SUPPORTED_LIGAND_SUFFIXES = {".sdf", ".smiles", ".smi", ".csv"}
LIGAND_MANIFEST_NAME = "ligand_state_manifest.csv"


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
    with zipfile.ZipFile(io.BytesIO(buf)) as zf:
        zf.extractall(dest_dir)
    return str(dest_dir)


def _is_hidden_or_macos_name(name: str) -> bool:
    return (
        not name
        or name in {"__MACOSX", ".DS_Store"}
        or name.startswith("._")
        or name.startswith(".")
    )


def _is_hidden_or_macos_path(rel_path: Path) -> bool:
    return any(_is_hidden_or_macos_name(part) for part in rel_path.parts)


def _is_generated_ligand_dir(name: str) -> bool:
    return (
        name.startswith("Ligands_TMP_SDF_")
        or (name.startswith("Ligands_") and "_PDBQT_" in name)
    )


def _safe_zip_member_path(member_name: str) -> Optional[Path]:
    name = (member_name or "").replace("\\", "/").strip("/")
    if not name:
        return None
    rel = Path(name)
    if rel.is_absolute() or any(part == ".." for part in rel.parts):
        return None
    return rel


def _resolve_unique_ligand_name(filename: str, seen_names: Dict[str, int]) -> tuple[str, Optional[str]]:
    path = Path(filename)
    stem = path.stem or "ligand"
    suffix = path.suffix
    candidate = f"{stem}{suffix}"
    key = candidate.lower()
    count = seen_names.get(key, 0) + 1
    seen_names[key] = count
    if count == 1:
        return candidate, None
    renamed = f"{stem}_{count:03d}{suffix}"
    return renamed, f"Renamed duplicate ligand file {candidate} to {renamed}"


def _copy_ligand_bytes(
    dest_dir: Path,
    dest_name: str,
    payload: bytes,
):
    out = dest_dir / dest_name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(payload)


def _normalize_ligand_upload_entries(
    entries: Iterable[tuple[Path, bytes]],
    dest_dir: Path,
) -> Dict[str, Any]:
    accepted_files: list[str] = []
    ignored_files: list[str] = []
    warnings: list[str] = []
    filetypes: set[str] = set()
    seen_names: Dict[str, int] = {}

    dest_dir.mkdir(parents=True, exist_ok=True)

    for rel_path, payload in entries:
        if _is_hidden_or_macos_path(rel_path):
            ignored_files.append(rel_path.as_posix())
            continue
        suffix = rel_path.suffix.lower()
        if suffix not in SUPPORTED_LIGAND_SUFFIXES:
            ignored_files.append(rel_path.as_posix())
            warnings.append(f"Ignored unsupported ligand file {rel_path.name}")
            continue
        flattened = len(rel_path.parts) > 1
        dest_name, rename_warning = _resolve_unique_ligand_name(rel_path.name, seen_names)
        if flattened:
            warnings.append(f"Flattened nested ligand path {rel_path.as_posix()} to Ligands/{dest_name}")
        if rename_warning:
            warnings.append(rename_warning)
        _copy_ligand_bytes(dest_dir, dest_name, payload)
        accepted_files.append(dest_name)
        filetypes.add(suffix)

    return {
        "accepted_files": accepted_files,
        "accepted_count": len(accepted_files),
        "ignored_files": ignored_files,
        "warnings": warnings,
        "filetypes": sorted(filetypes),
        "is_csv": len(filetypes) == 1 and ".csv" in filetypes and len(accepted_files) == 1,
        "ligands_root": str(dest_dir),
    }


def save_uploaded_ligand_zip(file_storage, dest_dir: Path) -> Dict[str, Any]:
    entries: list[tuple[Path, bytes]] = []
    buf = file_storage.read()
    with zipfile.ZipFile(io.BytesIO(buf)) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            rel = _safe_zip_member_path(info.filename)
            if rel is None:
                continue
            if (info.external_attr >> 16) & 0o170000 == 0o120000:
                continue
            entries.append((rel, zf.read(info)))
    return _normalize_ligand_upload_entries(entries, dest_dir)


def save_uploaded_ligand_folder(files, dest_dir: Path) -> Dict[str, Any]:
    entries: list[tuple[Path, bytes]] = []
    for file_storage in files:
        rel = _safe_zip_member_path(file_storage.filename or "")
        if rel is None:
            continue
        entries.append((rel, file_storage.read()))
    return _normalize_ligand_upload_entries(entries, dest_dir)


def normalize_ligand_tree(src_dir: Path, dest_dir: Path) -> Dict[str, Any]:
    accepted_files: list[str] = []
    ignored_files: list[str] = []
    warnings: list[str] = []
    filetypes: set[str] = set()
    seen_names: Dict[str, int] = {}
    copied_manifest = False

    dest_dir.mkdir(parents=True, exist_ok=True)

    for root, dirnames, filenames in os.walk(src_dir):
        root_path = Path(root)
        rel_root = root_path.relative_to(src_dir)
        dirnames[:] = [
            d for d in dirnames
            if not _is_hidden_or_macos_name(d) and not _is_generated_ligand_dir(d)
        ]
        for filename in sorted(filenames):
            rel_path = rel_root / filename
            if _is_hidden_or_macos_path(rel_path):
                ignored_files.append(rel_path.as_posix())
                continue
            src_path = root_path / filename
            if filename == LIGAND_MANIFEST_NAME:
                if copied_manifest:
                    warnings.append(f"Ignored duplicate ligand manifest {rel_path.as_posix()}")
                    continue
                _safe_copy(src_path, dest_dir / LIGAND_MANIFEST_NAME)
                copied_manifest = True
                continue
            suffix = src_path.suffix.lower()
            if suffix not in SUPPORTED_LIGAND_SUFFIXES:
                ignored_files.append(rel_path.as_posix())
                continue
            dest_name, rename_warning = _resolve_unique_ligand_name(src_path.name, seen_names)
            if rel_root != Path("."):
                warnings.append(f"Flattened nested ligand path {rel_path.as_posix()} to Ligands/{dest_name}")
            if rename_warning:
                warnings.append(rename_warning)
            _safe_copy(src_path, dest_dir / dest_name)
            accepted_files.append(dest_name)
            filetypes.add(suffix)

    return {
        "accepted_files": accepted_files,
        "accepted_count": len(accepted_files),
        "ignored_files": ignored_files,
        "warnings": warnings,
        "filetypes": sorted(filetypes),
        "is_csv": len(filetypes) == 1 and ".csv" in filetypes and len(accepted_files) == 1,
        "ligands_root": str(dest_dir),
    }


def move_dir(src: Path, dst: Path):
    if dst.exists():
        shutil.rmtree(dst)
    shutil.move(str(src), str(dst))


def _read_any_centers(csv_path: Path) -> list[dict]:
    rows: list[dict] = []
    if not csv_path.exists():
        return rows
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        fields = {(h or "").strip().upper(): h for h in (reader.fieldnames or [])}
        has_new = all(k in fields for k in ("PDB_ID", "X", "Y", "Z"))
        has_old = all(k in fields for k in ("RECEPTOR_PDBQT", "CENTER_X", "CENTER_Y", "CENTER_Z"))
        for row in reader:
            try:
                if has_new:
                    key = (row.get(fields["PDB_ID"]) or "").strip()
                    x = float(row.get(fields["X"], "0"))
                    y = float(row.get(fields["Y"], "0"))
                    z = float(row.get(fields["Z"], "0"))
                    s = float(row.get(fields.get("SIZE", "SIZE"), row.get("SIZE", "20")))
                elif has_old:
                    key = (row.get(fields["RECEPTOR_PDBQT"]) or "").strip()
                    x = float(row.get(fields["CENTER_X"], "0"))
                    y = float(row.get(fields["CENTER_Y"], "0"))
                    z = float(row.get(fields["CENTER_Z"], "0"))
                    s = float(row.get(fields.get("SIZE", "SIZE"), row.get("size", "20")))
                else:
                    continue
                if key:
                    rows.append({"PDB_ID": key, "X": x, "Y": y, "Z": z, "SIZE": s})
            except Exception:
                continue
    return rows


def _write_canonical_centers(csv_path: Path, rows: list[dict]):
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["PDB_ID", "X", "Y", "Z", "SIZE"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "PDB_ID": row["PDB_ID"],
                    "X": float(row["X"]),
                    "Y": float(row["Y"]),
                    "Z": float(row["Z"]),
                    "SIZE": float(row.get("SIZE", 20.0)),
                }
            )


def write_centers_csv_row(csv_path: Path, pdbqt_name: str, center_xyz, size: float):
    csv_path = Path(csv_path)
    rows = _read_any_centers(csv_path)
    cx, cy, cz = map(float, center_xyz)

    found = False
    for row in rows:
        if row["PDB_ID"] == pdbqt_name:
            row.update({"X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})
            found = True
            break
    if not found:
        rows.append({"PDB_ID": pdbqt_name, "X": cx, "Y": cy, "Z": cz, "SIZE": float(size)})

    _write_canonical_centers(csv_path, rows)


def normalize_centers_csv_to_canonical(src_csv: Path, dst_csv: Optional[Path] = None) -> Path:
    src_csv = Path(src_csv)
    dst_csv = Path(dst_csv) if dst_csv else src_csv
    rows = _read_any_centers(src_csv)
    if not rows and not src_csv.exists():
        _write_canonical_centers(dst_csv, [])
    else:
        _write_canonical_centers(dst_csv, rows)
    return dst_csv


def _safe_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    if dst.suffix in {".sh", ".py", ".lsf"}:
        try:
            mode = os.stat(dst).st_mode
            os.chmod(dst, mode | 0o111)
        except Exception:
            pass


def _resolve_runtime_source(repo_root: Path, requested_name: str) -> tuple[Optional[Path], str]:
    direct = repo_root / requested_name
    if direct.exists():
        return direct, requested_name
    for alias in RUNTIME_SOURCE_ALIASES.get(requested_name, []):
        src = repo_root / alias
        if src.exists():
            return src, requested_name
    return None, requested_name


def _stage_runtime_file_group(
    jobroot: Path,
    repo_root: Path,
    requested_names: list[str],
    warnings: list[str],
):
    for requested_name in requested_names:
        src, dst_name = _resolve_runtime_source(repo_root, requested_name)
        if src is None:
            warnings.append(f"Runtime file not found: {requested_name}")
            continue
        _safe_copy(src, jobroot / dst_name)


def _stage_runtime_dirs(jobroot: Path, repo_root: Path, warnings: list[str]):
    for rel in RUNTIME_TOOL_DIRS:
        src_dir = repo_root / rel
        if not src_dir.exists() or not src_dir.is_dir():
            warnings.append(f"Runtime directory not found: {rel}")
            continue
        dst_dir = jobroot / src_dir.name
        if dst_dir.exists():
            shutil.rmtree(dst_dir)
        shutil.copytree(src_dir, dst_dir)


def stage_runtime_assets(jobroot: Path, repo_root: Path, package_mode: str) -> list[str]:
    warnings: list[str] = []
    _stage_runtime_file_group(jobroot, repo_root, RUNTIME_ROOT_FILES, warnings)
    if package_mode in {"lsf", "joey_lsf", "mainak_lsf", "custom_lsf"}:
        _stage_runtime_file_group(jobroot, repo_root, RUNTIME_LSF_FILES, warnings)
    _stage_runtime_dirs(jobroot, repo_root, warnings)
    return warnings


def _copy_prepared_receptors(ws: Path, jobroot: Path) -> bool:
    prepared_sources = [ws / "Receptors_PDBQT", ws / "Receptors_PDBQT_Converted"]
    source = None
    for candidate in prepared_sources:
        if candidate.exists() and any(candidate.glob("*.pdbqt")):
            source = candidate
            break

    if source is None:
        raw_dir = ws / "Receptors"
        if raw_dir.exists():
            raw_pdbqts = list(raw_dir.glob("*.pdbqt"))
            if raw_pdbqts:
                target = jobroot / "Receptors"
                target.mkdir(parents=True, exist_ok=True)
                for src in raw_pdbqts:
                    _safe_copy(src, target / src.name)
                return True
        return False

    target = jobroot / "Receptors"
    target.mkdir(parents=True, exist_ok=True)
    for src in source.glob("*.pdbqt"):
        dst_name = src.name.replace(".converted.", ".")
        _safe_copy(src, target / dst_name)
    return True


def _copy_raw_receptors(ws: Path, jobroot: Path):
    raw_dir = ws / "Receptors"
    if raw_dir.exists():
        shutil.copytree(raw_dir, jobroot / "Receptors_PDB", dirs_exist_ok=True)
        return

    raw_snapshot = ws / "Receptors_PDB"
    if raw_snapshot.exists():
        shutil.copytree(raw_snapshot, jobroot / "Receptors_PDB", dirs_exist_ok=True)


def _chmod_scripts(root: Path):
    for path in root.rglob("*"):
        if path.is_file() and path.suffix in {".lsf", ".sh", ".py"}:
            try:
                mode = os.stat(path).st_mode
                os.chmod(path, mode | 0o111)
            except Exception:
                pass


def assemble_job_tree(ws: Path, rec_dir: Path, lig_dir: Path, package_mode: str = "portable") -> tuple[Path, list[str]]:
    jobroot = ws / "job"
    if jobroot.exists():
        shutil.rmtree(jobroot)
    jobroot.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []

    if lig_dir.exists():
        ligand_copy = normalize_ligand_tree(lig_dir, jobroot / "Ligands")
        warnings.extend(ligand_copy["warnings"])

    centers_exact = ws / "vina_centers.csv"
    candidate = None
    if centers_exact.exists():
        candidate = centers_exact
    else:
        candidates = sorted(ws.glob("vina_centers*.csv"), key=lambda p: p.stat().st_mtime)
        if candidates:
            candidate = candidates[-1]

    if candidate:
        temp_norm = candidate.parent / "._tmp_canonical_vina_centers.csv"
        normalize_centers_csv_to_canonical(candidate, temp_norm)
        shutil.copy2(temp_norm, jobroot / "vina_centers.csv")
        temp_norm.unlink(missing_ok=True)
    else:
        _write_canonical_centers(jobroot / "vina_centers.csv", [])

    _copy_prepared_receptors(ws, jobroot)
    _copy_raw_receptors(ws, jobroot)

    repo_root = Path(__file__).resolve().parent
    warnings.extend(stage_runtime_assets(jobroot, repo_root, package_mode))
    _chmod_scripts(jobroot)

    return jobroot, warnings


def zip_job_tree(jobroot: Path) -> Path:
    zpath = jobroot.with_suffix(".zip")
    if zpath.exists():
        zpath.unlink()
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in jobroot.rglob("*"):
            zf.write(path, path.relative_to(jobroot.parent))
    return zpath


def fetch_pdb_and_prep(pdb_code: str, dest_dir: Path, chains: str = "") -> Dict[str, Any]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    pdb_code = pdb_code.strip().lower()
    chains = (chains or "").replace(" ", "")

    attempts = [
        (f"https://files.rcsb.org/download/{pdb_code}.pdb", dest_dir / f"{pdb_code}.pdb", "pdb"),
        (f"https://files.rcsb.org/download/{pdb_code}.cif", dest_dir / f"{pdb_code}.cif", "cif"),
    ]
    raw_path: Optional[Path] = None
    source_format = ""
    errors: list[str] = []
    for url, candidate_path, fmt in attempts:
        try:
            with urllib.request.urlopen(url) as resp, open(candidate_path, "wb") as fout:
                fout.write(resp.read())
            raw_path = candidate_path
            source_format = fmt
            break
        except urllib.error.HTTPError as exc:
            errors.append(f"{candidate_path.name}: HTTP {exc.code}")
        except urllib.error.URLError as exc:
            errors.append(f"{candidate_path.name}: {exc.reason}")
    if raw_path is None:
        detail = "; ".join(errors) if errors else "no download attempts succeeded"
        raise ValueError(f"Unable to fetch receptor {pdb_code.upper()} from RCSB ({detail}).")

    if chains:
        if source_format != "pdb":
            raise ValueError(
                f"Receptor {pdb_code.upper()} is only available as mmCIF from RCSB, "
                "and chain filtering is currently supported only for PDB downloads. "
                "Fetch without chains, or upload a prepared single-chain receptor file."
            )
        wanted = {c.upper() for c in chains.split(",") if c}
        filtered = dest_dir / f"{pdb_code}_filtered.pdb"
        with open(raw_path) as fin, open(filtered, "w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    if line[21].upper() in wanted:
                        fout.write(line)
                else:
                    fout.write(line)
        filtered.replace(raw_path)

    return {"pdb_path": str(raw_path), "receptor_pdb": str(raw_path)}


_STD = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def write_cleaned_pdb(src: Path, out: Path, keep_residues: set[str]) -> None:
    with open(src) as fin, open(out, "w") as fout:
        for line in fin:
            if not line.startswith(("ATOM", "HETATM")):
                fout.write(line)
                continue
            rn = line[17:20].strip().upper()
            if rn in _STD or rn in keep_residues:
                fout.write(line)


def prep_receptor_to_pdbqt(
    single_pdb: Path,
    out_dir: Path,
    remove_het_codes: Iterable[str] | str = "all",
    remove_chains: Optional[Iterable[str]] = None,
    keep_clean_pdb: bool = True,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    script = Path(__file__).with_name("3a_PDB2PDBQTbatch.py")
    folder = single_pdb.parent
    files = [single_pdb.name]
    cfg = {
        files[0]: {
            "remove_het": ("all" if remove_het_codes == "all" else list(remove_het_codes)),
            "remove_chains": list(remove_chains or []),
        }
    }
    cfg_path = out_dir / "_perfile_config.json"
    cfg_path.write_text(json.dumps(cfg))

    args = [
        sys.executable,
        str(script),
        "--folder",
        str(folder),
        "--headless",
        "--mode",
        "per-file",
        "--per-file-config",
        str(cfg_path),
        "--files",
        files[0],
        "--output-dir",
        str(out_dir),
        "--backend",
        "auto",
        "--altloc",
        "collapse",
    ]
    if keep_clean_pdb:
        args.append("--keep-clean-pdb")

    proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out_dir / "_prep_log.txt").write_text(proc.stdout)
    return out_dir


def collect_receptor_tags(csv_path: Path) -> dict[str, str]:
    tags: dict[str, str] = {}
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return tags
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        fields = {(h or "").strip().upper(): h for h in (reader.fieldnames or [])}
        key_new = fields.get("PDB_ID")
        key_old = fields.get("RECEPTOR_PDBQT")
        tag_col = fields.get("TAG")
        for row in reader:
            if not tag_col:
                continue
            if key_new and row.get(key_new):
                tags[Path(row[key_new]).stem] = row.get(tag_col, "")
            elif key_old and row.get(key_old):
                tags[Path(row[key_old]).stem] = row.get(tag_col, "")
    return tags


def rename_centers_with_tags(jobroot: Path):
    root_csv = jobroot / "vina_centers.csv"
    _ = collect_receptor_tags(root_csv)
    # no destructive rename in this version
