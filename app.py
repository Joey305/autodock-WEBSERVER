# ==============================
# app.py — AutoDock-Vina PrepServer (finalized: correct converted detection, status flow)
# ==============================
from __future__ import annotations

import os, json, time, csv, subprocess, re, shutil, threading, uuid, zipfile, stat, signal, hmac
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional

from flask import (
    Flask, render_template, request, send_file, send_from_directory, jsonify, current_app, url_for, make_response
)
from flask_login import LoginManager, login_required
from werkzeug.utils import secure_filename

from models import db, User
from auth import auth_bp
from hpc_profiles import normalize_package_mode as normalize_hpc_package_mode, profile_for_mode
from lsf_templates import build_confgen_lsfs, build_vina_lsfs
from packager import (
    make_workspace, ensure_subdir, save_uploaded_zip,
    assemble_job_tree, zip_job_tree, fetch_pdb_and_prep, rename_centers_with_tags,
    save_uploaded_ligand_zip, save_uploaded_ligand_folder,
)
from runner_templates import build_portable_runners
from center_resolver import CenterResolutionError, resolve_center_from_file, resolve_xyz

try:
    import gemmi  # type: ignore
except Exception:
    gemmi = None

SITE_CONTACT_EMAIL = "jmschulz@med.miami.edu"
MAILTO_SUBJECT = "AutoDock-Vina PrepServer Question"
REPOSITORY_URL = "https://github.com/Joey305/autodock-WEBSERVER"
SCHURER_LAB_URL = "https://schurerlab.org"

TOOL_LINKS = [
    {
        "name": "Warhead Hunter",
        "url": "https://warheadhunter.com",
        "icon": "bi-crosshair",
        "role": "Solvent-exposed ligand atom and warhead/linker follow-up",
    },
    {
        "name": "PROTAC Builder",
        "url": "https://protacbuilder.com",
        "icon": "bi-diagram-3",
        "role": "Degrader design continuation and linker-recruiter-warhead assembly",
    },
    {
        "name": "E3 Ligandalyzer",
        "url": "https://e3ligandalyzer.com",
        "icon": "bi-bounding-box-circles",
        "role": "E3 recruiter and ligase-context exploration",
    },
    {
        "name": "V-LiSEMOD",
        "url": "https://vlisemod.com",
        "icon": "bi-grid",
        "role": "Viral ligand solvent-exposed moiety and PROTACability-style triage",
    },
    {
        "name": "Pymacs",
        "url": "https://github.com/schurerlab/Pymacs",
        "icon": "bi-github",
        "role": "Schürer Lab GitHub tool for molecular modeling and MD-related workflows",
    },
]

LAB_LINK = {
    "name": "Schürer Lab",
    "url": SCHURER_LAB_URL,
    "icon": "bi-flask",
    "role": "Laboratory context for the connected molecular design tool ecosystem",
}

CURATED_LIGAND_LIBRARY_MAP = {
    "phase4": {
        "filename": "chembl_phase4_approved_smiles.csv",
        "label": "ChEMBL Phase 4 approved drugs",
        "description": "Approved drugs in ChEMBL. Approval does not mean non-toxic in every setting or dose.",
    },
    "phase2": {
        "filename": "chembl_phase2_or_higher_smiles.csv",
        "label": "ChEMBL Phase 2 or higher compounds",
        "description": "Compounds that reached at least Phase 2 clinical testing. This does not mean they are non-toxic or safe for all uses.",
    },
}


def _env_bool(name: str, default: bool) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


def normalize_package_mode(values, default_mode: str = "portable", lsf_enabled: bool = True) -> str:
    return normalize_hpc_package_mode(values, default_mode=default_mode, lsf_enabled=lsf_enabled)


def infer_ligand_workflow(ligand_info: Dict[str, Any]) -> tuple[str, str, Optional[str]]:
    upload_mode = ligand_info.get("upload_mode")
    filename = (ligand_info.get("filename") or "").strip()
    ext = Path(filename).suffix.lower()
    filetypes = ligand_info.get("filetypes") or []
    if upload_mode in {"zip", "folder"} and set(filetypes) == {".csv"} and ligand_info.get("accepted_count") == 1:
        accepted = ligand_info.get("accepted_files") or []
        csv_name = accepted[0] if accepted else filename
        return "1", "csv", None

    if upload_mode == "single" and ext == ".csv":
        return "1", "csv", None
    if upload_mode == "single" and ext == ".sdf":
        return "3", "sdf", f"Ligands/{Path(filename).name}"
    if upload_mode == "single" and ext in {".smiles", ".smi"}:
        return "2", "smiles", None
    return "2", "sdf", None


STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "MSE",
}
STD_NT = {"A", "C", "G", "T", "U", "I", "DA", "DC", "DG", "DT", "DI", "RA", "RC", "RG", "RU"}
WATER_NAMES = {"HOH", "WAT", "H2O"}
COMMON_SUGARS = {"NAG", "BMA", "MAN", "GAL", "FUC", "NDG"}


def _is_true_het_residue(res, chain=None) -> bool:
    name = (res.name or "").strip().upper()
    if name in WATER_NAMES:
        return True
    et = getattr(res, "entity_type", None)
    if gemmi is not None and et in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Water):
        return True
    if name in STD_AA or name in STD_NT:
        return False
    if name in COMMON_SUGARS:
        return True
    return getattr(res, "het_flag", " ") != " "


def _iter_atoms_with_altloc_policy(res, policy: str = "collapse"):
    if policy == "all":
        for atom in res:
            yield atom
        return
    groups = {}
    for atom in res:
        key = atom.name.strip()
        alt = getattr(atom, "altloc", "") or ""
        occ = getattr(atom, "occ", 1.0)
        prev = groups.get(key)
        if prev is None:
            groups[key] = (atom, alt, occ)
            continue
        _, prev_alt, prev_occ = prev
        better = (occ > prev_occ) or (occ == prev_occ and (alt or " ") < (prev_alt or " "))
        if better:
            groups[key] = (atom, alt, occ)
    for atom, _, _ in groups.values():
        yield atom


def _format_pdb_atom_name(atom_name: str) -> str:
    name = atom_name.strip()
    return name[:4] if len(name) >= 4 else f"{name:>4}"


def _derive_element(atom) -> str:
    try:
        elem = atom.element.name.strip()
        if elem and elem != "X":
            return elem
    except Exception:
        pass
    letters = "".join(ch for ch in (atom.name or "") if ch.isalpha())
    if not letters:
        return ""
    return letters[:2].title() if len(letters) >= 2 and letters[1].islower() else letters[0].upper()


def _pdb_safe_char(value: Any, default: str = " ") -> str:
    text = ("" if value is None else str(value))
    if not text or text == "\x00":
        return default
    ch = text[0]
    return default if ch == "\x00" else ch


def _normalize_resname_for_pdb(resname: str, is_het: bool) -> tuple[str, Optional[str]]:
    original = (resname or "").strip().upper()
    if not original:
        return "UNK", None
    if len(original) <= 3:
        return original, None
    alias = original[:3]
    return (alias, alias) if is_het else (alias, None)


def _write_receptor_pdb_snapshot(
    src_path: Path,
    dst_path: Path,
    remove_hets: Optional[List[str]] = None,
    remove_chains: Optional[List[str]] = None,
    remove_all_hets: bool = False,
    altloc_mode: str = "collapse",
) -> Dict[str, str]:
    if gemmi is None:
        raise RuntimeError(
            "mmCIF receptor support requires the Python package `gemmi`, which is not installed in this environment."
        )
    structure = gemmi.read_structure(str(src_path))
    remove_het_set = {item.strip().upper() for item in (remove_hets or []) if item.strip()}
    remove_chain_set = {item.strip().upper() for item in (remove_chains or []) if item.strip()}
    alias_map: Dict[str, str] = {}
    serial = 1
    lines: List[str] = []
    wrote_any = False

    for model in structure:
        for chain in model:
            chain_name = (chain.name or "").strip()
            chain_id = (chain_name[:1] if chain_name else "A").upper()
            if chain_id in remove_chain_set:
                continue
            chain_wrote = False
            last_resname, last_resseq, last_icode = "UNK", 0, " "
            for res in chain:
                orig_resname = (res.name or "").strip().upper()
                is_het = _is_true_het_residue(res, chain)
                if is_het and remove_all_hets:
                    continue
                if is_het and orig_resname in remove_het_set:
                    continue
                pdb_resname, alias = _normalize_resname_for_pdb(orig_resname, is_het)
                if alias:
                    alias_map[orig_resname] = alias
                try:
                    resseq = res.seqid.num
                    icode = _pdb_safe_char(getattr(res.seqid, "icode", ""), " ")
                except Exception:
                    resseq, icode = 0, " "
                record = "HETATM" if is_het else "ATOM  "
                for atom in _iter_atoms_with_altloc_policy(res, policy=altloc_mode):
                    x, y, z = atom.pos.x, atom.pos.y, atom.pos.z
                    occ = getattr(atom, "occ", 1.00)
                    bfac = getattr(atom, "b_iso", 0.00)
                    altloc = _pdb_safe_char(getattr(atom, "altloc", ""), " ")
                    atom_name = _format_pdb_atom_name(atom.name)
                    element = _derive_element(atom)
                    lines.append(
                        f"{record:<6}{serial:>5} "
                        f"{atom_name:<4}{altloc:1}"
                        f"{pdb_resname:>3} "
                        f"{chain_id:1}{resseq:>4}{icode:1}"
                        f"   "
                        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                        f"{occ:>6.2f}{bfac:>6.2f}"
                        f"          {element:>2}  "
                    )
                    serial += 1
                    wrote_any = True
                    chain_wrote = True
                    last_resname, last_resseq, last_icode = pdb_resname, resseq, icode
            if chain_wrote:
                lines.append(f"{'TER':<6}{serial:>5}      {last_resname:>3} {chain_id:1}{last_resseq:>4}{last_icode:1}")
                serial += 1
    lines.append("END")
    if not wrote_any:
        raise RuntimeError(f"No atoms written from {src_path.name} after filtering.")
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    dst_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return alias_map


def _curated_ligand_source(library_key: str) -> Optional[Dict[str, Any]]:
    entry = CURATED_LIGAND_LIBRARY_MAP.get((library_key or "").strip().lower())
    if not entry:
        return None
    path = Path(current_app.static_folder or "static") / "data" / entry["filename"]
    if not path.exists():
        return None
    headers: List[str] = []
    try:
        with path.open("r", newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            headers = next(reader, []) or []
    except Exception:
        headers = []
    smiles_col = "canonical_smiles" if "canonical_smiles" in headers else ""
    id_candidates = ["pref_name", "molecule_chembl_id", "compound_id", "id", "name"]
    id_col = next((col for col in id_candidates if col in headers), "")
    return {
        "key": library_key.strip().lower(),
        "path": path,
        "headers": headers,
        "suggested_smiles_col": smiles_col,
        "suggested_id_col": id_col,
        **entry,
    }


def _csv_upload_metadata(lig_dir: Path, result: Dict[str, Any]) -> Dict[str, Any]:
    if not result.get("is_csv"):
        return {}
    accepted = result.get("accepted_files") or []
    if len(accepted) != 1:
        return {}
    csv_path = lig_dir / accepted[0]
    if not csv_path.exists():
        return {}
    try:
        with csv_path.open("r", newline="", encoding="utf-8") as fh:
            headers = next(csv.reader(fh), []) or []
    except Exception:
        headers = []
    smiles_col = "canonical_smiles" if "canonical_smiles" in headers else ""
    id_candidates = ["pref_name", "molecule_chembl_id", "compound_id", "id", "name"]
    id_col = next((col for col in id_candidates if col in headers), "")
    return {
        "headers": headers,
        "suggested_smiles_col": smiles_col,
        "suggested_id_col": id_col,
    }

# ---------- Config ----------
class Config:
    SECRET_KEY = os.getenv("PORTAL_SECRET", "dev-secret-change-me")
    SQLALCHEMY_DATABASE_URI = os.getenv("PORTAL_DB", "sqlite:///portal.db")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    TMP_ROOT = os.getenv("PORTAL_TMP", "/tmp/autodock_prep")
    MAX_RECEPTORS = int(os.getenv("PORTAL_MAX_RECEPTORS", "5"))
    DEFAULT_ENV_LINE = os.getenv("PORTAL_ENV_LINE", "")
    PUBLIC_EMAIL = os.getenv("PORTAL_PUBLIC_EMAIL", "public@autodock.local")
    PUBLIC_MODE = _env_bool("PUBLIC_MODE", True)
    ENABLE_AUTH = _env_bool("ENABLE_AUTH", False)
    ENABLE_LSF_PACKAGE = _env_bool("ENABLE_LSF_PACKAGE", True)
    DEFAULT_PACKAGE_MODE = os.getenv("DEFAULT_PACKAGE_MODE", "portable").strip().lower() or "portable"
    VTS_RUNNER_ROOT = os.getenv("VTS_RUNNER_ROOT", "")
    VTS_RUNNER_TOKEN = os.getenv("VTS_RUNNER_TOKEN", "")
    VTS_RUNNER_CONDA_ENV = os.getenv("VTS_RUNNER_CONDA_ENV", "docking")
    VTS_RUNNER_ENV_LINE = os.getenv("VTS_RUNNER_ENV_LINE", "")
    VTS_RUNNER_TIMEOUT_SECONDS = int(os.getenv("VTS_RUNNER_TIMEOUT_SECONDS", "1800"))
    VTS_RUNNER_MAX_ZIP_MB = int(os.getenv("VTS_RUNNER_MAX_ZIP_MB", "1024"))


class PublicUser:
    """Minimal compatibility user used when the app runs without login."""

    is_authenticated = True
    is_active = True
    is_anonymous = False
    is_public = True
    role = "public"

    def __init__(self):
        self.email = os.getenv("PORTAL_PUBLIC_EMAIL", "public@autodock.local")
        self.id = None

    def get_id(self):
        return None

login_manager = LoginManager()
login_manager.anonymous_user = PublicUser

APP_ROOT = Path(__file__).resolve().parent

def create_app() -> Flask:
    app = Flask(__name__)
    app.config.from_object(Config)
    Path(app.config["TMP_ROOT"]).mkdir(parents=True, exist_ok=True)
    if not app.config["VTS_RUNNER_ROOT"]:
        app.config["VTS_RUNNER_ROOT"] = str(Path(app.config["TMP_ROOT"]) / "vts_runs")
    Path(app.config["VTS_RUNNER_ROOT"]).mkdir(parents=True, exist_ok=True)

    db.init_app(app)
    with app.app_context():
        db.create_all()

    login_manager.init_app(app)

    @login_manager.user_loader
    def _load(uid):
        return db.session.get(User, int(uid))

    app.register_blueprint(auth_bp)

    def _public_email() -> str:
        return current_app.config["PUBLIC_EMAIL"]

    def _public_name() -> str:
        return _public_email().split("@")[0] or "public"

    @app.context_processor
    def inject_site_links():
        mailto = f"mailto:{SITE_CONTACT_EMAIL}?subject={MAILTO_SUBJECT.replace(' ', '%20')}"
        nav_links = [
            {"name": "Home", "endpoint": "home"},
            {"name": "Build", "endpoint": "build"},
            {"name": "Results", "endpoint": "results"},
            {"name": "Workflow", "endpoint": "workflow"},
            {"name": "About", "endpoint": "about"},
            {"name": "Docs", "endpoint": "documentation"},
            {"name": "GitHub", "url": REPOSITORY_URL, "external": True},
        ]
        return {
            "contact_email": SITE_CONTACT_EMAIL,
            "contact_mailto": mailto,
            "repository_url": REPOSITORY_URL,
            "tool_links": TOOL_LINKS,
            "lab_link": LAB_LINK,
            "ecosystem_links": [*TOOL_LINKS, LAB_LINK],
            "nav_links": nav_links,
            "public_mode": current_app.config.get("PUBLIC_MODE", True),
        }

    def _resolve_casefold_path(base_dir: Path, relative_path: Path) -> Optional[Path]:
        current = base_dir
        if not current.exists():
            return None
        for part in relative_path.parts:
            exact = current / part
            if exact.exists():
                current = exact
                continue
            try:
                children = {child.name.lower(): child for child in current.iterdir()}
            except Exception:
                return None
            match = children.get(part.lower())
            if match is None:
                return None
            current = match
        return current.resolve()

    def _resolve_workspace_file(ws: Path, rel: str) -> Optional[Path]:
        rel_path = Path(rel)
        if rel_path.is_absolute() or ".." in rel_path.parts:
            return None

        first = rel_path.parts[0] if rel_path.parts else ""
        remainder = Path(*rel_path.parts[1:]) if len(rel_path.parts) > 1 else Path()
        candidate_dirs: List[Path] = []

        if first == "Receptors":
            for name in ("Receptors", "Receptors_PDB", "Receptors_PDBQT", "Receptors_PDBQT_Converted"):
                candidate_dirs.append(ws / name)
        elif first:
            candidate_dirs.append(ws / first)
        else:
            candidate_dirs.append(ws)

        for base_dir in candidate_dirs:
            resolved = _resolve_casefold_path(base_dir, remainder)
            if resolved and str(resolved).startswith(str(ws.resolve())) and resolved.exists():
                return resolved

        resolved = _resolve_casefold_path(ws, rel_path)
        if resolved and str(resolved).startswith(str(ws.resolve())) and resolved.exists():
            return resolved
        return None

    def _ensure_receptor_pdb_snapshot(ws: Path, receptor_rel: str, *, altloc_mode: str = "collapse") -> Path:
        rel_path = Path(receptor_rel)
        stem = rel_path.stem
        src = _resolve_workspace_file(ws, receptor_rel)
        if src is None and rel_path.suffix.lower() == ".pdb":
            for ext in (".cif", ".mmcif"):
                candidate = _resolve_workspace_file(ws, str(Path("Receptors") / f"{stem}{ext}"))
                if candidate is not None:
                    src = candidate
                    break
        if src is None:
            raise FileNotFoundError(f"Receptor source not found for {receptor_rel}")
        if src.suffix.lower() not in {".cif", ".mmcif"}:
            return src
        dst = ensure_subdir(ws, "Receptors_PDB") / f"{stem}.pdb"
        if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
            return dst
        alias_map = _write_receptor_pdb_snapshot(src, dst, altloc_mode=altloc_mode)
        if alias_map:
            dst.with_suffix(".aliases.json").write_text(json.dumps(alias_map, indent=2, sort_keys=True), encoding="utf-8")
        return dst

    def _path_within(base: Path, candidate: Path) -> bool:
        try:
            return str(candidate.resolve()).startswith(str(base.resolve()))
        except Exception:
            return False

    vts_jobs_lock = threading.Lock()

    def _vts_runner_root() -> Path:
        root = Path(current_app.config["VTS_RUNNER_ROOT"])
        root.mkdir(parents=True, exist_ok=True)
        return root

    def _vts_token_ok() -> bool:
        expected = current_app.config.get("VTS_RUNNER_TOKEN", "")
        if not expected:
            return True
        supplied = (
            request.args.get("token")
            or request.form.get("token")
            or request.headers.get("X-VTS-Runner-Token")
            or ""
        )
        return hmac.compare_digest(str(supplied), str(expected))

    def _vts_status_path(run_dir: Path) -> Path:
        return run_dir / "status.json"

    def _write_vts_status(run_dir: Path, status: Dict[str, Any]) -> None:
        status_path = _vts_status_path(run_dir)
        tmp_path = status_path.with_suffix(".json.tmp")
        tmp_path.write_text(json.dumps(status, indent=2, sort_keys=True), encoding="utf-8")
        tmp_path.replace(status_path)

    def _read_vts_status(run_dir: Path) -> Dict[str, Any]:
        status_path = _vts_status_path(run_dir)
        if not status_path.exists():
            return {}
        try:
            return json.loads(status_path.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _update_vts_status(run_dir: Path, **updates: Any) -> Dict[str, Any]:
        with vts_jobs_lock:
            status = _read_vts_status(run_dir)
            status.update(updates)
            status["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
            _write_vts_status(run_dir, status)
            return status

    def _tail_text(path: Path, max_chars: int = 16000) -> str:
        if not path.exists():
            return ""
        try:
            data = path.read_bytes()
        except Exception:
            return ""
        return data[-max_chars:].decode("utf-8", errors="replace")

    def _safe_extract_vts_zip(zip_path: Path, extract_dir: Path, max_uncompressed_bytes: int) -> None:
        if not zipfile.is_zipfile(zip_path):
            raise ValueError("Upload must be a valid .zip file.")

        base = extract_dir.resolve()
        total = 0
        with zipfile.ZipFile(zip_path) as zf:
            members = zf.infolist()
            if not members:
                raise ValueError("The zip archive is empty.")
            for info in members:
                name = info.filename.replace("\\", "/")
                if not name or name.startswith("__MACOSX/") or "/__MACOSX/" in name:
                    continue
                if name.startswith("/") or re.match(r"^[A-Za-z]:/", name):
                    raise ValueError(f"Unsafe absolute path in zip: {info.filename}")
                parts = Path(name).parts
                if any(part in {"", ".", ".."} for part in parts):
                    raise ValueError(f"Unsafe relative path in zip: {info.filename}")
                mode = (info.external_attr >> 16) & 0o170000
                if mode == stat.S_IFLNK:
                    raise ValueError(f"Symbolic links are not allowed in zip uploads: {info.filename}")
                total += info.file_size
                if total > max_uncompressed_bytes:
                    raise ValueError("The expanded zip is too large for the VTS runner.")
                target = (extract_dir / name).resolve()
                try:
                    target.relative_to(base)
                except ValueError as exc:
                    raise ValueError(f"Unsafe extraction path in zip: {info.filename}") from exc

            for info in members:
                name = info.filename.replace("\\", "/")
                if name and not name.startswith("__MACOSX/") and "/__MACOSX/" not in name:
                    zf.extract(info, extract_dir)

    def _detect_vts_package_root(extract_dir: Path) -> Path:
        visible_children = [
            p for p in extract_dir.iterdir()
            if p.name not in {".DS_Store", "__MACOSX"} and not p.name.startswith("._")
        ]
        dirs = [p for p in visible_children if p.is_dir()]
        files = [p for p in visible_children if p.is_file()]
        if len(dirs) == 1 and not files:
            return dirs[0]
        return extract_dir

    def _find_vts_runner_script(package_root: Path) -> Path:
        candidates = [
            p for p in package_root.rglob("*.sh")
            if p.is_file()
            and not any(part == "__MACOSX" for part in p.parts)
            and p.name.startswith("run")
            and ("linux" in p.name.lower() or "macos" in p.name.lower())
        ]
        if not candidates:
            raise ValueError("No Linux/macOS runner script was found. Expected a file like run_*_macos_linux.sh.")
        candidates.sort(key=lambda p: (len(p.relative_to(package_root).parts), p.name.lower()))
        return candidates[0]

    def _vts_command() -> List[str]:
        bootstrap = r'''
set -e
if [ -n "${VTS_RUNNER_ENV_LINE:-}" ]; then
  eval "$VTS_RUNNER_ENV_LINE"
elif [ -n "${VTS_RUNNER_CONDA_ENV:-}" ]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$VTS_RUNNER_CONDA_ENV"
  elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    . "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate "$VTS_RUNNER_CONDA_ENV"
  elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    . "$HOME/anaconda3/etc/profile.d/conda.sh"
    conda activate "$VTS_RUNNER_CONDA_ENV"
  else
    echo "[BOB-VTS] Conda environment requested but conda was not found; continuing with current Python."
  fi
fi
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python" ]; then
  mkdir -p "$PWD/.bob_vts_bin"
  ln -sf "$CONDA_PREFIX/bin/python" "$PWD/.bob_vts_bin/python3"
  export PATH="$PWD/.bob_vts_bin:$CONDA_PREFIX/bin:$PATH"
fi
echo "[BOB-VTS] CWD=$(pwd)"
echo "[BOB-VTS] Runner=$VTS_SCRIPT"
echo "[BOB-VTS] Python=$(command -v python3 || command -v python || true)"
python3 -V 2>/dev/null || python -V 2>/dev/null || true
exec bash "$VTS_SCRIPT"
'''
        return ["bash", "-lc", bootstrap]

    def _run_vts_job(run_dir: Path, package_root: Path, script_path: Path) -> None:
        log_path = run_dir / "bob_vts_terminal_output.txt"
        started = time.time()
        _update_vts_status(
            run_dir,
            status="running",
            started_at=time.strftime("%Y-%m-%d %H:%M:%S"),
            script=str(script_path.relative_to(package_root)),
        )

        env = os.environ.copy()
        env.update(
            {
                "PYTHONUNBUFFERED": "1",
                "MPLBACKEND": env.get("MPLBACKEND", "Agg"),
                "VTS_NUM_THREADS": env.get("VTS_NUM_THREADS", "1"),
                "OMP_NUM_THREADS": env.get("OMP_NUM_THREADS", "1"),
                "MKL_NUM_THREADS": env.get("MKL_NUM_THREADS", "1"),
                "OPENBLAS_NUM_THREADS": env.get("OPENBLAS_NUM_THREADS", "1"),
                "VTS_SCRIPT": str(script_path),
                "VTS_RUNNER_CONDA_ENV": str(app.config.get("VTS_RUNNER_CONDA_ENV", "")),
                "VTS_RUNNER_ENV_LINE": str(app.config.get("VTS_RUNNER_ENV_LINE", "")),
            }
        )

        return_code: Optional[int] = None
        timed_out = False
        error = ""
        try:
            script_path.chmod(script_path.stat().st_mode | 0o755)
            with log_path.open("w", encoding="utf-8", errors="replace") as log_file:
                proc = subprocess.Popen(
                    _vts_command(),
                    cwd=package_root,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    env=env,
                    start_new_session=True,
                    text=True,
                )
                try:
                    return_code = proc.wait(timeout=int(app.config["VTS_RUNNER_TIMEOUT_SECONDS"]))
                except subprocess.TimeoutExpired:
                    timed_out = True
                    os.killpg(proc.pid, signal.SIGTERM)
                    try:
                        return_code = proc.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(proc.pid, signal.SIGKILL)
                        return_code = proc.wait(timeout=10)
        except Exception as exc:
            error = str(exc)
            return_code = 1 if return_code is None else return_code

        elapsed = round(time.time() - started, 2)
        status_name = "succeeded" if return_code == 0 and not timed_out and not error else "failed"
        metadata = {
            "status": status_name,
            "return_code": return_code,
            "timed_out": timed_out,
            "elapsed_seconds": elapsed,
            "script": str(script_path.relative_to(package_root)),
            "error": error,
            "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        try:
            (package_root / "bob_vts_run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
            shutil.copyfile(log_path, package_root / "bob_vts_terminal_output.txt")
            zip_stem = secure_filename(f"{package_root.name}_BOB_VTS_results") or "BOB_VTS_results"
            result_zip = run_dir / f"{zip_stem}.zip"
            if result_zip.exists():
                result_zip.unlink()
            shutil.make_archive(str(result_zip.with_suffix("")), "zip", root_dir=package_root.parent, base_dir=package_root.name)
            metadata["result_zip"] = str(result_zip)
        except Exception as exc:
            status_name = "failed"
            metadata["status"] = status_name
            metadata["error"] = f"{error}; result zip failed: {exc}".strip("; ")

        _update_vts_status(run_dir, **metadata)

    def _public_vts_status(run_id: str, run_dir: Path) -> Dict[str, Any]:
        status = _read_vts_status(run_dir)
        if not status:
            return {}
        result_zip_value = status.get("result_zip", "")
        result_zip = Path(result_zip_value) if result_zip_value else None
        payload = {
            "run_id": run_id,
            "status": status.get("status", "unknown"),
            "original_filename": status.get("original_filename", ""),
            "package_name": status.get("package_name", ""),
            "script": status.get("script", ""),
            "return_code": status.get("return_code"),
            "timed_out": bool(status.get("timed_out", False)),
            "elapsed_seconds": status.get("elapsed_seconds"),
            "started_at": status.get("started_at", ""),
            "completed_at": status.get("completed_at", ""),
            "updated_at": status.get("updated_at", ""),
            "error": status.get("error", ""),
            "log_tail": _tail_text(run_dir / "bob_vts_terminal_output.txt"),
        }
        if result_zip and result_zip.exists():
            download_args = {"run_id": run_id}
            if current_app.config.get("VTS_RUNNER_TOKEN"):
                token = request.args.get("token") or request.headers.get("X-VTS-Runner-Token") or ""
                download_args["token"] = token
            payload["download_url"] = url_for("bob_vts_download", **download_args)
        return payload

    # ---------- PAGES ----------
    @app.get("/robots.txt")
    def robots_txt():
        return send_from_directory(app.static_folder, "robots.txt", mimetype="text/plain")

    @app.get("/sitemap.xml")
    def sitemap_xml():
        return send_from_directory(app.static_folder, "sitemap.xml", mimetype="application/xml")

    @app.get("/llms.txt")
    def llms_txt():
        return send_from_directory(app.static_folder, "llms.txt", mimetype="text/markdown")

    @app.get("/llms-full.txt")
    def llms_full_txt():
        return send_from_directory(app.static_folder, "llms-full.txt", mimetype="text/markdown")

    @app.get("/")
    def home():
        return render_template("home.html")

    @app.get("/build")
    def build():
        return render_template(
            "build.html",
            max_receptors=app.config["MAX_RECEPTORS"],
            enable_lsf_package=app.config["ENABLE_LSF_PACKAGE"],
            default_package_mode=normalize_package_mode(
                {},
                default_mode=app.config["DEFAULT_PACKAGE_MODE"],
                lsf_enabled=app.config["ENABLE_LSF_PACKAGE"],
            ),
        )

    @app.get("/about")
    def about():
        return render_template("about.html")

    @app.get("/workflow")
    def workflow():
        return render_template("workflow.html")

    @app.get("/results")
    def results():
        return render_template("results.html", example_project_url=url_for("docking_visualization_example"))

    @app.get("/modules")
    def modules():
        return render_template("modules.html")

    @app.get("/documentation")
    def documentation():
        return render_template("documentation.html")

    @app.get("/README.md")
    def readme_markdown():
        return send_file(APP_ROOT / "README.md", mimetype="text/markdown")

    @app.get("/contact")
    def contact():
        return render_template(
            "contact.html",
            repository_url=REPOSITORY_URL
        )
    
    @app.get("/docking-preparation-problems")
    def docking_problems():
        return render_template("docking_problems.html")

    @app.get("/open-source")
    def open_source():
        return render_template("open_source.html")

    @app.get("/api")
    def api_docs():
        return render_template("api.html")

    @app.get("/troubleshooting")
    def troubleshooting():
        return render_template("troubleshooting.html")

    @app.get("/BOB-VTS.html")
    @app.get("/bob-vts.html")
    @login_required
    def bob_vts_page():
        if not _vts_token_ok():
            return ("not found", 404)
        response = make_response(
            render_template(
                "BOB-VTS.html",
                token_required=bool(current_app.config.get("VTS_RUNNER_TOKEN")),
                supplied_token=request.args.get("token", ""),
                max_zip_mb=current_app.config["VTS_RUNNER_MAX_ZIP_MB"],
                timeout_minutes=round(int(current_app.config["VTS_RUNNER_TIMEOUT_SECONDS"]) / 60),
            )
        )
        response.headers["X-Robots-Tag"] = "noindex, nofollow"
        return response

    @app.post("/api/bob-vts/runs")
    @login_required
    def bob_vts_start():
        if not _vts_token_ok():
            return jsonify({"error": "not found"}), 404

        max_zip_bytes = int(current_app.config["VTS_RUNNER_MAX_ZIP_MB"]) * 1024 * 1024
        if request.content_length and request.content_length > max_zip_bytes:
            return jsonify({"error": f"Upload is larger than {current_app.config['VTS_RUNNER_MAX_ZIP_MB']} MB."}), 413

        upload = request.files.get("package")
        if not upload or not upload.filename:
            return jsonify({"error": "Drop or choose a VTS .zip package first."}), 400
        if Path(upload.filename).suffix.lower() != ".zip":
            return jsonify({"error": "Only .zip packages are accepted."}), 400

        run_id = uuid.uuid4().hex[:16]
        run_dir = _vts_runner_root() / run_id
        upload_dir = run_dir / "upload"
        extract_dir = run_dir / "extracted"
        upload_dir.mkdir(parents=True)
        extract_dir.mkdir(parents=True)

        original_name = secure_filename(upload.filename) or "vts_package.zip"
        upload_path = upload_dir / original_name
        upload.save(upload_path)

        try:
            _safe_extract_vts_zip(upload_path, extract_dir, max_zip_bytes)
            package_root = _detect_vts_package_root(extract_dir)
            script_path = _find_vts_runner_script(package_root)
        except ValueError as exc:
            shutil.rmtree(run_dir, ignore_errors=True)
            return jsonify({"error": str(exc)}), 400

        _write_vts_status(
            run_dir,
            {
                "run_id": run_id,
                "status": "queued",
                "original_filename": original_name,
                "package_name": package_root.name,
                "package_root": str(package_root),
                "script": str(script_path.relative_to(package_root)),
                "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            },
        )

        thread = threading.Thread(target=_run_vts_job, args=(run_dir, package_root, script_path), daemon=True)
        thread.start()

        payload = _public_vts_status(run_id, run_dir)
        status_args = {"run_id": run_id}
        if current_app.config.get("VTS_RUNNER_TOKEN"):
            status_args["token"] = request.form.get("token", "")
        payload["status_url"] = url_for("bob_vts_status", **status_args)
        return jsonify(payload), 202

    @app.get("/api/bob-vts/runs/<run_id>")
    @login_required
    def bob_vts_status(run_id: str):
        if not _vts_token_ok():
            return jsonify({"error": "not found"}), 404
        if not re.fullmatch(r"[a-f0-9]{16}", run_id):
            return jsonify({"error": "Invalid run id."}), 400
        run_dir = _vts_runner_root() / run_id
        payload = _public_vts_status(run_id, run_dir)
        if not payload:
            return jsonify({"error": "Run not found."}), 404
        return jsonify(payload)

    @app.get("/api/bob-vts/runs/<run_id>/download")
    @login_required
    def bob_vts_download(run_id: str):
        if not _vts_token_ok():
            return ("not found", 404)
        if not re.fullmatch(r"[a-f0-9]{16}", run_id):
            return ("invalid run id", 400)
        run_dir = _vts_runner_root() / run_id
        status = _read_vts_status(run_dir)
        result_zip = Path(status.get("result_zip", "")) if status.get("result_zip") else None
        if not result_zip or not result_zip.exists():
            return ("result zip not found", 404)
        try:
            result_zip.resolve().relative_to(run_dir.resolve())
        except ValueError:
            return ("result zip not found", 404)
        return send_file(result_zip, as_attachment=True, download_name=result_zip.name)

    @app.get("/viz/project")
    @login_required
    def docking_visualization_project():
        jobname = request.args.get("jobname", "").strip()
        rel = request.args.get("rel", "").strip()
        if not jobname or not rel:
            return ("missing jobname or rel", 400)
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 404)
        manifest_path = _resolve_workspace_file(ws, rel)
        if manifest_path is None or manifest_path.name != "manifest.json":
            return ("manifest not found", 404)
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            return ("manifest is invalid", 400)

        entries = []
        raw_entries = manifest.get("entries", [])
        for idx, entry in enumerate(raw_entries):
            viewer_rel = (entry.get("viewer_file") or "").strip()
            viewer_path = _resolve_workspace_file(ws, str(Path(rel).parent / viewer_rel))
            if viewer_path is None or not _path_within(ws, viewer_path):
                continue
            entries.append(
                {
                    **entry,
                    "entry_index": idx,
                    "viewer_url": url_for(
                        "api_wsinline",
                        jobname=jobname,
                        rel=str((Path(rel).parent / viewer_rel).as_posix()),
                    ),
                    "project_view_url": url_for(
                        "docking_visualization_project",
                        jobname=jobname,
                        rel=rel,
                        entry=idx,
                        view="full",
                    ),
                }
            )
        selected_index = request.args.get("entry", default=0, type=int)
        if not entries:
            selected_index = 0
        else:
            selected_index = max(0, min(selected_index, len(entries) - 1))
        standalone_mode = request.args.get("view", "").strip().lower() == "full"
        selected_entry = entries[selected_index] if entries else None
        return render_template(
            "docking_visualization_project.html",
            manifest=manifest,
            entries=entries,
            selected_entry=selected_entry,
            selected_index=selected_index,
            standalone_mode=standalone_mode,
        )

    @app.get("/viz/example")
    def docking_visualization_example():
        example_root = Path(current_app.static_folder or "static") / "data" / "examples" / "dr7"
        manifest_path = example_root / "Docking_HTML_Viz_Project" / "manifest.json"
        if not manifest_path.exists():
            return ("example manifest not found", 404)
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            return ("example manifest is invalid", 400)

        base_dir = manifest_path.parent
        entries = []
        raw_entries = manifest.get("entries", [])
        for idx, entry in enumerate(raw_entries):
            viewer_rel = (entry.get("viewer_file") or "").strip()
            if not viewer_rel:
                continue
            viewer_path = (base_dir / viewer_rel).resolve()
            try:
                viewer_path.relative_to(base_dir.resolve())
            except ValueError:
                continue
            if not viewer_path.exists():
                continue
            viewer_static_rel = viewer_path.relative_to(Path(current_app.static_folder or "static")).as_posix()
            entries.append(
                {
                    **entry,
                    "entry_index": idx,
                    "viewer_url": url_for("static", filename=viewer_static_rel),
                    "project_view_url": url_for(
                        "docking_visualization_example",
                        entry=idx,
                        view="full",
                    ),
                }
            )
        selected_index = request.args.get("entry", default=0, type=int)
        if not entries:
            selected_index = 0
        else:
            selected_index = max(0, min(selected_index, len(entries) - 1))
        standalone_mode = request.args.get("view", "").strip().lower() == "full"
        selected_entry = entries[selected_index] if entries else None
        return render_template(
            "docking_visualization_project.html",
            manifest=manifest,
            entries=entries,
            selected_entry=selected_entry,
            selected_index=selected_index,
            standalone_mode=standalone_mode,
        )
    

    # ---------- STATE HELPERS ----------
    def _ws(jobname: str) -> Path:
        return Path(current_app.config["TMP_ROOT"]) / jobname

    def _load_state(ws: Path) -> Dict[str, Any]:
        s = ws / "_state.json"
        if not s.exists():
            return {"receptors": [], "centers_csv": "vina_centers.csv",
                    "ligands_uploaded": False, "ligand_info": {}, "prep_job": None}
        return json.loads(s.read_text())

    def _save_state(ws: Path, obj: Dict[str, Any]):
        (ws / "_state.json").write_text(json.dumps(obj, indent=2))

    # ---------- CENTERS CSV HELPERS (now canonical: PDB_ID,X,Y,Z,SIZE) ----------
    def _centers_csv_path(ws: Path, st: Dict[str, Any]) -> Path:
        """
        Ensure a centers CSV exists. Canonical headers are:
            PDB_ID,X,Y,Z,SIZE
        """
        name = st.get("centers_csv") or "vina_centers.csv"
        p = ws / name
        if not p.exists():
            with p.open("w", newline="") as f:
                csv.writer(f).writerow(["PDB_ID", "X", "Y", "Z", "SIZE"])
        return p

    def _read_centers(ws: Path, st: Dict[str, Any]) -> Dict[str, Tuple[float, float, float, float]]:
        """
        Read centers as a dict: {PDB_ID: (X,Y,Z,SIZE)}.
        Accepts either the new schema (PDB_ID,X,Y,Z,SIZE) or the legacy schema
        (receptor_pdbqt, center_x, center_y, center_z, size).
        """
        p = _centers_csv_path(ws, st)
        out: Dict[str, Tuple[float, float, float, float]] = {}
        with p.open() as f:
            r = csv.DictReader(f)
            # Normalize header names
            fields = { (h or "").strip().upper(): h for h in (r.fieldnames or []) }

            has_new = all(k in fields for k in ("PDB_ID","X","Y","Z"))
            has_old = all(k in fields for k in ("RECEPTOR_PDBQT","CENTER_X","CENTER_Y","CENTER_Z"))

            for row in r:
                try:
                    if has_new:
                        key = row.get(fields["PDB_ID"], "").strip()
                        x = float(row.get(fields["X"], "0"))
                        y = float(row.get(fields["Y"], "0"))
                        z = float(row.get(fields["Z"], "0"))
                        s = float(row.get(fields.get("SIZE","SIZE"), row.get("SIZE", "20") ))
                    elif has_old:
                        key = row.get(fields["RECEPTOR_PDBQT"], "").strip()
                        x = float(row.get(fields["CENTER_X"], "0"))
                        y = float(row.get(fields["CENTER_Y"], "0"))
                        z = float(row.get(fields["CENTER_Z"], "0"))
                        s = float(row.get(fields.get("SIZE","SIZE"), row.get("size", "20")))
                    else:
                        # Unknown schema; skip row
                        continue
                    if key:
                        out[key] = (x, y, z, s)
                except Exception:
                    continue
        return out

    def _write_centers(ws: Path, st: Dict[str, Any], mapping: Dict[str, Tuple[float, float, float, float]]):
        """
        Write centers using the canonical schema: PDB_ID,X,Y,Z,SIZE
        """
        p = _centers_csv_path(ws, st)
        with p.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["PDB_ID", "X", "Y", "Z", "SIZE"])
            for k, (x, y, z, s) in mapping.items():
                w.writerow([k, x, y, z, s])

    def _upsert_center_row(ws: Path, st: Dict[str, Any], receptor_pdbqt: str,
                           center: Tuple[float, float, float], size: float):
        """
        Upsert by PDB_ID (we use the PDBQT filename key). Always writes canonical schema.
        """
        mapping = _read_centers(ws, st)
        mapping[receptor_pdbqt] = (float(center[0]), float(center[1]), float(center[2]), float(size))
        _write_centers(ws, st, mapping)

    def _clean_pdb(
        in_path: Path,
        out_path: Path,
        remove_hets: list[str],
        remove_chains: list[str],
        remove_all_hets: bool = False,
        altloc_mode: str = "collapse",
    ):
        """
        Clean PDB file:
        - Remove selected HETATMs, chains, or all HETATMs.
        - Handle altLoc (collapse=first conformer only).
        - Renumber atoms sequentially (required for valid CONECT).
        - Drop or fix CONECT records to only reference surviving atoms.
        """
        serial_map = {}
        seen_altlocs = {}
        new_lines = []
        new_serial = 1

        with open(in_path) as fin:
            for line in fin:
                rec = line[:6].strip().upper()

                if rec in ("ATOM", "HETATM"):
                    resname = line[17:20].strip().upper()
                    chain   = line[21].strip().upper()
                    altloc  = line[16].strip()
                    old_serial = int(line[6:11])

                    if remove_all_hets and rec == "HETATM":
                        continue
                    if resname in remove_hets:
                        continue
                    if chain in remove_chains:
                        continue

                    if altloc and altloc_mode == "collapse":
                        key = (line[12:16].strip(), resname, chain)
                        if key in seen_altlocs:
                            continue
                        seen_altlocs[key] = True

                    serial_map[old_serial] = new_serial
                    line = line[:6] + f"{new_serial:5d}" + line[11:]
                    new_lines.append(line)
                    new_serial += 1

                elif rec == "CONECT":
                    continue
                else:
                    new_lines.append(line)

        # rebuild CONECT
        with open(in_path) as fin:
            for line in fin:
                if not line.startswith("CONECT"):
                    continue
                refs = [line[i:i+5] for i in range(6, len(line), 5) if line[i:i+5].strip()]
                refs = [int(r) for r in refs if r.strip()]
                new_refs = [serial_map[r] for r in refs if r in serial_map]
                if not new_refs:
                    continue
                base_old = int(line[6:11])
                if base_old not in serial_map:
                    continue
                base_new = serial_map[base_old]
                new_line = f"CONECT{base_new:5d}" + "".join(f"{r:5d}" for r in new_refs) + "\n"
                new_lines.append(new_line)

        with open(out_path, "w") as fout:
            fout.writelines(new_lines)





    
    # ---------- PDB/PDBQT first-model stripper ----------
    def _write_first_model_only(src_path: Path, dst_path: Path) -> None:
        """
        Copy src_path to dst_path but keep only the first MODEL...ENDMDL block.
        If no MODEL lines exist, copy the whole file. Works for .pdb and .pdbqt.
        """
        in_model = False
        got_model = False
        buf = []
        with open(src_path, "r") as fin:
            for line in fin:
                # PDB / PDBQT are column-based but MODEL/ENDMDL are at column 1
                rec = line[:6].strip().upper()
                if rec == "MODEL" and not got_model:
                    in_model = True
                    got_model = True
                    buf.append(line)
                    continue
                if rec == "ENDMDL" and in_model:
                    buf.append(line)
                    in_model = False
                    break  # stop after first model
                if in_model:
                    buf.append(line)

        if got_model:
            # We captured a MODEL block
            dst_path.write_text("".join(buf))
        else:
            # No MODEL section; copy the whole file
            dst_path.write_text(src_path.read_text())

    
    
    
    # ---------- WORKSPACE ----------
    @app.post("/api/workspace")
    @login_required
    def api_workspace():
        stamp = time.strftime("%m-%d-%Y-%H-%M-%S")
        uname = _public_name()
        jobname = f"{stamp}-{uname}"
        ws = make_workspace(_ws(jobname))
        _save_state(ws, {"receptors": [], "centers_csv": "vina_centers.csv",
                         "ligands_uploaded": False, "ligand_info": {}, "prep_job": None})
        return jsonify({"jobname": jobname, "workspace": str(ws)})

    # ---------- RECEPTOR INGEST ----------
    @app.post("/api/receptors/upload")
    @login_required
    def api_receptors_upload():
        jobname = request.form.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        mode = request.form.get("mode")
        f = request.files.get("file")
        if not f:
            return ("no file", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        added: List[str] = []
        if mode == "single":
            out = rec_dir / Path(f.filename).name
            out.parent.mkdir(parents=True, exist_ok=True)
            f.save(out)
            if out.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                added.append(str(Path("Receptors") / out.name))
        elif mode == "zip":
            save_uploaded_zip(f, rec_dir)
            for p in sorted(rec_dir.rglob("*")):
                if p.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                    added.append(str(Path("Receptors") / p.relative_to(rec_dir)))
        else:
            return ("bad mode", 400)

        st = _load_state(ws)
        have = {r["rel"] for r in st["receptors"]}
        for rel in added:
            if rel not in have and len(st["receptors"]) < current_app.config["MAX_RECEPTORS"]:
                st["receptors"].append({"rel": rel, "display": Path(rel).name, "status": "new"})
        _save_state(ws, st)
        return jsonify({"count": len(st["receptors"]), "receptors": st["receptors"]})

    @app.post("/api/receptors/fetch")
    @login_required
    def api_receptors_fetch():
        jobname = request.form.get("jobname", "")
        pdbid = (request.form.get("pdb") or "").strip()
        chains = (request.form.get("chains") or "").strip()
        if not pdbid:
            return ("missing pdb", 400)
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        try:
            out = fetch_pdb_and_prep(pdbid, rec_dir, chains=chains)
        except ValueError as exc:
            return (str(exc), 400)
        rel = str(Path("Receptors") / Path(out["pdb_path"]).name)

        st = _load_state(ws)
        have = {r["rel"] for r in st["receptors"]}
        if rel not in have and len(st["receptors"]) < current_app.config["MAX_RECEPTORS"]:
            st["receptors"].append({"rel": rel, "display": Path(rel).name, "status": "new"})
        _save_state(ws, st)
        return jsonify({"rel": rel, "count": len(st["receptors"])})

    @app.get("/api/receptors/list")
    @login_required
    def api_receptors_list():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        return jsonify(_load_state(ws))

    # Serve a file from workspace (for NGL)
    @app.get("/api/wsfile")
    @login_required
    def api_wsfile():
        jobname = request.args.get("jobname", "")
        rel = request.args.get("rel", "")
        ws = _ws(jobname)
        if not jobname or not rel or not ws.exists():
            return ("missing", 400)
        p = _resolve_workspace_file(ws, rel)
        if p is None and rel.lower().startswith("receptors/") and rel.lower().endswith(".pdb"):
            try:
                p = _ensure_receptor_pdb_snapshot(ws, rel)
            except Exception:
                p = None
        if p is None:
            return ("not found", 404)
        return send_file(p)

    @app.get("/api/wsinline")
    @login_required
    def api_wsinline():
        jobname = request.args.get("jobname", "")
        rel = request.args.get("rel", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 404)
        p = _resolve_workspace_file(ws, rel)
        if p is None and rel.lower().startswith("receptors/") and rel.lower().endswith(".pdb"):
            try:
                p = _ensure_receptor_pdb_snapshot(ws, rel)
            except Exception:
                p = None
        if p is None or not _path_within(ws, p):
            return ("not found", 404)
        return send_file(p, as_attachment=False, download_name=p.name)

    # ---------- HET COUNTS ----------
    @app.get("/api/het_counts")
    @login_required
    def api_het_counts():
        jobname = request.args.get("jobname","")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        st = _load_state(ws)

        std = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
               "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
               "SEC","PYL","MSE"}
        counts: Dict[str,int] = {}
        for rec in st.get("receptors", []):
            rel = rec.get("rel","")
            p = (ws / rel).resolve()
            if not p.exists() or not p.is_file():
                continue
            with open(p) as f:
                for line in f:
                    if line.startswith(("ATOM","HETATM")):
                        rn = line[17:20].strip().upper()
                        if rn not in std and line.startswith("HETATM"):
                            counts[rn] = counts.get(rn, 0) + 1
        return jsonify({"het_counts": counts})

    # ---------- CENTER ----------
    @app.post("/api/receptor/center")
    @login_required
    def api_receptor_center():
        data = request.get_json(force=True)
        jobname = data.get("jobname")
        rel = data.get("rel")
        size = float(data.get("size") or 20.0)
        method = (data.get("method") or "").lower()
        ws = _ws(jobname)
        src = _resolve_workspace_file(ws, rel or "")
        if src is None:
            return ("not found", 404)

        center: Optional[Tuple[float, float, float]] = None

        def _avg(pts: List[Tuple[float, float, float]]):
            n = max(1, len(pts))
            return (sum(p[0] for p in pts) / n,
                    sum(p[1] for p in pts) / n,
                    sum(p[2] for p in pts) / n)

        if method == "xyz":
            arr = data.get("center") or []
            if len(arr) == 3:
                center = (float(arr[0]), float(arr[1]), float(arr[2]))

        elif method == "selection_atoms":
            pts = [(float(x), float(y), float(z)) for x, y, z in (data.get("atoms") or [])]
            if pts:
                center = _avg(pts)

        else:
            payload = {"method": method, "size": size}
            if method == "ligand":
                payload["ligand"] = data.get("lig_resname")
                payload["chain"] = data.get("chain")
                payload["resi"] = data.get("resi")
            elif method == "residue":
                payload["resname"] = data.get("resname")
                payload["chain"] = data.get("chain")
                payload["resi"] = data.get("resi")
            try:
                result = resolve_center_from_file(src, payload)
                center = tuple(result["center"])
            except CenterResolutionError as exc:
                return (exc.message, exc.status_code)

        if center is None:
            return ("could not compute center", 400)

        pdbqt_name = Path(rel).name
        pdbqt_name = (pdbqt_name.replace(".pdb", ".pdbqt")
                                  .replace(".cif", ".pdbqt")
                                  .replace(".mmcif", ".pdbqt")
                                  .replace(".ent", ".pdbqt"))
        st = _load_state(ws)
        _upsert_center_row(ws, st, pdbqt_name, center, size)

        # status becomes CENTERED here; only becomes PREPPED after 3a
        for r in st["receptors"]:
            if r["rel"] == rel:
                r["status"] = "centered"
        _save_state(ws, st)

        # all receptors centered?
        expected = []
        for r in st.get("receptors", []):
            nm = Path(r["rel"]).name
            nm = nm.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt")
            expected.append(nm)
        csv_map = _read_centers(ws, st)
        all_centered = (len(expected) > 0 and all(n in csv_map for n in expected))

        return jsonify({"center": center, "size": size,
                        "csv": _centers_csv_path(ws, st).name,
                        "all_centered": all_centered})

    # ---------- helpers for converted detection ----------
    def _stem(p: Path) -> str:
        # basename without suffix; "HDGF_7hg9.converted.pdbqt" -> "HDGF_7hg9.converted"
        return p.name[:-len(p.suffix)] if p.suffix else p.name

    def _expected_stem(rel_name: str) -> str:
        # "HDGF_7hg9.pdb" -> "HDGF_7hg9"
        return Path(rel_name).with_suffix("").name

    def _file_matches_receptor(receptor_rel: str, pdbqt_file: Path) -> bool:
        """
        True if pdbqt filename looks like it was generated from the receptor.
        Accepts NAME.pdbqt, NAME.converted.pdbqt, NAME.anything.pdbqt
        where NAME is the receptor base (without original extension).
        """
        base = _expected_stem(receptor_rel)
        return _stem(pdbqt_file).startswith(base)

    def _mark_prepped_from_output(ws: Path, st: Dict[str, Any], out_dir: Path):
        """Update in-memory state if corresponding PDBQT exists."""
        changed = False
        files = list(out_dir.glob("*.pdbqt"))
        for r in st.get("receptors", []):
            if r.get("status") == "prepped":
                continue
            for f in files:
                if _file_matches_receptor(Path(r["rel"]).name, f):
                    r["status"] = "prepped"
                    changed = True
                    break
        if changed:
            _save_state(ws, st)

    # ---------- SUMMARY ----------
    @app.get("/api/summary")
    @login_required
    def api_summary():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        st = _load_state(ws)

        out_dir = (ensure_subdir(ws, "Receptors").parent / "Receptors_PDBQT").resolve()
        if out_dir.exists():
            _mark_prepped_from_output(ws, st, out_dir)
            st = _load_state(ws)  # re-read after potential update

        total = len(st.get("receptors", []))
        centered = sum(1 for r in st.get("receptors", []) if r.get("status") in ("centered","prepped"))
        prepped  = sum(1 for r in st.get("receptors", []) if r.get("status") == "prepped")
        csv_map = _read_centers(ws, st)
        csv_rows = len(csv_map)

        expected = []
        for r in st.get("receptors", []):
            nm = Path(r["rel"]).name
            nm = nm.replace(".pdb", ".pdbqt").replace(".cif", ".pdbqt").replace(".mmcif", ".pdbqt").replace(".ent", ".pdbqt")
            expected.append(nm)

        converted = [p.name for p in sorted(out_dir.glob("*.pdbqt"))] if out_dir.exists() else []
        centers = [
            {"receptor_pdbqt": key, "center": [x, y, z], "size": size}
            for key, (x, y, z, size) in sorted(csv_map.items())
        ]

        all_receptors_centered = (total > 0 and centered == total)
        have_rows_for_all = all(n in csv_map for n in expected) if expected else False

        return jsonify({
            "jobname": jobname,
            "receptors_total": total,
            "receptors_centered": centered,
            "receptors_prepped": prepped,
            "centers_csv": str(_centers_csv_path(ws, st)),
            "centers_rows": csv_rows,
            "expected_pdbqt": expected,
            "have_rows_for_all": have_rows_for_all,
            "conversion_out_dir": str(out_dir),
            "converted_count": len(converted),
            "converted_list": converted,
            "all_receptors_centered": all_receptors_centered,
            "centers": centers,
        })

    # ---------- 3a CONVERSION (batch) ----------


    @app.post("/api/prep/start")
    @login_required
    def api_prep_start():
        f = request.form
        jobname = f.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        st = _load_state(ws)
        recs = st.get("receptors", [])

        if not recs:
            return ("Add receptors first.", 400)

        csv_map = _read_centers(ws, st)
        expected = [Path(r["rel"]).name for r in recs]
        expected = [
            n.replace(".pdb", ".pdbqt")
             .replace(".cif", ".pdbqt")
             .replace(".mmcif", ".pdbqt")
             .replace(".ent", ".pdbqt")
            for n in expected
        ]

        if not all(n in csv_map for n in expected):
            return ("Centers CSV missing one or more receptors. Save a center for each receptor first.", 400)

        remove_hets = []
        remove_all_hets = False

        if f.get("remove_het_csv"):
            if f.get("remove_het_csv").lower() == "all":
                remove_all_hets = True
            else:
                remove_hets = [
                    x.strip().upper()
                    for x in f.get("remove_het_csv").split(",")
                    if x.strip()
                ]

        remove_chains = [
            x.strip().upper()
            for x in f.get("remove_chains_csv", "").split(",")
            if x.strip()
        ]

        altloc_mode = f.get("altloc", "collapse")

        out_dir = (rec_dir.parent / "Receptors_PDBQT").resolve()
        log_path = (ws / "prep3a.log").resolve()
        out_dir.mkdir(exist_ok=True)

        with open(log_path, "w") as logf:
            for r in recs:
                in_path = _ensure_receptor_pdb_snapshot(ws, r["rel"], altloc_mode=altloc_mode)
                cleaned_path = ws / f"{Path(r['rel']).stem}_clean.pdb"
                out_path = out_dir / Path(r["rel"]).with_suffix(".pdbqt").name

                try:
                    _clean_pdb(
                        in_path,
                        cleaned_path,
                        remove_hets,
                        remove_chains,
                        remove_all_hets,
                        altloc_mode,
                    )
                except Exception as exc:
                    logf.write(f"[ERROR] Failed to clean {in_path}: {exc}\n")
                    st["prep_job"] = {
                        "pid": None,
                        "log": str(log_path),
                        "out_dir": str(out_dir),
                    }
                    _save_state(ws, st)
                    return jsonify({
                        "ok": False,
                        "error": "receptor_cleaning_failed",
                        "message": f"Failed to clean receptor {Path(r['rel']).name}.",
                        "details": {
                            "receptor": r.get("rel"),
                            "input": str(in_path),
                            "log": str(log_path),
                            "exception": str(exc),
                        },
                    }), 500

                cmd = ["obabel", str(cleaned_path), "-O", str(out_path), "-xr"]
                logf.write(f"Running: {' '.join(cmd)}\n")

                try:
                    result = subprocess.run(
                        cmd,
                        stdout=logf,
                        stderr=subprocess.STDOUT,
                    )
                except FileNotFoundError:
                    logf.write("[ERROR] obabel executable was not found on PATH.\n")
                    st["prep_job"] = {
                        "pid": None,
                        "log": str(log_path),
                        "out_dir": str(out_dir),
                    }
                    _save_state(ws, st)
                    return jsonify({
                        "ok": False,
                        "error": "obabel_missing",
                        "message": "Open Babel CLI executable `obabel` is not installed or not available on PATH.",
                        "details": {
                            "cmd": cmd,
                            "log": str(log_path),
                            "hint": "Install Open Babel on Heroku using an Aptfile or another buildpack/system dependency method.",
                        },
                    }), 500
                except Exception as exc:
                    logf.write(f"[ERROR] Failed to run obabel: {exc}\n")
                    st["prep_job"] = {
                        "pid": None,
                        "log": str(log_path),
                        "out_dir": str(out_dir),
                    }
                    _save_state(ws, st)
                    return jsonify({
                        "ok": False,
                        "error": "obabel_execution_failed",
                        "message": "Open Babel execution failed unexpectedly.",
                        "details": {
                            "cmd": cmd,
                            "log": str(log_path),
                            "exception": str(exc),
                        },
                    }), 500

                if result.returncode != 0 or not out_path.exists():
                    logf.write(f"[ERROR] Failed to convert {in_path}\n")
                    st["prep_job"] = {
                        "pid": None,
                        "log": str(log_path),
                        "out_dir": str(out_dir),
                    }
                    _save_state(ws, st)
                    return jsonify({
                        "ok": False,
                        "error": "receptor_conversion_failed",
                        "message": f"Failed to convert receptor {Path(r['rel']).name} to PDBQT.",
                        "details": {
                            "cmd": cmd,
                            "returncode": result.returncode,
                            "input": str(in_path),
                            "cleaned": str(cleaned_path),
                            "output": str(out_path),
                            "log": str(log_path),
                        },
                    }), 500

        st["prep_job"] = {
            "pid": None,
            "log": str(log_path),
            "out_dir": str(out_dir),
        }
        _save_state(ws, st)

        return jsonify({
            "done": True,
            "out_dir": str(out_dir),
            "log": str(log_path),
        })

    
    # ---------- 3a CONVERSION (single) ----------
    @app.post("/api/prep/start_one")
    @login_required
    def api_prep_start_one():
        f = request.form
        jobname = f.get("jobname", "")
        rel = f.get("rel", "")
        ws = _ws(jobname)

        if not ws.exists() or not rel:
            return ("workspace or rel missing", 400)

        rec_dir = ensure_subdir(ws, "Receptors")
        st = _load_state(ws)

        csv_map = _read_centers(ws, st)
        name = Path(rel).name
        expected = (
            name.replace(".pdb", ".pdbqt")
                .replace(".cif", ".pdbqt")
                .replace(".mmcif", ".pdbqt")
                .replace(".ent", ".pdbqt")
        )

        if expected not in csv_map:
            return ("Save a center for this receptor first.", 400)

        remove_hets = []
        remove_all_hets = False

        if f.get("remove_het_csv"):
            if f.get("remove_het_csv").lower() == "all":
                remove_all_hets = True
            else:
                remove_hets = [
                    x.strip().upper()
                    for x in f.get("remove_het_csv").split(",")
                    if x.strip()
                ]

        remove_chains = [
            x.strip().upper()
            for x in f.get("remove_chains_csv", "").split(",")
            if x.strip()
        ]

        altloc_mode = f.get("altloc", "collapse")

        out_dir = (rec_dir.parent / "Receptors_PDBQT").resolve()
        log_path = (ws / "prep3a.log").resolve()
        out_dir.mkdir(exist_ok=True)

        in_path = _ensure_receptor_pdb_snapshot(ws, rel, altloc_mode=altloc_mode)
        cleaned_path = ws / f"{Path(name).stem}_clean.pdb"
        out_path = out_dir / Path(name).with_suffix(".pdbqt").name

        try:
            _clean_pdb(
                in_path,
                cleaned_path,
                remove_hets,
                remove_chains,
                remove_all_hets,
                altloc_mode,
            )
        except Exception as exc:
            with open(log_path, "a") as logf:
                logf.write(f"[ERROR] Failed to clean {in_path}: {exc}\n")

            st["prep_job"] = {
                "pid": None,
                "log": str(log_path),
                "out_dir": str(out_dir),
            }
            _save_state(ws, st)

            return jsonify({
                "ok": False,
                "error": "receptor_cleaning_failed",
                "message": f"Failed to clean receptor {name}.",
                "details": {
                    "receptor": rel,
                    "input": str(in_path),
                    "log": str(log_path),
                    "exception": str(exc),
                },
            }), 500

        with open(log_path, "a") as logf:
            cmd = ["obabel", str(cleaned_path), "-O", str(out_path), "-xr"]
            logf.write(f"Running: {' '.join(cmd)}\n")

            try:
                result = subprocess.run(
                    cmd,
                    stdout=logf,
                    stderr=subprocess.STDOUT,
                )
            except FileNotFoundError:
                logf.write("[ERROR] obabel executable was not found on PATH.\n")
                st["prep_job"] = {
                    "pid": None,
                    "log": str(log_path),
                    "out_dir": str(out_dir),
                }
                _save_state(ws, st)

                return jsonify({
                    "ok": False,
                    "error": "obabel_missing",
                    "message": "Open Babel CLI executable `obabel` is not installed or not available on PATH.",
                    "details": {
                        "cmd": cmd,
                        "log": str(log_path),
                        "hint": "Install Open Babel on Heroku using an Aptfile or another buildpack/system dependency method.",
                    },
                }), 500
            except Exception as exc:
                logf.write(f"[ERROR] Failed to run obabel: {exc}\n")
                st["prep_job"] = {
                    "pid": None,
                    "log": str(log_path),
                    "out_dir": str(out_dir),
                }
                _save_state(ws, st)

                return jsonify({
                    "ok": False,
                    "error": "obabel_execution_failed",
                    "message": "Open Babel execution failed unexpectedly.",
                    "details": {
                        "cmd": cmd,
                        "log": str(log_path),
                        "exception": str(exc),
                    },
                }), 500

            if result.returncode != 0 or not out_path.exists():
                logf.write(f"[ERROR] Failed to convert {in_path}\n")
                st["prep_job"] = {
                    "pid": None,
                    "log": str(log_path),
                    "out_dir": str(out_dir),
                }
                _save_state(ws, st)

                return jsonify({
                    "ok": False,
                    "error": "receptor_conversion_failed",
                    "message": f"Failed to convert receptor {name} to PDBQT.",
                    "details": {
                        "cmd": cmd,
                        "returncode": result.returncode,
                        "input": str(in_path),
                        "cleaned": str(cleaned_path),
                        "output": str(out_path),
                        "log": str(log_path),
                    },
                }), 500

        st["prep_job"] = {
            "pid": None,
            "log": str(log_path),
            "out_dir": str(out_dir),
        }
        _save_state(ws, st)

        return jsonify({
            "done": True,
            "out": str(out_path),
            "log": str(log_path),
        })



    # ---------- 3a STATUS ----------
    @app.get("/api/prep/status")
    @login_required
    def api_prep_status():
        jobname = request.args.get("jobname", "")
        ws = _ws(jobname)
        st = _load_state(ws)
        info = st.get("prep_job") or {}
        if not info:
            return jsonify({"running": False, "done": False, "log": ""})

        pid = info.get("pid")
        try:
            os.kill(pid, 0)
            running = True
        except Exception:
            running = False

        lp = Path(info.get("log", ""))
        log = lp.read_text()[-8000:] if lp.exists() else ""
        out_dir = Path(info.get("out_dir", "")).resolve()

        # Determine "done": each receptor must have at least one matching *.pdbqt
        st_now = _load_state(ws)
        receptor_names = [Path(r["rel"]).name for r in st_now.get("receptors", [])]
        files = list(out_dir.glob("*.pdbqt")) if out_dir.exists() else []
        matched = 0
        for rel in receptor_names:
            if any(_file_matches_receptor(rel, f) for f in files):
                matched += 1
        done = (matched == len(receptor_names) and matched > 0)

        # Mark prepped immediately when matching files appear
        if files:
            _mark_prepped_from_output(ws, st_now, out_dir)

        return jsonify({"running": running, "done": done, "log": log})

    # ---------- LIGANDS ----------
    @app.post("/api/ligands/upload")
    @login_required
    def api_lig_upload():
        jobname = request.form.get("jobname", "")
        mode = request.form.get("mode")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        lig_dir = ensure_subdir(ws, "Ligands")
        result: Dict[str, Any]

        if mode == "single":
            f = request.files.get("file")
            if not f:
                return ("no file", 400)
            result = save_uploaded_ligand_folder([f], lig_dir)
            upload_name = Path(f.filename).name
        elif mode == "zip":
            f = request.files.get("file")
            if not f:
                return ("no file", 400)
            result = save_uploaded_ligand_zip(f, lig_dir)
            upload_name = Path(f.filename).name
        elif mode == "folder":
            files = request.files.getlist("files")
            if not files:
                return ("no files", 400)
            result = save_uploaded_ligand_folder(files, lig_dir)
            raw_folder = request.form.get("folder_name", "").strip()
            upload_name = raw_folder or Path((files[0].filename or "").split("/", 1)[0]).name or "Ligands"
        else:
            return ("bad mode", 400)

        if not result["accepted_count"]:
            return jsonify({
                "ok": False,
                "error": "No supported ligand files were found. Upload .sdf, .smiles, .smi, or .csv files.",
                **result,
            }), 400

        csv_meta = _csv_upload_metadata(lig_dir, result)
        st = _load_state(ws)
        st["ligands_uploaded"] = True
        st["ligand_info"] = {
            "upload_mode": mode,
            "filename": upload_name,
            "accepted_files": result["accepted_files"],
            "accepted_count": result["accepted_count"],
            "ignored_files": result["ignored_files"],
            "warnings": result["warnings"],
            "is_csv": result["is_csv"],
            "filetypes": result["filetypes"],
            "ligands_root": result["ligands_root"],
            **csv_meta,
        }
        _save_state(ws, st)
        return jsonify({"ok": True, **st["ligand_info"]})

    @app.post("/api/ligands/curated")
    @login_required
    def api_lig_curated():
        jobname = request.form.get("jobname", "")
        library_key = request.form.get("library", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)
        curated = _curated_ligand_source(library_key)
        if not curated:
            return jsonify({"ok": False, "error": "Curated ligand library not found."}), 404

        lig_dir = ensure_subdir(ws, "Ligands")
        source_path = curated["path"]
        target_path = lig_dir / source_path.name
        shutil.copy2(source_path, target_path)

        result = {
            "accepted_files": [target_path.name],
            "accepted_count": 1,
            "ignored_files": [],
            "warnings": [],
            "is_csv": True,
            "filetypes": [".csv"],
            "ligands_root": str(lig_dir),
            "headers": curated["headers"],
            "suggested_smiles_col": curated["suggested_smiles_col"],
            "suggested_id_col": curated["suggested_id_col"],
            "curated_library": {
                "key": curated["key"],
                "label": curated["label"],
                "description": curated["description"],
            },
        }

        st = _load_state(ws)
        st["ligands_uploaded"] = True
        st["ligand_info"] = {
            "upload_mode": "curated",
            "filename": source_path.name,
            **result,
        }
        _save_state(ws, st)
        return jsonify({"ok": True, **st["ligand_info"]})

    # ---------- BUILD ----------
    @app.post("/api/build")
    @login_required
    def api_build():
        f = request.form
        jobname = f.get("jobname", "")
        ws = _ws(jobname)
        if not ws.exists():
            return ("workspace missing", 400)

        package_mode = normalize_package_mode(
            f,
            default_mode=current_app.config["DEFAULT_PACKAGE_MODE"],
            lsf_enabled=current_app.config["ENABLE_LSF_PACKAGE"],
        )
        if package_mode in {"joey_lsf", "mainak_lsf", "custom_lsf"} and not current_app.config["ENABLE_LSF_PACKAGE"]:
            return ("lsf packaging is disabled", 400)

        st = _load_state(ws)
        if not st["receptors"]:
            return ("add receptors first", 400)

        # Verify centers CSV covers all expected receptors
        csv_map = _read_centers(ws, st)
        expected = [Path(r["rel"]).name for r in st["receptors"]]
        expected = [
            n.replace(".pdb", ".pdbqt")
            .replace(".cif", ".pdbqt")
            .replace(".mmcif", ".pdbqt")
            .replace(".ent", ".pdbqt")
            for n in expected
        ]
        if not all(n in csv_map for n in expected):
            return ("centers csv is incomplete", 400)

        if not st.get("ligands_uploaded"):
            return ("upload ligands first", 400)

        ligand_info = st.get("ligand_info") or {}

        # Build the job tree (copies real scripts, receptors, ligands, etc.)
        jobroot, warnings = assemble_job_tree(ws, ws / "Receptors", ws / "Ligands", package_mode=package_mode)
        lsf_dir = jobroot  # write .lsf files directly into jobroot

        # ---------- Normalize inputs from the form (no frontend changes required) ----------
        def _get_int(*keys, default):
            for k in keys:
                v = f.get(k, "").strip()
                if v:
                    try:
                        return int(v)
                    except Exception:
                        pass
            return int(default)

        def _get_str(*keys, default=""):
            for k in keys:
                v = f.get(k)
                if v is not None and v.strip() != "":
                    return v.strip()
            return default

        # ConfGen parameters
        poses_conf        = _get_int("confgen_poses", "poses_conf", "poses", default=64)
        inferred_lig_mode, inferred_lig_filetype, inferred_single_sdf = infer_ligand_workflow(ligand_info)
        lig_mode          = _get_str("lig_mode", default=inferred_lig_mode)  # "1" CSV, "2" folder, "3" single SDF
        lig_filetype      = _get_str("lig_filetype", "filetype", default=inferred_lig_filetype)
        csv_smiles_col    = _get_str("csv_smiles_col", "smiles_col", default="")
        csv_id_col        = _get_str("csv_id_col", "id_col", default="")
        single_sdf_rel    = _get_str("single_sdf_rel", "single_sdf", default=inferred_single_sdf or "")

        if lig_mode == "1" and not csv_smiles_col:
            return ("select a CSV SMILES column before building", 400)

        # Vina parameters
        poses_vina        = _get_int("vina_poses", "poses_vina", default=20)
        profile = profile_for_mode(
            package_mode,
            {
                **f,
                "lsf_email": _get_str("lsf_email", "email", default=""),
                "notify_begin": f.get("notify_begin"),
                "notify_end": f.get("notify_end"),
                "queue": _get_str("queue", default="general"),
                "project": _get_str("project", default=""),
                "workers": _get_int("workers", "ncores", default=16),
                "mem_per_core": _get_int("mem_per_core", "mem", default=2000),
                "confgen_walltime": _get_str("confgen_walltime", "walltime_confgen", default="48:00"),
                "vina_walltime": _get_str("vina_walltime", "walltime_vina", "walltime", default="96:00"),
                "conda_sh": _get_str("conda_sh", default=""),
                "conda_env": _get_str("conda_env", default=""),
                "vina_path": _get_str("vina_path", "vina_exe", default=""),
                "python_command": _get_str("python_command", default='$(command -v python3 || command -v python)'),
                "setup_commands": f.get("setup_commands") or _get_str("env_line", default=current_app.config["DEFAULT_ENV_LINE"]),
            },
        )

        if package_mode in {"joey_lsf", "mainak_lsf", "custom_lsf"}:
            build_confgen_lsfs(
                jobroot, lsf_dir, profile=profile, poses=poses_conf,
                lig_mode=lig_mode, lig_filetype=lig_filetype,
                csv_smiles_col=csv_smiles_col, csv_id_col=csv_id_col,
                single_sdf_rel=(single_sdf_rel or None)
            )

            build_vina_lsfs(
                jobroot, lsf_dir, profile=profile, poses=poses_vina
            )
        else:
            build_portable_runners(jobroot)

        # Optional tag-based rename (no-op unless TAG column exists)
        rename_centers_with_tags(jobroot)

        # Zip the job folder
        z = zip_job_tree(jobroot)
        return jsonify(
            {
                "zip": str(z),
                "download_url": url_for("download", path=str(z)),
                "package_mode": package_mode,
                "warnings": warnings,
            }
        )

    # ---------- VERSIONED HEADLESS API ----------
    def _v1_ok(data: Optional[Dict[str, Any]] = None, warnings: Optional[List[str]] = None, status: int = 200):
        return jsonify({"ok": True, "data": data or {}, "warnings": warnings or []}), status

    def _v1_error(error: str, message: str, status: int = 400, details: Optional[Dict[str, Any]] = None):
        return jsonify({"ok": False, "error": error, "message": message, "details": details or {}}), status

    def _json_payload() -> Dict[str, Any]:
        if request.is_json:
            return request.get_json(silent=True) or {}
        return dict(request.form.items())

    def _sanitize_workspace_name(value: str) -> str:
        value = re.sub(r"[^A-Za-z0-9._-]+", "-", (value or "").strip()).strip(".-_")
        return value[:80] or _public_name()

    def _new_jobname(requested: str = "", reuse: bool = False) -> Tuple[str, Path, bool]:
        safe = _sanitize_workspace_name(requested)
        if requested:
            candidate = safe
            ws = _ws(candidate)
            if reuse and ws.exists():
                return candidate, ws, True
            if not ws.exists():
                return candidate, make_workspace(ws), False
        stamp = time.strftime("%m-%d-%Y-%H-%M-%S")
        base = f"{stamp}-{safe}"
        candidate = base
        idx = 2
        while _ws(candidate).exists():
            candidate = f"{base}-{idx}"
            idx += 1
        return candidate, make_workspace(_ws(candidate)), False

    def _initial_state() -> Dict[str, Any]:
        return {"receptors": [], "centers_csv": "vina_centers.csv",
                "ligands_uploaded": False, "ligand_info": {}, "prep_job": None}

    def _receptor_pdbqt_name(rel: str) -> str:
        return (Path(rel).name.replace(".pdb", ".pdbqt")
                .replace(".cif", ".pdbqt")
                .replace(".mmcif", ".pdbqt")
                .replace(".ent", ".pdbqt"))

    def _resolve_receptor_for_api(ws: Path, st: Dict[str, Any], payload: Dict[str, Any]) -> Tuple[Optional[str], Optional[Path]]:
        requested = (payload.get("receptor") or payload.get("rel") or "").strip()
        receptors = st.get("receptors", [])
        rel = ""
        if requested:
            for rec in receptors:
                if requested in {rec.get("rel"), rec.get("display"), Path(rec.get("rel", "")).name}:
                    rel = rec.get("rel", "")
                    break
            rel = rel or requested
        elif len(receptors) == 1:
            rel = receptors[0].get("rel", "")
        if not rel:
            return None, None
        return rel, _resolve_workspace_file(ws, rel)

    def _register_receptors(ws: Path, st: Dict[str, Any], added: List[str]) -> Dict[str, Any]:
        have = {r["rel"] for r in st["receptors"]}
        for rel in added:
            if rel not in have and len(st["receptors"]) < current_app.config["MAX_RECEPTORS"]:
                st["receptors"].append({"rel": rel, "display": Path(rel).name, "status": "new"})
                have.add(rel)
        _save_state(ws, st)
        return {"count": len(st["receptors"]), "receptors": st["receptors"]}

    def _summary_data(jobname: str, ws: Path) -> Dict[str, Any]:
        st = _load_state(ws)
        out_dir = (ensure_subdir(ws, "Receptors").parent / "Receptors_PDBQT").resolve()
        if out_dir.exists():
            _mark_prepped_from_output(ws, st, out_dir)
            st = _load_state(ws)
        total = len(st.get("receptors", []))
        centered = sum(1 for r in st.get("receptors", []) if r.get("status") in ("centered", "prepped"))
        prepped = sum(1 for r in st.get("receptors", []) if r.get("status") == "prepped")
        csv_map = _read_centers(ws, st)
        expected = [_receptor_pdbqt_name(r["rel"]) for r in st.get("receptors", [])]
        converted = [p.name for p in sorted(out_dir.glob("*.pdbqt"))] if out_dir.exists() else []
        centers = [
            {"receptor_pdbqt": key, "center": [x, y, z], "size": size}
            for key, (x, y, z, size) in sorted(csv_map.items())
        ]
        zips = sorted(ws.glob("*.zip"), key=lambda p: p.stat().st_mtime, reverse=True)
        return {
            "jobname": jobname,
            "workspace": str(ws),
            "receptors": st.get("receptors", []),
            "receptors_total": total,
            "receptors_centered": centered,
            "receptors_prepped": prepped,
            "ligands_uploaded": bool(st.get("ligands_uploaded")),
            "ligand_info": st.get("ligand_info") or {},
            "centers_csv": str(_centers_csv_path(ws, st)),
            "centers_rows": len(csv_map),
            "expected_pdbqt": expected,
            "have_rows_for_all": all(n in csv_map for n in expected) if expected else False,
            "prep_running": False,
            "prep_done": prepped == total and total > 0,
            "package_ready": bool(zips),
            "artifacts": [
                {"name": p.name, "path": str(p), "download_url": url_for("download", path=str(p))}
                for p in zips
            ],
            "conversion_out_dir": str(out_dir),
            "converted_count": len(converted),
            "converted_list": converted,
            "centers": centers,
        }

    def _resolve_center_payload(jobname: str, payload: Dict[str, Any]) -> Tuple[Dict[str, Any], Optional[str], Optional[Path]]:
        ws = _ws(jobname)
        if not ws.exists():
            raise CenterResolutionError("workspace_missing", f"Workspace {jobname} does not exist.", status_code=404)
        st = _load_state(ws)
        method = (payload.get("method") or "xyz").strip().lower()
        rel, src = _resolve_receptor_for_api(ws, st, payload)
        if method != "xyz":
            if src is None:
                raise CenterResolutionError("receptor_not_found", "A matching receptor file was not found in the workspace.", {"receptor": payload.get("receptor")}, 404)
            result = resolve_center_from_file(src, payload)
        else:
            result = resolve_center_from_file(src, payload) if src is not None else resolve_xyz(payload)
        result["jobname"] = jobname
        if rel:
            result["receptor_rel"] = rel
        return result, rel, src

    @app.get("/api/v1/health")
    @login_required
    def api_v1_health():
        return _v1_ok({
            "service": "autodock-vina-prepserver",
            "api_version": "v1",
            "public_mode": current_app.config.get("PUBLIC_MODE", True),
            "tmp_root": current_app.config["TMP_ROOT"],
        })

    @app.post("/api/v1/workspaces")
    @login_required
    def api_v1_workspaces_create():
        payload = _json_payload()
        jobname, ws, reused = _new_jobname(
            payload.get("workspace_name") or payload.get("jobname") or "",
            bool(payload.get("reuse")),
        )
        if not reused:
            _save_state(ws, _initial_state())
        return _v1_ok({"jobname": jobname, "workspace": str(ws), "reused": reused}, status=201 if not reused else 200)

    @app.get("/api/v1/workspaces/<jobname>")
    @login_required
    def api_v1_workspace_get(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        return _v1_ok({"jobname": jobname, "workspace": str(ws), "state": _load_state(ws)})

    @app.get("/api/v1/workspaces/<jobname>/summary")
    @login_required
    def api_v1_workspace_summary(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        return _v1_ok(_summary_data(jobname, ws))

    @app.post("/api/v1/workspaces/<jobname>/receptors/upload")
    @login_required
    def api_v1_receptors_upload(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        mode = request.form.get("mode") or "single"
        rec_dir = ensure_subdir(ws, "Receptors")
        added: List[str] = []
        if mode == "folder":
            files = request.files.getlist("files")
            if not files:
                return _v1_error("missing_file", "No files were uploaded.", 400)
            for file_storage in files:
                out = rec_dir / Path(file_storage.filename or "").name
                file_storage.save(out)
                if out.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                    added.append(str(Path("Receptors") / out.name))
        else:
            f = request.files.get("file")
            if not f:
                return _v1_error("missing_file", "No receptor file was uploaded.", 400)
            if mode == "zip":
                save_uploaded_zip(f, rec_dir)
                for p in sorted(rec_dir.rglob("*")):
                    if p.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                        added.append(str(Path("Receptors") / p.relative_to(rec_dir)))
            elif mode == "single":
                out = rec_dir / Path(f.filename).name
                f.save(out)
                if out.suffix.lower() in {".pdb", ".pdbqt", ".cif", ".mmcif", ".ent"}:
                    added.append(str(Path("Receptors") / out.name))
            else:
                return _v1_error("bad_mode", "Receptor upload mode must be single, zip, or folder.", 400)
        return _v1_ok(_register_receptors(ws, _load_state(ws), added))

    @app.post("/api/v1/workspaces/<jobname>/receptors/fetch")
    @login_required
    def api_v1_receptors_fetch(jobname: str):
        payload = _json_payload()
        pdbid = (payload.get("pdb_id") or payload.get("pdb") or "").strip()
        chains = payload.get("chains") or ""
        if isinstance(chains, list):
            chains = ",".join(str(c) for c in chains)
        chains = str(chains).strip()
        if not pdbid:
            return _v1_error("missing_pdb_id", "Provide pdb_id.", 400)
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        try:
            out = fetch_pdb_and_prep(pdbid, ensure_subdir(ws, "Receptors"), chains=chains)
        except ValueError as exc:
            return _v1_error("receptor_fetch_failed", str(exc), 400)
        rel = str(Path("Receptors") / Path(out["pdb_path"]).name)
        data = _register_receptors(ws, _load_state(ws), [rel])
        data["rel"] = rel
        return _v1_ok(data)

    @app.get("/api/v1/workspaces/<jobname>/receptors")
    @login_required
    def api_v1_receptors_list(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        return _v1_ok({"receptors": _load_state(ws).get("receptors", [])})

    @app.post("/api/v1/workspaces/<jobname>/centers/resolve")
    @login_required
    def api_v1_centers_resolve(jobname: str):
        try:
            result, _rel, _src = _resolve_center_payload(jobname, _json_payload())
            return _v1_ok(result)
        except CenterResolutionError as exc:
            return _v1_error(exc.error, exc.message, exc.status_code, exc.details)

    @app.post("/api/v1/workspaces/<jobname>/centers/save")
    @login_required
    def api_v1_centers_save(jobname: str):
        ws = _ws(jobname)
        try:
            payload = _json_payload()
            result, rel, _src = _resolve_center_payload(jobname, payload)
            st = _load_state(ws)
            if not rel:
                rel, _src = _resolve_receptor_for_api(ws, st, payload)
            if not rel:
                return _v1_error("receptor_required", "Saving a center requires receptor when the workspace has zero or multiple receptors.", 400)
            _upsert_center_row(ws, st, _receptor_pdbqt_name(rel), tuple(result["center"]), result["size"])
            for receptor in st.get("receptors", []):
                if receptor.get("rel") == rel:
                    receptor["status"] = "centered"
            _save_state(ws, st)
            result["csv"] = _centers_csv_path(ws, st).name
            return _v1_ok(result)
        except CenterResolutionError as exc:
            return _v1_error(exc.error, exc.message, exc.status_code, exc.details)

    @app.get("/api/v1/workspaces/<jobname>/centers")
    @login_required
    def api_v1_centers_list(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        st = _load_state(ws)
        centers = [
            {"receptor_pdbqt": key, "center": [x, y, z], "size": size}
            for key, (x, y, z, size) in _read_centers(ws, st).items()
        ]
        return _v1_ok({"centers": centers, "csv": str(_centers_csv_path(ws, st))})

    @app.post("/api/v1/workspaces/<jobname>/prep/start")
    @login_required
    def api_v1_prep_start(jobname: str):
        payload = _json_payload()
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        st = _load_state(ws)
        recs = st.get("receptors", [])
        if not recs:
            return _v1_error("no_receptors", "Add receptors first.", 400)
        csv_map = _read_centers(ws, st)
        expected = [_receptor_pdbqt_name(r["rel"]) for r in recs]
        if not all(n in csv_map for n in expected):
            return _v1_error("centers_incomplete", "Save a center for each receptor first.", 400, {"expected_pdbqt": expected})
        remove_hets = []
        remove_all_hets = False
        raw_hets = payload.get("remove_het_csv") or payload.get("remove_hets") or ""
        if isinstance(raw_hets, list):
            remove_hets = [str(x).strip().upper() for x in raw_hets if str(x).strip()]
        elif str(raw_hets).lower() == "all":
            remove_all_hets = True
        elif raw_hets:
            remove_hets = [x.strip().upper() for x in str(raw_hets).split(",") if x.strip()]
        raw_chains = payload.get("remove_chains_csv") or payload.get("remove_chains") or ""
        remove_chains = [str(x).strip().upper() for x in raw_chains] if isinstance(raw_chains, list) else [x.strip().upper() for x in str(raw_chains).split(",") if x.strip()]
        altloc_mode = payload.get("altloc", "collapse")
        rec_dir = ensure_subdir(ws, "Receptors")
        out_dir = (rec_dir.parent / "Receptors_PDBQT").resolve()
        log_path = (ws / "prep3a.log").resolve()
        out_dir.mkdir(exist_ok=True)
        with open(log_path, "w") as logf:
            for r in recs:
                in_path = _ensure_receptor_pdb_snapshot(ws, r["rel"], altloc_mode=altloc_mode)
                cleaned_path = ws / f"{Path(r['rel']).stem}_clean.pdb"
                out_path = out_dir / Path(r["rel"]).with_suffix(".pdbqt").name
                _clean_pdb(in_path, cleaned_path, remove_hets, remove_chains, remove_all_hets, altloc_mode)
                cmd = ["obabel", str(cleaned_path), "-O", str(out_path), "-xr"]
                logf.write(f"Running: {' '.join(cmd)}\n")

                try:
                    result = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)
                except FileNotFoundError:
                    logf.write("[ERROR] obabel executable was not found on PATH.\n")
                    st["prep_job"] = {"pid": None, "log": str(log_path), "out_dir": str(out_dir)}
                    _save_state(ws, st)
                    return _v1_error(
                        "obabel_missing",
                        "Open Babel CLI executable `obabel` is not installed or not available on PATH.",
                        500,
                        {
                            "cmd": cmd,
                            "log": str(log_path),
                            "hint": "Install Open Babel on Heroku using an Aptfile or another buildpack/system dependency method."
                        }
                    )

                if result.returncode != 0 or not out_path.exists():
                    logf.write(f"[ERROR] Failed to convert {in_path}\n")
                    st["prep_job"] = {"pid": None, "log": str(log_path), "out_dir": str(out_dir)}
                    _save_state(ws, st)
                    return _v1_error(
                        "receptor_conversion_failed",
                        f"Failed to convert receptor {Path(r['rel']).name} to PDBQT.",
                        500,
                        {
                            "cmd": cmd,
                            "returncode": result.returncode,
                            "log": str(log_path)
                        }
                    )
        st["prep_job"] = {"pid": None, "log": str(log_path), "out_dir": str(out_dir)}
        _save_state(ws, st)
        return _v1_ok({"jobname": jobname, "done": True, "out_dir": str(out_dir), "log": str(log_path)})

    @app.get("/api/v1/workspaces/<jobname>/prep/status")
    @login_required
    def api_v1_prep_status(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        st = _load_state(ws)
        info = st.get("prep_job") or {}
        if not info:
            return _v1_ok({"jobname": jobname, "running": False, "done": False, "log": ""})
        lp = Path(info.get("log", ""))
        log = lp.read_text(errors="replace")[-8000:] if lp.exists() else ""
        out_dir = Path(info.get("out_dir", "")).resolve()
        files = list(out_dir.glob("*.pdbqt")) if out_dir.exists() else []
        receptor_names = [Path(r["rel"]).name for r in st.get("receptors", [])]
        matched = sum(1 for rel in receptor_names if any(_file_matches_receptor(rel, f) for f in files))
        done = matched == len(receptor_names) and matched > 0
        if files:
            _mark_prepped_from_output(ws, st, out_dir)
        return _v1_ok({"jobname": jobname, "running": False, "done": done, "matched_receptors": matched, "log": log})

    @app.post("/api/v1/workspaces/<jobname>/ligands/upload")
    @login_required
    def api_v1_ligands_upload(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        mode = request.form.get("mode") or "single"
        lig_dir = ensure_subdir(ws, "Ligands")
        if mode == "single":
            f = request.files.get("file")
            if not f:
                return _v1_error("missing_file", "No ligand file was uploaded.", 400)
            result = save_uploaded_ligand_folder([f], lig_dir)
            upload_name = Path(f.filename).name
        elif mode == "zip":
            f = request.files.get("file")
            if not f:
                return _v1_error("missing_file", "No ligand ZIP was uploaded.", 400)
            result = save_uploaded_ligand_zip(f, lig_dir)
            upload_name = Path(f.filename).name
        elif mode == "folder":
            files = request.files.getlist("files")
            if not files:
                return _v1_error("missing_file", "No ligand files were uploaded.", 400)
            result = save_uploaded_ligand_folder(files, lig_dir)
            upload_name = request.form.get("folder_name", "").strip() or "Ligands"
        else:
            return _v1_error("bad_mode", "Ligand upload mode must be single, zip, or folder.", 400)
        if not result["accepted_count"]:
            return _v1_error("no_supported_ligands", "No supported ligand files were found. Upload .sdf, .smiles, .smi, or .csv files.", 400, result)
        csv_meta = _csv_upload_metadata(lig_dir, result)
        st = _load_state(ws)
        st["ligands_uploaded"] = True
        st["ligand_info"] = {"upload_mode": mode, "filename": upload_name, **result, **csv_meta}
        _save_state(ws, st)
        return _v1_ok(st["ligand_info"], warnings=result.get("warnings", []))

    @app.post("/api/v1/workspaces/<jobname>/ligands/curated")
    @login_required
    def api_v1_ligands_curated(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        library_key = request.form.get("library", "")
        curated = _curated_ligand_source(library_key)
        if not curated:
            return _v1_error("curated_library_missing", "Curated ligand library not found.", 404)

        lig_dir = ensure_subdir(ws, "Ligands")
        source_path = curated["path"]
        target_path = lig_dir / source_path.name
        shutil.copy2(source_path, target_path)

        result = {
            "accepted_files": [target_path.name],
            "accepted_count": 1,
            "ignored_files": [],
            "warnings": [],
            "is_csv": True,
            "filetypes": [".csv"],
            "ligands_root": str(lig_dir),
            "headers": curated["headers"],
            "suggested_smiles_col": curated["suggested_smiles_col"],
            "suggested_id_col": curated["suggested_id_col"],
            "curated_library": {
                "key": curated["key"],
                "label": curated["label"],
                "description": curated["description"],
            },
        }

        st = _load_state(ws)
        st["ligands_uploaded"] = True
        st["ligand_info"] = {"upload_mode": "curated", "filename": source_path.name, **result}
        _save_state(ws, st)
        return _v1_ok(st["ligand_info"], warnings=result.get("warnings", []))

    @app.get("/api/v1/workspaces/<jobname>/ligands")
    @login_required
    def api_v1_ligands_list(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        lig_dir = ws / "Ligands"
        files = [str(p.relative_to(lig_dir)) for p in sorted(lig_dir.rglob("*")) if p.is_file()] if lig_dir.exists() else []
        return _v1_ok({"ligands": files, "ligand_info": _load_state(ws).get("ligand_info") or {}})

    @app.post("/api/v1/workspaces/<jobname>/build")
    @login_required
    def api_v1_build(jobname: str):
        payload = _json_payload()
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        package_opts = payload.get("package") if isinstance(payload.get("package"), dict) else payload
        package_mode = normalize_package_mode(package_opts, current_app.config["DEFAULT_PACKAGE_MODE"], current_app.config["ENABLE_LSF_PACKAGE"])
        st = _load_state(ws)
        if not st.get("receptors"):
            return _v1_error("no_receptors", "Add receptors first.", 400)
        csv_map = _read_centers(ws, st)
        expected = [_receptor_pdbqt_name(r["rel"]) for r in st["receptors"]]
        if not all(n in csv_map for n in expected):
            return _v1_error("centers_incomplete", "Centers CSV is incomplete.", 400, {"expected_pdbqt": expected})
        if not st.get("ligands_uploaded"):
            return _v1_error("ligands_missing", "Upload ligands before building a package.", 400)
        ligand_info = st.get("ligand_info") or {}
        jobroot, warnings = assemble_job_tree(ws, ws / "Receptors", ws / "Ligands", package_mode=package_mode)
        def _int_value(*keys, default):
            for key in keys:
                if package_opts.get(key) not in (None, ""):
                    try:
                        return int(package_opts.get(key))
                    except Exception:
                        pass
            return int(default)
        def _str_value(*keys, default=""):
            for key in keys:
                if package_opts.get(key) not in (None, ""):
                    return str(package_opts.get(key)).strip()
            return default
        poses_conf = _int_value("confgen_poses", "poses_conf", "poses", default=64)
        poses_vina = _int_value("vina_poses", "poses_vina", default=20)
        inferred_lig_mode, inferred_lig_filetype, inferred_single_sdf = infer_ligand_workflow(ligand_info)
        lig_mode = _str_value("lig_mode", default=inferred_lig_mode)
        lig_filetype = _str_value("lig_filetype", "filetype", default=inferred_lig_filetype)
        csv_smiles_col = _str_value("csv_smiles_col", "smiles_col", default="")
        csv_id_col = _str_value("csv_id_col", "id_col", default="")
        single_sdf_rel = _str_value("single_sdf_rel", "single_sdf", default=inferred_single_sdf or "")
        if lig_mode == "1" and not csv_smiles_col:
            return _v1_error("csv_smiles_column_required", "Select a CSV SMILES column before building.", 400)
        profile = profile_for_mode(
            package_mode,
            {
                **package_opts,
                "lsf_email": _str_value("lsf_email", "email", default=""),
                "notify_begin": package_opts.get("notify_begin"),
                "notify_end": package_opts.get("notify_end"),
                "queue": _str_value("queue", default="general"),
                "project": _str_value("project", default=""),
                "workers": _int_value("workers", "ncores", default=16),
                "mem_per_core": _int_value("mem_per_core", "mem", default=2000),
                "confgen_walltime": _str_value("confgen_walltime", "walltime_confgen", default="48:00"),
                "vina_walltime": _str_value("vina_walltime", "walltime_vina", "walltime", default="96:00"),
                "conda_sh": _str_value("conda_sh", default=""),
                "conda_env": _str_value("conda_env", default=""),
                "vina_path": _str_value("vina_path", "vina_exe", default=""),
                "python_command": _str_value("python_command", default='$(command -v python3 || command -v python)'),
                "setup_commands": package_opts.get("setup_commands") or _str_value("env_line", default=current_app.config["DEFAULT_ENV_LINE"]),
            },
        )
        if package_mode in {"joey_lsf", "mainak_lsf", "custom_lsf"}:
            build_confgen_lsfs(
                jobroot, jobroot, profile=profile, poses=poses_conf,
                lig_mode=lig_mode,
                lig_filetype=lig_filetype,
                csv_smiles_col=csv_smiles_col,
                csv_id_col=csv_id_col,
                single_sdf_rel=(single_sdf_rel or None),
            )
            build_vina_lsfs(
                jobroot, jobroot, profile=profile, poses=poses_vina,
            )
        else:
            build_portable_runners(jobroot)
        rename_centers_with_tags(jobroot)
        z = zip_job_tree(jobroot)
        return _v1_ok({"zip": str(z), "download_url": url_for("download", path=str(z)), "package_mode": package_mode}, warnings=warnings)

    @app.get("/api/v1/workspaces/<jobname>/artifacts")
    @login_required
    def api_v1_artifacts(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        zips = sorted(ws.glob("*.zip"), key=lambda p: p.stat().st_mtime, reverse=True)
        return _v1_ok({"artifacts": [{"name": p.name, "path": str(p), "download_url": url_for("download", path=str(p))} for p in zips]})

    @app.get("/api/v1/workspaces/<jobname>/download")
    @login_required
    def api_v1_download(jobname: str):
        ws = _ws(jobname)
        if not ws.exists():
            return _v1_error("workspace_missing", f"Workspace {jobname} does not exist.", 404)
        requested = request.args.get("path", "")
        if requested:
            p = Path(requested)
            if not p.is_absolute():
                p = ws / requested
        else:
            zips = sorted(ws.glob("*.zip"), key=lambda item: item.stat().st_mtime, reverse=True)
            p = zips[0] if zips else Path()
        try:
            resolved = p.resolve()
            if not str(resolved).startswith(str(ws.resolve())) or not resolved.exists():
                return _v1_error("artifact_not_found", "No matching downloadable artifact was found.", 404)
            return send_file(resolved, as_attachment=True, download_name=resolved.name)
        except Exception:
            return _v1_error("artifact_not_found", "No matching downloadable artifact was found.", 404)

    @app.post("/api/v1/headless/package")
    @login_required
    def api_v1_headless_package():
        return _v1_error(
            "staged_workflow_required",
            "Use the staged /api/v1 workspace, receptor, center, ligand, prep, and build endpoints for this release.",
            501,
            {"docs": "/documentation", "reason": "Full one-call orchestration is deferred to avoid hiding prep/build failures."},
        )


    # ---------- DOWNLOAD ----------
    @app.get("/download")
    @login_required
    def download():
        path = request.args.get("path", "")
        p = Path(path)
        if not p.exists():
            return ("not found", 404)
        return send_file(p, as_attachment=True, download_name=p.name)

    return app


if __name__ == "__main__":
    app = create_app()
    app.run(host="0.0.0.0", port=5050, debug=True)
