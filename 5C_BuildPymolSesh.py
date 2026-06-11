#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5D_BuildPymolSesh_FINAL.py

Final PyMOL session/export builder for the updated docking workflow.

What this fixes
---------------
1) New ligand naming is collapsed correctly:
      DR8_p01_t03_c002  -> base ligand DR8
      DR7_p01_t04_c008  -> base ligand DR7
      4CI3_Y70_p01_t01_c007 -> base ligand 4CI3_Y70

   Therefore "Top 5 poses per ligand" means Top 5 per ORIGINAL ligand/base,
   not Top 5 for every generated protomer/tautomer/conformer variant.

2) Uses the OG PyMOL/ligand reconstruction behavior that kept ligands intact:
      PDBQT pose coordinates -> Open Babel mol -> RDKit MCS rigid fit of the
      generated/reference SDF chemistry onto the docked pose.

3) Still understands the newer folder/CSV provenance columns and generated
   SDF layout from the conformer-generation workflow.

Typical command
---------------
python 5D_BuildPymolSesh_FINAL.py \
  --csv Docking_Results_..._vina_docking_scores_sorted.csv \
  --mode per_ligand \
  --top 5 \
  --receptor-roots Receptors \
  --non-interactive

Audit selection only, without PyMOL/OpenBabel/RDKit reconstruction:
python 5D_BuildPymolSesh_FINAL.py \
  --csv Docking_Results_..._vina_docking_scores_sorted.csv \
  --mode per_ligand \
  --top 5 \
  --receptor-roots Receptors \
  --dry-run-selection \
  --non-interactive
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

# RDKit is required for ligand reconstruction, but not for --dry-run-selection.
try:
    from rdkit import Chem
    from rdkit.Chem import rdFMCS, rdMolAlign
    from rdkit.Chem.rdmolfiles import SDWriter
    HAVE_RDKIT = True
except Exception:
    Chem = None  # type: ignore
    rdFMCS = None  # type: ignore
    rdMolAlign = None  # type: ignore
    SDWriter = None  # type: ignore
    HAVE_RDKIT = False

# PyMOL is launched only when needed.
cmd = None


# ---------------------------------------------------------------------------
# Naming / parsing helpers
# ---------------------------------------------------------------------------

# New conformer-generation naming, anchored at the END so bases can contain _.
# Examples:
#   DR8_p01_t03_c002
#   4CI3_Y70_p01_t01_c007
#   Some_Base_Name_p02_t01_c064
NEW_VARIANT_RE = re.compile(
    r"^(?P<base>.+?)_p(?P<protomer>\d+)_t(?P<tautomer>\d+)(?:_c(?P<conformer>\d+))?$",
    re.IGNORECASE,
)

# Older naming seen in OG workflow.
OLD_POSE_SUFFIX_RE = re.compile(r"(.+?)_pose\d+$", re.IGNORECASE)
DOUBLE_UNDERSCORE_POSE_RE = re.compile(r"__pose\d+$", re.IGNORECASE)


def sanitize_filename(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")
    return cleaned or "item"


def sanitize_pymol_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", (name or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    if not cleaned or not cleaned[0].isalpha():
        cleaned = "obj_" + (cleaned or "item")
    return cleaned[:120]


def norm_base(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def strip_known_pose_suffix(text: str) -> str:
    value = (text or "").strip()
    value = DOUBLE_UNDERSCORE_POSE_RE.sub("", value)
    m = OLD_POSE_SUFFIX_RE.match(value)
    if m:
        return m.group(1)
    return value


def parse_ligand_name(name: str) -> Dict[str, str]:
    """Return stable base + variant metadata for both old and new naming."""
    raw = Path((name or "").strip()).stem
    raw = strip_known_pose_suffix(raw)

    m = NEW_VARIANT_RE.match(raw)
    if not m:
        return {
            "LigandVariant": raw,
            "LigandBase": raw,
            "ProtomerTag": "",
            "TautomerTag": "",
            "ConformerTag": "",
            "StateTag": "",
        }

    base = m.group("base")
    p = "p" + m.group("protomer")
    t = "t" + m.group("tautomer")
    c = "c" + m.group("conformer") if m.group("conformer") else ""
    tags = [p, t] + ([c] if c else [])
    return {
        "LigandVariant": raw,
        "LigandBase": base,
        "ProtomerTag": p,
        "TautomerTag": t,
        "ConformerTag": c,
        "StateTag": "_".join(tags),
    }


def infer_ligand_variant(row: Dict[str, Any]) -> str:
    """Find the most specific generated ligand variant from a CSV row."""
    for key in ("LigandVariant", "LigandDir"):
        value = (row.get(key) or "").strip()
        if value:
            return Path(value).name

    out = (row.get("OutFile") or "").strip()
    if out:
        parent = Path(out).parent.name
        if parent and parent not in {".", ""}:
            return parent

    value = (row.get("Ligand") or "").strip()
    if value:
        return Path(value).name

    return ""


def binding_column(row: Dict[str, Any]) -> Optional[str]:
    for key in ("Binding_Affinity", "Binding_Affinity_kcal_per_mol", "Affinity", "Score"):
        if key in row:
            return key
    return None


def safe_float(value: Any, default: float = 1e9) -> float:
    try:
        return float(value)
    except Exception:
        return default


# ---------------------------------------------------------------------------
# CLI / input helpers
# ---------------------------------------------------------------------------


def build_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build corrected PyMOL sessions/PDB/SDF/PDBQT exports from Vina score CSVs."
    )
    parser.add_argument("--csv", help="Docking score CSV to use")
    parser.add_argument(
        "--mode",
        choices=["per_ligand", "per_receptor"],
        help=(
            "per_ligand = top N poses per original/base ligand per receptor; "
            "per_receptor = top N original/base ligands per receptor, best pose only"
        ),
    )
    parser.add_argument("--top", type=int, help="Top N for the chosen mode")
    parser.add_argument("--receptor-roots", nargs="+", help="One or more folders containing receptors")
    parser.add_argument("--outdir", help="Output root directory; default=current working directory")
    parser.add_argument("--obabel-bin", help="Path to obabel executable; default=OBABEL_BIN or PATH")
    parser.add_argument(
        "--hydrogen-mode",
        choices=["none", "nonpolar", "all"],
        default="none",
        help="Hydrogen handling after OG-style fitting. Default none matches the clean OG behavior.",
    )
    parser.add_argument(
        "--allow-reference-fallback",
        action="store_true",
        help="If exact generated SDF is missing, allow base/original SDF fallback.",
    )
    parser.add_argument(
        "--include-previous-outputs-as-reference",
        action="store_true",
        help="Also search previous 11_AA_* output folders for reference SDFs. Default excludes them.",
    )
    parser.add_argument(
        "--min-mcs-fraction",
        type=float,
        default=0.80,
        help="Minimum heavy-atom MCS fraction accepted when fitting reference SDF to docked pose. Default 0.80.",
    )
    parser.add_argument(
        "--dry-run-selection",
        action="store_true",
        help="Only write the Top-N selection audit; do not run Open Babel, RDKit fitting, or PyMOL.",
    )
    parser.add_argument("--non-interactive", action="store_true", help="Require all needed arguments; do not prompt")
    return parser.parse_args()


def list_and_select(prompt: str, options: Sequence[str]) -> str:
    print("\n📁 " + prompt)
    for idx, opt in enumerate(options, 1):
        print(f"{idx}. {opt}")
    raw = input("🔢 Select: ").strip()
    if not raw.isdigit():
        raise SystemExit("❌ Invalid selection.")
    idx = int(raw) - 1
    if idx < 0 or idx >= len(options):
        raise SystemExit("❌ Selection out of range.")
    return options[idx]


def list_and_select_multiple(prompt: str, options: Sequence[str]) -> List[str]:
    print("\n📁 " + prompt)
    for idx, opt in enumerate(options, 1):
        print(f"{idx}. {opt}")
    print("\n🔢 Select one or more by number, comma-separated, or type 'all'.")
    raw = input("🔢 Select: ").strip()
    if not raw:
        raise SystemExit("❌ No selection provided.")
    if raw.lower() == "all":
        return list(options)

    chosen: List[str] = []
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
            chosen.append(value)
            seen.add(value)
    if not chosen:
        raise SystemExit("❌ No valid selection made.")
    return chosen


def prompt_for_inputs(args: argparse.Namespace, cwd: Path) -> Tuple[Path, str, int, List[Path], Path]:
    csv_path = Path(args.csv).expanduser() if args.csv else None
    mode = args.mode
    top_n = args.top
    receptor_roots = [Path(p).expanduser() for p in (args.receptor_roots or [])]
    out_root = Path(args.outdir).expanduser() if args.outdir else cwd

    if args.non_interactive:
        missing = []
        if not csv_path:
            missing.append("--csv")
        if not mode:
            missing.append("--mode")
        if not top_n:
            missing.append("--top")
        if not receptor_roots and not args.dry_run_selection:
            missing.append("--receptor-roots")
        if missing:
            raise SystemExit("❌ Non-interactive mode requires: " + ", ".join(missing))
    else:
        if not csv_path:
            csv_files = sorted([f.name for f in cwd.glob("*.csv")])
            if not csv_files:
                raise SystemExit("❌ No CSV files found in current directory.")
            csv_path = cwd / list_and_select("Available CSV files:", csv_files)

        if not mode:
            print("\n⚙️ Selection mode:")
            print("1. Top N poses per original/base ligand")
            print("2. Top N ligands total per receptor")
            mode_choice = input("🔢 Choose [1/2]: ").strip()
            if mode_choice == "1":
                mode = "per_ligand"
            elif mode_choice == "2":
                mode = "per_receptor"
            else:
                raise SystemExit("❌ Invalid choice.")

        if not top_n:
            top_n = int(input("🔝 How many top entries? (e.g., 5): ").strip())

        if not receptor_roots and not args.dry_run_selection:
            rec_dirs = sorted([d.name for d in cwd.iterdir() if d.is_dir() and "Receptor" in d.name])
            if not rec_dirs:
                raise SystemExit("❌ No receptor folders found in current directory.")
            selected = list_and_select_multiple("Receptor folders (choose one or more):", rec_dirs)
            receptor_roots = [cwd / item for item in selected]

    assert csv_path is not None
    assert mode is not None
    assert top_n is not None
    return csv_path.resolve(), mode, int(top_n), [p.resolve() for p in receptor_roots], out_root.resolve()


# ---------------------------------------------------------------------------
# CSV loading and Top-N selection
# ---------------------------------------------------------------------------


def detect_csv_style(rows: Sequence[Dict[str, Any]]) -> str:
    if not rows:
        return "simple"
    usable = 0
    for row in rows:
        out = (row.get("OutFile") or "").strip()
        rr = (row.get("ResultsRoot") or "").strip()
        rd = (row.get("ReceptorDir") or "").strip()
        ld = (row.get("LigandDir") or "").strip()
        if out or (rr and rd and ld):
            usable += 1
    return "provenance" if usable >= max(1, int(0.10 * len(rows))) else "simple"


def normalize_score_row(row: Dict[str, str]) -> Dict[str, Any]:
    cooked: Dict[str, Any] = dict(row)
    variant = infer_ligand_variant(cooked)
    parsed = parse_ligand_name(variant)

    # Preserve existing fields but set canonical new ones.
    cooked["LigandVariant"] = parsed["LigandVariant"] or variant
    cooked["LigandBase"] = parsed["LigandBase"] or strip_known_pose_suffix(cooked.get("Ligand", ""))
    cooked["ProtomerTag"] = parsed["ProtomerTag"]
    cooked["TautomerTag"] = parsed["TautomerTag"]
    cooked["ConformerTag"] = parsed["ConformerTag"]
    cooked["StateTag"] = parsed["StateTag"]
    if not cooked.get("Ligand"):
        cooked["Ligand"] = cooked["LigandBase"]
    return cooked


def load_csv_rows(csv_path: Path) -> Tuple[List[Dict[str, Any]], str]:
    if not csv_path.exists():
        raise SystemExit(f"❌ CSV not found: {csv_path}")
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [normalize_score_row(row) for row in reader]
    if not rows:
        raise SystemExit("❌ CSV is empty.")

    required = ["Receptor", "Pose"]
    missing = [key for key in required if key not in rows[0]]
    if missing:
        raise SystemExit(f"❌ CSV missing required columns: {missing}")

    bind_col = binding_column(rows[0])
    if bind_col:
        rows.sort(key=lambda row: safe_float(row.get(bind_col), 1e9))

    return rows, detect_csv_style(rows)


def select_rows(rows: Sequence[Dict[str, Any]], top_n: int, mode: str) -> Dict[str, List[Dict[str, Any]]]:
    """Correct Top-N logic using LigandBase, not LigandVariant."""
    grouped: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    if mode == "per_ligand":
        counts: Dict[Tuple[str, str], int] = defaultdict(int)
        for row in rows:
            receptor = (row.get("Receptor") or "").strip()
            base = (row.get("LigandBase") or "").strip()
            if not receptor or not base:
                continue
            key = (receptor, base)
            if counts[key] >= top_n:
                continue
            grouped[receptor].append(row)
            counts[key] += 1
        return grouped

    if mode == "per_receptor":
        seen_bases: Dict[str, set] = defaultdict(set)
        for row in rows:
            receptor = (row.get("Receptor") or "").strip()
            base = (row.get("LigandBase") or "").strip()
            if not receptor or not base:
                continue
            if base in seen_bases[receptor] or len(seen_bases[receptor]) >= top_n:
                continue
            grouped[receptor].append(row)
            seen_bases[receptor].add(base)
        return grouped

    raise SystemExit(f"❌ Unknown mode: {mode}")


def print_selection_summary(grouped: Dict[str, List[Dict[str, Any]]], mode: str, top_n: int) -> None:
    print("\n✅ Corrected Top-N selection summary")
    print(f"   Mode: {mode}")
    print(f"   Top N: {top_n}")
    for receptor, rows in grouped.items():
        counts = Counter(row.get("LigandBase", "") for row in rows)
        print(f"   - {receptor}: {len(rows)} selected rows across {len(counts)} base ligand(s)")
        for base, count in sorted(counts.items()):
            print(f"       {base}: {count}")


# ---------------------------------------------------------------------------
# Docked PDBQT and reference SDF resolution
# ---------------------------------------------------------------------------


def resolve_outfile_unified(row: Dict[str, Any], cwd: Path) -> Optional[Path]:
    out = (row.get("OutFile") or "").strip()
    if out:
        p = Path(out)
        if p.is_file():
            return p.resolve()
        if p.is_absolute() and not p.exists():
            # User may have moved the whole job folder; rebuild from the useful tail.
            for n_parts in (4, 5, 6):
                if len(p.parts) >= n_parts:
                    tail = Path(*p.parts[-n_parts:])
                    candidate = cwd / tail
                    if candidate.is_file():
                        return candidate.resolve()

    rr = (row.get("ResultsRoot") or "").strip()
    rd = (row.get("ReceptorDir") or row.get("Receptor") or "").strip()
    ld = (row.get("LigandDir") or row.get("LigandVariant") or "").strip()
    if rr and rd and ld:
        candidate = cwd / rr / rd / ld / "out.pdbqt"
        if candidate.is_file():
            return candidate.resolve()

    receptor = (row.get("Receptor") or "").strip()
    variant = (row.get("LigandVariant") or "").strip()
    ligand_field = (row.get("Ligand") or "").strip()
    names_to_try = [x for x in [variant, ligand_field] if x]

    for root in cwd.glob("Docking_Results*"):
        if not root.is_dir():
            continue
        for name in names_to_try:
            candidate = root / receptor / name / "out.pdbqt"
            if candidate.is_file():
                return candidate.resolve()

    # Last-resort exact parent-name search. This prevents accidentally selecting
    # every variant under a base ligand.
    for root in cwd.glob("Docking_Results*"):
        if not root.is_dir():
            continue
        for p in root.rglob("out.pdbqt"):
            if p.parent.name in names_to_try and (not receptor or receptor in str(p)):
                return p.resolve()

    return None


def _is_hidden_or_junk_path(path: Path) -> bool:
    return any(part.startswith(".") or part == "__MACOSX" or part.startswith("._") for part in path.parts)


def _is_previous_output_path(path: Path) -> bool:
    return any(part in {"11_AA_Ligands", "11_AA_PDBOUTPUT", "11_AA_PDBQT"} for part in path.parts)


class ReferenceResolver:
    """Find generated/reference SDFs for ligand variants."""

    def __init__(self, cwd: Path, include_previous_outputs_as_reference: bool = False):
        self.cwd = cwd.resolve()
        self.include_previous_outputs_as_reference = include_previous_outputs_as_reference
        self.tmp_roots = sorted([p.resolve() for p in self.cwd.glob("Ligands_TMP_SDF_*") if p.is_dir()])
        self.exact_index: Dict[str, List[Path]] = defaultdict(list)
        self.base_index: Dict[str, List[Path]] = defaultdict(list)
        self.manifest_rows: Dict[str, List[Tuple[Path, Dict[str, str]]]] = defaultdict(list)
        self._index()

    def _should_index_sdf(self, path: Path) -> bool:
        if _is_hidden_or_junk_path(path):
            return False
        if not self.include_previous_outputs_as_reference and _is_previous_output_path(path):
            return False
        return True

    def _candidate_roots(self) -> List[Path]:
        roots: List[Path] = []
        roots.extend(self.tmp_roots)
        for p in self.cwd.iterdir():
            if not p.is_dir():
                continue
            name = p.name
            if "Ligands" in name and "PDBQT" not in name and not name.startswith("."):
                roots.append(p.resolve())
        if self.include_previous_outputs_as_reference:
            for name in ("11_AA_Ligands", "11_AA_PDBOUTPUT"):
                p = self.cwd / name
                if p.is_dir():
                    roots.append(p.resolve())
        unique: List[Path] = []
        seen = set()
        for root in roots:
            if str(root) not in seen:
                unique.append(root)
                seen.add(str(root))
        return unique

    def _index(self) -> None:
        for root in self._candidate_roots():
            for sdf in root.rglob("*.sdf"):
                if not self._should_index_sdf(sdf):
                    continue
                sdf = sdf.resolve()
                parsed = parse_ligand_name(sdf.stem)
                self.exact_index[sdf.stem].append(sdf)
                self.base_index[parsed["LigandBase"]].append(sdf)

        for manifest in self.cwd.rglob("ligand_state_manifest.csv"):
            if _is_hidden_or_junk_path(manifest):
                continue
            if not self.include_previous_outputs_as_reference and _is_previous_output_path(manifest):
                continue
            try:
                with manifest.open(newline="", encoding="utf-8") as handle:
                    for row in csv.DictReader(handle):
                        variant = (
                            row.get("LigandVariant")
                            or row.get("Variant")
                            or row.get("Ligand")
                            or row.get("PDBQTStem")
                            or row.get("OutputStem")
                            or ""
                        ).strip()
                        if not variant:
                            # Try to infer from an SDF filename in the manifest row.
                            for key in ("TmpSDFFile", "SDFFile", "SDFPath", "SourceSDF", "GeneratedSDF"):
                                value = (row.get(key) or "").strip()
                                if value:
                                    variant = Path(value).stem
                                    break
                        if variant:
                            self.manifest_rows[variant].append((manifest.resolve(), row))
            except Exception:
                continue

    def _manifest_candidate_paths(self, manifest_path: Path, row: Dict[str, str]) -> Iterable[Path]:
        manifest_dir = manifest_path.parent
        keys = [
            "TmpSDFFile",
            "SDFFile",
            "SDFPath",
            "TmpSDFPath",
            "GeneratedSDF",
            "SourceSDF",
            "SourceInput",
        ]
        for key in keys:
            value = (row.get(key) or "").strip()
            if not value:
                continue
            raw = Path(value)
            if raw.is_absolute():
                yield raw
            else:
                yield manifest_dir / raw
                yield self.cwd / raw
                for root in self.tmp_roots:
                    yield root / raw
                    yield root / raw.name

    def _prefer(self, paths: Sequence[Path]) -> Optional[Path]:
        if not paths:
            return None
        paths = [p.resolve() for p in paths if p.exists() and self._should_index_sdf(p)]
        if not paths:
            return None
        # Prefer exact generated TMP SDFs, then shorter paths.
        return sorted(paths, key=lambda p: (0 if "Ligands_TMP_SDF_" in str(p) else 1, len(str(p))))[0]

    def resolve(self, row: Dict[str, Any], allow_fallback: bool = False) -> Dict[str, Any]:
        variant = (row.get("LigandVariant") or "").strip()
        base = (row.get("LigandBase") or parse_ligand_name(variant)["LigandBase"]).strip()
        parsed = parse_ligand_name(variant)
        checked: List[str] = []

        def pack(path: Optional[Path], method: str, status: str, warning: str = "", fallback: bool = False) -> Dict[str, Any]:
            return {
                "reference_sdf": path.resolve() if path and path.exists() else None,
                "reference_lookup_method": method,
                "reference_lookup_status": status,
                "reference_lookup_warning": warning,
                "is_exact_variant_match": bool(path and path.exists() and path.stem == variant and not fallback),
                "is_fallback_reference": fallback,
                "candidate_count": len(checked),
                "candidate_paths_checked": checked[:],
            }

        # 1) Manifest exact paths.
        for manifest_path, manifest_row in self.manifest_rows.get(variant, []):
            for candidate in self._manifest_candidate_paths(manifest_path, manifest_row):
                checked.append(str(candidate))
                if candidate.exists():
                    return pack(candidate, "manifest_exact", "exact_reference_sdf")

        # 2) Exact filename/stem search.
        exact = self._prefer(self.exact_index.get(variant, []))
        if exact:
            checked.extend([str(p) for p in self.exact_index.get(variant, [])])
            return pack(exact, "exact_filename_search", "exact_reference_sdf")

        # 3) Structured TMP inference: base/base_pXX_tYY/variant.sdf
        if parsed.get("StateTag"):
            state_no_conf = "_".join([x for x in [parsed.get("ProtomerTag"), parsed.get("TautomerTag")] if x])
            for root in self.tmp_roots:
                for candidate in [
                    root / base / f"{base}_{state_no_conf}" / f"{variant}.sdf",
                    root / base / f"{variant}.sdf",
                    root / f"{variant}.sdf",
                ]:
                    checked.append(str(candidate))
                    if candidate.exists():
                        return pack(candidate, "structured_tmp_inference", "exact_reference_sdf")

        # 4) Optional base/original fallback. This can still work with OG MCS fitting,
        # but it is intentionally marked as fallback in the audit.
        if allow_fallback:
            source_input = (row.get("SourceInput") or "").strip()
            if source_input:
                for candidate in [self.cwd / source_input, self.cwd / "Ligands" / Path(source_input).name]:
                    checked.append(str(candidate))
                    if candidate.exists():
                        return pack(
                            candidate,
                            "source_input_fallback",
                            "fallback_reference_sdf",
                            "Exact generated SDF missing; used SourceInput/base SDF fallback.",
                            fallback=True,
                        )
            fallback_path = self._prefer(self.base_index.get(base, []))
            if fallback_path:
                checked.extend([str(p) for p in self.base_index.get(base, [])])
                return pack(
                    fallback_path,
                    "base_filename_fallback",
                    "fallback_reference_sdf",
                    "Exact generated SDF missing; used base ligand SDF fallback.",
                    fallback=True,
                )

        return pack(None, "not_found", "missing_reference_sdf", "No exact generated SDF found for this ligand variant.")


# ---------------------------------------------------------------------------
# PDBQT pose extraction and ligand reconstruction
# ---------------------------------------------------------------------------


def obabel_path(cli_path: Optional[str]) -> str:
    found = cli_path or os.environ.get("OBABEL_BIN") or shutil.which("obabel")
    if not found:
        raise SystemExit("❌ Cannot find Open Babel `obabel`. Install it or set OBABEL_BIN/--obabel-bin.")
    return found


def extract_selected_pose_pdbqt(source_pdbqt: Path, pose_index: int, dest_pdbqt: Path) -> Tuple[Path, str]:
    """Write just the selected MODEL block from a Vina out.pdbqt."""
    text = source_pdbqt.read_text(errors="replace")
    lines = text.splitlines()
    models: List[Tuple[int, List[str]]] = []
    current_num: Optional[int] = None
    current_lines: List[str] = []

    for line in lines:
        if line.startswith("MODEL"):
            if current_lines:
                models.append((current_num or len(models) + 1, current_lines[:]))
            current_lines = [line]
            m = re.search(r"MODEL\s+(\d+)", line)
            current_num = int(m.group(1)) if m else len(models) + 1
            continue
        if current_lines:
            current_lines.append(line)
            if line.startswith("ENDMDL"):
                models.append((current_num or len(models) + 1, current_lines[:]))
                current_lines = []
                current_num = None

    if current_lines:
        models.append((current_num or len(models) + 1, current_lines[:]))

    dest_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    if not models:
        shutil.copy2(source_pdbqt, dest_pdbqt)
        return dest_pdbqt, "No MODEL blocks found; copied whole PDBQT."

    selected = next((block for number, block in models if number == pose_index), None)
    warning = ""
    if selected is None:
        number, selected = models[0]
        warning = f"Requested Vina pose {pose_index} not found; used MODEL {number}."

    dest_pdbqt.write_text("\n".join(selected) + "\n")
    return dest_pdbqt, warning


def load_reference_mol(reference_sdf: Optional[Path]) -> Optional[Any]:
    if not HAVE_RDKIT or reference_sdf is None:
        return None
    try:
        supplier = Chem.SDMolSupplier(str(reference_sdf), sanitize=False, removeHs=False)
        for mol in supplier:
            if mol is not None:
                return mol
    except Exception:
        return None
    return None


def pdbqt_to_mol(pdbqt_file: Path, obabel_bin: str) -> Any:
    if not HAVE_RDKIT:
        raise RuntimeError("RDKit is required for ligand reconstruction.")
    with tempfile.TemporaryDirectory(prefix="pdbqt_pose_") as tmpdir:
        mol_file = Path(tmpdir) / "pose.mol"
        subprocess.run(
            [obabel_bin, str(pdbqt_file), "-O", str(mol_file)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        mol = Chem.MolFromMolFile(str(mol_file), sanitize=False, removeHs=False)
        if mol is None:
            raise RuntimeError(f"Open Babel produced no molecule for {pdbqt_file}")
        return mol


def _noH_with_map(mol: Any) -> Tuple[Any, List[int]]:
    idx_map = []
    for idx, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() > 1:
            idx_map.append(idx)
    no_h = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    return no_h, idx_map


def _heavy_count(mol: Any) -> int:
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)


def rigid_fit_by_mcs(coords_full: Any, ref_full: Any, min_mcs_fraction: float = 0.80) -> Tuple[Optional[Any], str, int, str]:
    """
    OG-style reconstruction: rigidly align the entire reference SDF chemistry onto
    the docked coordinates using heavy-atom MCS.

    Returns: fitted_mol, status, map_size, warning
    """
    if coords_full is None or ref_full is None:
        return None, "missing_molecule", 0, "Missing docked or reference molecule."
    if coords_full.GetNumConformers() == 0 or ref_full.GetNumConformers() == 0:
        return None, "missing_conformer", 0, "Missing conformers; cannot align."

    coords_no_h, coords_map = _noH_with_map(coords_full)
    ref_no_h, ref_map = _noH_with_map(ref_full)
    ref_heavy = ref_no_h.GetNumAtoms()
    coords_heavy = coords_no_h.GetNumAtoms()
    min_heavy = min(ref_heavy, coords_heavy)

    if min_heavy <= 0:
        return None, "no_heavy_atoms", 0, "No heavy atoms found."

    mcs = rdFMCS.FindMCS([ref_no_h, coords_no_h], timeout=10)
    if mcs.canceled or not mcs.smartsString:
        return None, "mcs_failed", 0, "MCS alignment canceled or failed."

    patt = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = list(ref_no_h.GetSubstructMatch(patt))
    coords_match = list(coords_no_h.GetSubstructMatch(patt))

    if len(ref_match) != len(coords_match) or len(ref_match) == 0:
        return None, "mcs_match_mismatch", 0, "MCS atom-match mismatch."

    map_size = len(ref_match)
    frac = map_size / float(min_heavy)
    if frac < min_mcs_fraction:
        return (
            None,
            "mcs_too_small",
            map_size,
            f"MCS covered {map_size}/{min_heavy} heavy atoms ({frac:.2%}), below threshold {min_mcs_fraction:.2%}.",
        )

    atom_map = []
    for i_ref_no_h, i_coords_no_h in zip(ref_match, coords_match):
        ref_idx_full = ref_map[i_ref_no_h]
        coords_idx_full = coords_map[i_coords_no_h]
        atom_map.append((ref_idx_full, coords_idx_full))

    fitted = Chem.Mol(ref_full)
    try:
        rdMolAlign.AlignMol(fitted, coords_full, atomMap=atom_map)
    except Exception as exc:
        return None, "alignmol_failed", map_size, f"AlignMol failed: {exc}"

    # Match OG: remove hydrogens after alignment to avoid disconnected/weird H artifacts.
    fitted = Chem.RemoveHs(Chem.Mol(fitted), sanitize=False)
    try:
        Chem.SanitizeMol(
            fitted,
            sanitizeOps=(
                Chem.SanitizeFlags.SANITIZE_PROPERTIES
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
            ),
        )
    except Exception:
        pass

    warning = ""
    if ref_heavy != coords_heavy:
        warning = f"Reference heavy atoms ({ref_heavy}) != docked heavy atoms ({coords_heavy}); OG-style MCS fit used."
    elif map_size != ref_heavy:
        warning = f"MCS used {map_size}/{ref_heavy} heavy atoms; fitted reference kept intact."

    return fitted, "og_mcs_rigid_fit", map_size, warning


def apply_hydrogen_mode(mol: Any, mode: str) -> Any:
    """Default 'none' matches the OG output. Optional modes are convenience only."""
    if mode == "none":
        return Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    if mode == "all":
        return Chem.AddHs(Chem.Mol(mol), addCoords=True)
    if mode == "nonpolar":
        no_h = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
        with_h = Chem.AddHs(Chem.Mol(no_h), addCoords=True)
        edit = Chem.RWMol(with_h)
        remove = []
        for atom in edit.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbors = atom.GetNeighbors()
            if not neighbors or neighbors[0].GetAtomicNum() != 6:
                remove.append(atom.GetIdx())
        for idx in sorted(remove, reverse=True):
            edit.RemoveAtom(idx)
        return edit.GetMol()
    return mol


def write_sdf(path: Path, mol: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = SDWriter(str(path))
    writer.SetKekulize(False)
    writer.write(mol)
    writer.close()


def write_mol_file(path: Path, mol: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(Chem.MolToMolBlock(mol, kekulize=False))


# ---------------------------------------------------------------------------
# Receptor/PyMOL helpers
# ---------------------------------------------------------------------------


def find_receptor_file_in_root(receptors_root: Path, rec_name: str) -> Optional[Path]:
    rec_clean = rec_name.replace(".pdbqt", "").replace(".converted", "")
    stems = {rec_name, rec_clean}
    for root, _, files in os.walk(receptors_root):
        for filename in files:
            if filename.startswith("._"):
                continue
            p = Path(root) / filename
            if p.stem in stems or p.name.startswith(rec_clean) or p.name.startswith(rec_name):
                if p.suffix.lower() in {".pdb", ".pdbqt", ".mol2", ".cif", ".ent"} or str(p).endswith(".converted.pdbqt"):
                    return p.resolve()
    return None


def find_receptor_file(receptor_roots: Sequence[Path], rec_name: str) -> Optional[Path]:
    for root in receptor_roots:
        hit = find_receptor_file_in_root(root, rec_name)
        if hit is not None:
            return hit
    return None


def strip_conect_lines(pdb_path: Path) -> None:
    try:
        lines = [ln for ln in pdb_path.read_text(errors="replace").splitlines() if not ln.startswith("CONECT")]
        pdb_path.write_text("\n".join(lines) + "\n")
    except Exception as exc:
        print(f"⚠️ Could not strip CONECT from {pdb_path}: {exc}")


def launch_pymol() -> bool:
    global cmd
    try:
        import pymol as _pymol
        try:
            _pymol.pymol_argv = ["pymol", "-cq"]
            _pymol.finish_launching()
        except Exception:
            pass
        from pymol import cmd as pymol_cmd
        cmd = pymol_cmd
        cmd.reinitialize()
        try:
            cmd.feedback("disable", "all", "warnings")
            cmd.set("pdb_conect_all", 0)
            cmd.set("pdb_conect_nodup", 1)
        except Exception:
            pass
        return True
    except Exception as exc:
        print(f"⚠️ PyMOL is not available; PSE and complex PDB exports will be skipped. ({exc})")
        return False


def safe_group(group_name: str, member: str) -> None:
    if cmd is None:
        return
    try:
        cmd.group(group_name, member)
    except Exception:
        # Grouping is cosmetic; never let it break exports.
        pass


def load_receptor_once(pymol_enabled: bool, receptor_file: Path, receptor: str, loaded: Dict[str, str]) -> str:
    rec_obj = sanitize_pymol_name("obj_" + receptor)
    if not pymol_enabled:
        return rec_obj
    if receptor in loaded:
        return loaded[receptor]
    cmd.load(str(receptor_file), rec_obj)
    try:
        cmd.hide("everything", rec_obj)
        cmd.show("cartoon", rec_obj)
        cmd.color("gray70", rec_obj)
    except Exception:
        pass
    safe_group(sanitize_pymol_name("grp_" + receptor), rec_obj)
    loaded[receptor] = rec_obj
    return rec_obj


# ---------------------------------------------------------------------------
# Audit writers
# ---------------------------------------------------------------------------


AUDIT_COLUMNS = [
    "timestamp",
    "csv_style",
    "selection_mode",
    "receptor",
    "rank_for_receptor",
    "base_ligand",
    "ligand_variant",
    "protomer_tag",
    "tautomer_tag",
    "conformer_tag",
    "state_tag",
    "vina_pose",
    "binding_affinity",
    "dock_pdbqt",
    "selected_pose_pdbqt",
    "reference_sdf",
    "reference_lookup_method",
    "reference_lookup_status",
    "reference_lookup_warning",
    "is_exact_variant_match",
    "is_fallback_reference",
    "fit_method",
    "fit_status",
    "mcs_map_size",
    "fit_warning",
    "hydrogen_mode",
    "receptor_file",
    "saved_complex_pdb",
    "saved_corrected_sdf",
    "saved_raw_pdbqt",
    "pymol_object",
    "pymol_group",
    "skipped_reason",
]


def write_audit(path: Path, rows: Sequence[Dict[str, Any]], columns: Sequence[str] = AUDIT_COLUMNS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------


def main() -> None:
    args = build_args()
    cwd = Path(".").resolve()
    csv_path, mode, top_n, receptor_roots, out_root = prompt_for_inputs(args, cwd)

    out_dir_pdb = out_root / "11_AA_PDBOUTPUT"
    out_dir_sdf = out_root / "11_AA_Ligands"
    out_dir_pdbqt = out_root / "11_AA_PDBQT"
    out_dir_pdb.mkdir(parents=True, exist_ok=True)
    out_dir_sdf.mkdir(parents=True, exist_ok=True)
    out_dir_pdbqt.mkdir(parents=True, exist_ok=True)

    rows, csv_style = load_csv_rows(csv_path)
    grouped = select_rows(rows, top_n, mode)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    print(f"\n📄 Reading: {csv_path}")
    print(f"🔎 Detected CSV style: {csv_style}")
    print_selection_summary(grouped, mode, top_n)

    if args.dry_run_selection:
        audit_rows: List[Dict[str, Any]] = []
        for receptor, selected_rows in grouped.items():
            for rank, row in enumerate(selected_rows, 1):
                audit_rows.append(
                    {
                        "timestamp": timestamp,
                        "csv_style": csv_style,
                        "selection_mode": mode,
                        "receptor": receptor,
                        "rank_for_receptor": rank,
                        "base_ligand": row.get("LigandBase", ""),
                        "ligand_variant": row.get("LigandVariant", ""),
                        "protomer_tag": row.get("ProtomerTag", ""),
                        "tautomer_tag": row.get("TautomerTag", ""),
                        "conformer_tag": row.get("ConformerTag", ""),
                        "state_tag": row.get("StateTag", ""),
                        "vina_pose": row.get("Pose", ""),
                        "binding_affinity": row.get("Binding_Affinity", row.get("Binding_Affinity_kcal_per_mol", "")),
                        "dock_pdbqt": str(resolve_outfile_unified(row, cwd) or ""),
                    }
                )
        audit_csv = out_dir_pdb / f"selection_audit_{timestamp}.csv"
        write_audit(audit_csv, audit_rows)
        print("\n🎉 Dry-run selection complete.")
        print(f"🧾 Selection audit: {audit_csv}")
        return

    if not HAVE_RDKIT:
        raise SystemExit("❌ RDKit is required. Install it in the docking env before running this script.")

    obabel_bin = obabel_path(args.obabel_bin)
    resolver = ReferenceResolver(cwd, include_previous_outputs_as_reference=args.include_previous_outputs_as_reference)
    pymol_enabled = launch_pymol()
    loaded_receptors: Dict[str, str] = {}
    audit_rows = []

    print("\n🔬 Building corrected ligand exports and PyMOL session...")

    for receptor, selected_rows in grouped.items():
        rec_name = receptor.replace(".pdbqt", "").replace(".converted", "")
        receptor_file = find_receptor_file(receptor_roots, rec_name)
        if receptor_file is None:
            print(f"⚠️ Receptor file not found for {rec_name}; skipping receptor.")
            continue

        rec_obj = load_receptor_once(pymol_enabled, receptor_file, rec_name, loaded_receptors)
        receptor_group = sanitize_pymol_name("grp_" + rec_name)

        for rank_for_receptor, row in enumerate(selected_rows, 1):
            ligand_variant = row.get("LigandVariant", "")
            ligand_base = row.get("LigandBase", "")
            pose_index = int(safe_float(row.get("Pose") or 1, 1))
            binding = row.get("Binding_Affinity", row.get("Binding_Affinity_kcal_per_mol", ""))

            prefix = sanitize_filename(
                f"{rec_name}_Top{rank_for_receptor:02d}_{ligand_base}_{ligand_variant}_vinaPose{pose_index}"
            )
            selected_pose_pdbqt = out_dir_pdbqt / f"{prefix}_selected.pdbqt"
            corrected_sdf = out_dir_sdf / f"{prefix}_corrected_docked.sdf"
            complex_pdb = out_dir_pdb / f"{prefix}_complex.pdb"

            audit: Dict[str, Any] = {
                "timestamp": timestamp,
                "csv_style": csv_style,
                "selection_mode": mode,
                "receptor": rec_name,
                "rank_for_receptor": rank_for_receptor,
                "base_ligand": ligand_base,
                "ligand_variant": ligand_variant,
                "protomer_tag": row.get("ProtomerTag", ""),
                "tautomer_tag": row.get("TautomerTag", ""),
                "conformer_tag": row.get("ConformerTag", ""),
                "state_tag": row.get("StateTag", ""),
                "vina_pose": pose_index,
                "binding_affinity": binding,
                "receptor_file": str(receptor_file),
                "selected_pose_pdbqt": str(selected_pose_pdbqt),
                "saved_raw_pdbqt": str(selected_pose_pdbqt),
                "saved_corrected_sdf": "",
                "saved_complex_pdb": "",
                "pymol_object": "",
                "pymol_group": "",
                "hydrogen_mode": args.hydrogen_mode,
                "skipped_reason": "",
            }

            dock_pdbqt = resolve_outfile_unified(row, cwd)
            audit["dock_pdbqt"] = str(dock_pdbqt or "")
            if dock_pdbqt is None or not dock_pdbqt.exists():
                audit["skipped_reason"] = "missing_docked_pdbqt"
                print(f"⚠️ Missing docked PDBQT for {rec_name}/{ligand_variant}; skipping.")
                audit_rows.append(audit)
                continue

            try:
                _, pose_warning = extract_selected_pose_pdbqt(dock_pdbqt, pose_index, selected_pose_pdbqt)
            except Exception as exc:
                audit["skipped_reason"] = "pose_extraction_failed"
                audit["fit_warning"] = str(exc)
                print(f"⚠️ Could not extract pose for {ligand_variant}: {exc}")
                audit_rows.append(audit)
                continue

            resolution = resolver.resolve(row, allow_fallback=args.allow_reference_fallback)
            reference_sdf = resolution.get("reference_sdf")
            audit.update(
                {
                    "reference_sdf": str(reference_sdf or ""),
                    "reference_lookup_method": resolution.get("reference_lookup_method", ""),
                    "reference_lookup_status": resolution.get("reference_lookup_status", ""),
                    "reference_lookup_warning": resolution.get("reference_lookup_warning", ""),
                    "is_exact_variant_match": resolution.get("is_exact_variant_match", False),
                    "is_fallback_reference": resolution.get("is_fallback_reference", False),
                }
            )

            try:
                coords_mol = pdbqt_to_mol(selected_pose_pdbqt, obabel_bin)
            except Exception as exc:
                audit["skipped_reason"] = "pdbqt_to_mol_failed"
                audit["fit_warning"] = str(exc)
                print(f"⚠️ PDBQT conversion failed for {ligand_variant}: {exc}")
                audit_rows.append(audit)
                continue

            ref_mol = load_reference_mol(reference_sdf)
            fitted_mol = None
            fit_status = ""
            fit_method = ""
            mcs_size = 0
            fit_warning = pose_warning or ""

            if ref_mol is not None:
                fitted_mol, fit_status, mcs_size, warning = rigid_fit_by_mcs(
                    coords_mol, ref_mol, min_mcs_fraction=args.min_mcs_fraction
                )
                fit_method = "og_mcs_rigid_fit"
                if warning:
                    fit_warning = (fit_warning + " " + warning).strip()

            if fitted_mol is None:
                # Fallback to the direct coordinates from PDBQT conversion. This still
                # keeps the actual docked geometry rather than fabricating atoms.
                fitted_mol = Chem.RemoveHs(Chem.Mol(coords_mol), sanitize=False)
                if not fit_status:
                    fit_status = "coords_only_fallback"
                fit_method = fit_method or "coords_only_from_pdbqt"
                if ref_mol is None:
                    fit_warning = (fit_warning + " No reference SDF was available; used PDBQT coordinates only.").strip()

            final_mol = apply_hydrogen_mode(fitted_mol, args.hydrogen_mode)
            try:
                write_sdf(corrected_sdf, final_mol)
                audit["saved_corrected_sdf"] = str(corrected_sdf)
            except Exception as exc:
                audit["skipped_reason"] = "write_sdf_failed"
                audit["fit_warning"] = (fit_warning + " " + str(exc)).strip()
                print(f"⚠️ Could not write corrected SDF for {ligand_variant}: {exc}")
                audit_rows.append(audit)
                continue

            ligand_obj = sanitize_pymol_name(f"lig_{rec_name}_Top{rank_for_receptor:02d}_{ligand_variant}_pose{pose_index}")
            ligand_group = sanitize_pymol_name(f"grp_{rec_name}_{ligand_base}")

            if pymol_enabled:
                try:
                    cmd.load(str(corrected_sdf), ligand_obj)
                    cmd.hide("everything", ligand_obj)
                    cmd.show("sticks", ligand_obj)
                    safe_group(ligand_group, ligand_obj)
                    safe_group(receptor_group, ligand_group)
                    cmd.save(str(complex_pdb), selection=f"({rec_obj}) or ({ligand_obj})")
                    strip_conect_lines(complex_pdb)
                    audit["saved_complex_pdb"] = str(complex_pdb)
                    audit["pymol_object"] = ligand_obj
                    audit["pymol_group"] = ligand_group
                except Exception as exc:
                    audit["skipped_reason"] = "pymol_export_failed"
                    fit_warning = (fit_warning + " " + str(exc)).strip()
                    print(f"⚠️ PyMOL export failed for {ligand_variant}: {exc}")

            audit["fit_method"] = fit_method
            audit["fit_status"] = fit_status
            audit["mcs_map_size"] = mcs_size
            audit["fit_warning"] = fit_warning
            audit_rows.append(audit)

    session_path = out_dir_pdb / f"Top{top_n}_{mode}_Docking_Results_{timestamp}.pse"
    if pymol_enabled:
        try:
            cmd.save(str(session_path))
        except Exception as exc:
            print(f"⚠️ Could not save PyMOL session: {exc}")

    audit_csv = out_dir_pdb / f"top_hits_audit_{timestamp}.csv"
    write_audit(audit_csv, audit_rows)

    print("\n🎉 Done!")
    if pymol_enabled:
        print(f"🧪 PyMOL session: {session_path}")
    print(f"📦 PDB complexes: {out_dir_pdb}")
    print(f"📦 Corrected SDFs: {out_dir_sdf}")
    print(f"📦 Raw selected PDBQT: {out_dir_pdbqt}")
    print(f"🧾 Audit CSV: {audit_csv}")


if __name__ == "__main__":
    main()
