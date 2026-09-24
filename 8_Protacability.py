#!/usr/bin/env python3
from __future__ import annotations

"""
8_Protacability.py
==================
Consensus binding-mode and solvent-exposure analysis for AutoDock Vina ensembles.

Designed for the docking workflow that produces:
  Docking_Results_*/
      docking_parameters.json
      docking_jobs.csv
      <receptor>/<ligand_variant>/out.pdbqt

and the full Step-4 CSV written beside the Docking_Results_* directory:
  Docking_Results_*_<timestamp>_vina_docking_scores_sorted.csv

V1 capabilities
---------------
1. Read the COMPLETE Step-4 pose ensemble (never the Step-4C TOP-N subset).
2. Group poses by Receptor + LigandBase.
3. Conservatively split incompatible PDBQT atom mappings into MappingGroups.
4. Cluster poses using receptor-frame heavy-atom RMSD WITHOUT ligand superposition.
5. Choose a real medoid pose for every binding-mode cluster.
6. Report pose/job/conformer/top-pose support for every cluster.
7. Calculate heavy-atom relative SASA consensus for selected clusters using a
   Shrake-Rupley-style 1.4 A solvent probe.
8. Write auditable CSV/JSON output and representative medoid PDBQT files.

Important V1 interpretation
---------------------------
- Cluster support is docking-ensemble support, NOT thermodynamic probability.
- RMSD currently assumes identical PDBQT heavy-atom ordering inside a MappingGroup.
  Poses with incompatible signatures are never silently mixed.
- SASA is calculated from the atoms present in the docking PDBQT/receptor files.
  Relative SASA (bound/free for the same pose) is the primary exposure metric.
- Chemical derivatizability, persistent protein interactions, and linker-vector
  clearance are intentionally reserved for the next version.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import numpy as np
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"ERROR: NumPy is required for Step 8: {exc}")

try:
    from rdkit.ML.Cluster import Butina
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"ERROR: RDKit is required for Butina clustering: {exc}")

try:
    from Bio.PDB.SASA import ATOMIC_RADII, KDTree
except Exception as exc:  # pragma: no cover
    raise SystemExit(f"ERROR: Biopython Bio.PDB.SASA is required for SASA analysis: {exc}")


# -----------------------------------------------------------------------------
# Constants / atom typing
# -----------------------------------------------------------------------------

VINA_RESULT_RE = re.compile(r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)")
MODEL_RE = re.compile(r"^MODEL\s+(\d+)")

# Common AutoDock/Vina PDBQT atom types -> chemical element.
AD_TYPE_TO_ELEMENT = {
    "H": "H", "HD": "H", "HS": "H",
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    "MG": "MG", "MN": "MN", "ZN": "ZN", "CA": "CA", "FE": "FE",
    "CU": "CU", "NI": "NI", "CO": "CO", "NA+": "NA", "K": "K",
}

# Fallback radii for uncommon elements not present in Bio.PDB.SASA.ATOMIC_RADII.
FALLBACK_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "B": 1.92, "SI": 2.10,
}


# -----------------------------------------------------------------------------
# Data models
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class AtomRecord:
    serial: int
    name: str
    element: str
    atom_type: str
    coord: Tuple[float, float, float]


@dataclass
class ParsedModel:
    model_number: int
    atoms: List[AtomRecord]
    lines: List[str]
    affinity: Optional[float] = None


@dataclass
class PoseRecord:
    receptor: str
    ligand_base: str
    ligand_variant: str
    pose: int
    affinity: float
    outfile: Path
    protomer_tag: str = ""
    tautomer_tag: str = ""
    conformer_tag: str = ""
    state_tag: str = ""
    seed: str = ""
    receptor_file: Optional[Path] = None
    model: Optional[ParsedModel] = None
    mapping_group: str = ""
    cluster_id: int = 0
    rmsd_to_medoid: Optional[float] = None

    @property
    def heavy_atoms(self) -> List[AtomRecord]:
        if self.model is None:
            return []
        return [a for a in self.model.atoms if a.element.upper() != "H"]

    @property
    def heavy_coords(self) -> np.ndarray:
        return np.asarray([a.coord for a in self.heavy_atoms], dtype=np.float32)

    @property
    def job_key(self) -> str:
        return str(self.outfile.resolve())

    @property
    def conformer_key(self) -> str:
        pieces = [self.state_tag.strip(), self.conformer_tag.strip()]
        pieces = [p for p in pieces if p]
        if pieces:
            return "|".join(pieces)
        if self.ligand_variant:
            return self.ligand_variant
        return self.job_key


@dataclass
class ClusterResult:
    receptor: str
    ligand_base: str
    mapping_group: str
    cluster_id: int
    member_indices: Tuple[int, ...]
    medoid_index: int
    pose_count: int
    pose_support_ensemble_pct: float
    pose_support_mapping_pct: float
    job_count: int
    job_support_pct: float
    top_pose_job_count: int
    top_pose_job_support_pct: float
    conformer_count: int
    conformer_support_pct: float
    variant_count: int
    variant_support_pct: float
    best_affinity: float
    median_affinity: float
    mean_affinity: float
    mean_rmsd_to_medoid: float
    max_rmsd_to_medoid: float
    sasa_analyzed: bool = False
    sasa_pose_count: int = 0


@dataclass
class ReceptorContext:
    path: Path
    atoms: List[AtomRecord]
    coords: np.ndarray
    elements: List[str]
    expanded_radii: np.ndarray
    kdtree: object


# -----------------------------------------------------------------------------
# Small UI / path helpers
# -----------------------------------------------------------------------------

def find_result_dirs(base: Path) -> List[Path]:
    return [p.resolve() for p in sorted(base.glob("Docking_Results_*")) if p.is_dir()]


def choose_one(items: Sequence[Path], title: str) -> Path:
    print(f"\n{title}")
    for i, item in enumerate(items, 1):
        print(f" [{i}] {item.name}")
    while True:
        raw = input("Enter index: ").strip()
        try:
            idx = int(raw)
            if 1 <= idx <= len(items):
                return items[idx - 1]
        except Exception:
            pass
        print("Invalid selection. Try again.")


def latest_step4_csv(results_dir: Path) -> Optional[Path]:
    pattern = f"{results_dir.name}_*_vina_docking_scores_sorted.csv"
    candidates = [p for p in results_dir.parent.glob(pattern) if p.is_file()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: (p.stat().st_mtime, p.name))


def resolve_maybe_relative(path_text: str, root: Path) -> Path:
    p = Path(path_text).expanduser()
    if p.is_absolute():
        return p.resolve()
    return (root / p).resolve()


def safe_component(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")
    return text or "X"


def fmt_pct(value: float) -> float:
    return round(float(value), 4)


# -----------------------------------------------------------------------------
# PDB/PDBQT/MOL2 parsing
# -----------------------------------------------------------------------------

def normalize_element(element: str) -> str:
    e = re.sub(r"[^A-Za-z]", "", (element or "").strip()).upper()
    if not e:
        return "C"
    if len(e) >= 2 and e[:2] in {"CL", "BR", "MG", "MN", "ZN", "FE", "CU", "NI", "CO", "NA", "CA", "SI"}:
        return e[:2]
    return e[0]


def infer_element(atom_name: str, atom_type: str = "", pdb_element: str = "") -> str:
    if pdb_element.strip():
        return normalize_element(pdb_element)

    at = (atom_type or "").strip().upper()
    if at in AD_TYPE_TO_ELEMENT:
        return AD_TYPE_TO_ELEMENT[at]
    if at:
        letters = re.sub(r"[^A-Za-z]", "", at).upper()
        if letters in AD_TYPE_TO_ELEMENT:
            return AD_TYPE_TO_ELEMENT[letters]
        if letters.startswith("CL"):
            return "CL"
        if letters.startswith("BR"):
            return "BR"
        if letters:
            return normalize_element(letters)

    name = re.sub(r"^[0-9]+", "", atom_name or "").strip()
    return normalize_element(name)


def parse_pdb_like_atom(line: str) -> Optional[AtomRecord]:
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        serial = int(line[6:11].strip() or 0)
    except Exception:
        serial = 0
    name = line[12:16].strip() or f"A{serial}"

    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        toks = line.split()
        if len(toks) < 8:
            return None
        # Typical PDBQT token layout has XYZ around tokens 5-7; this fallback
        # scans for the first plausible run of three floats after the atom name.
        xyz = None
        for start in range(4, max(5, len(toks) - 2)):
            try:
                vals = tuple(float(toks[start + k]) for k in range(3))
            except Exception:
                continue
            xyz = vals
            break
        if xyz is None:
            return None
        x, y, z = xyz

    toks = line.split()
    atom_type = toks[-1] if toks else ""
    pdb_element = line[76:78].strip() if len(line) >= 78 else ""
    # In PDBQT, columns 77-78 may not be a true PDB element. Prefer AD atom type.
    if atom_type.upper() in AD_TYPE_TO_ELEMENT or re.fullmatch(r"[A-Za-z]{1,3}", atom_type or ""):
        element = infer_element(name, atom_type=atom_type)
    else:
        element = infer_element(name, pdb_element=pdb_element)

    return AtomRecord(serial, name, element, atom_type, (x, y, z))


def parse_pdbqt_models(path: Path) -> Dict[int, ParsedModel]:
    """Parse all Vina MODEL blocks. If MODEL records are absent, return model 1."""
    models: Dict[int, ParsedModel] = {}
    current_num: Optional[int] = None
    current_lines: List[str] = []
    current_atoms: List[AtomRecord] = []
    current_affinity: Optional[float] = None
    saw_model = False

    def flush() -> None:
        nonlocal current_num, current_lines, current_atoms, current_affinity
        if not current_lines and not current_atoms:
            return
        number = current_num if current_num is not None else (max(models.keys(), default=0) + 1)
        models[number] = ParsedModel(number, list(current_atoms), list(current_lines), current_affinity)
        current_num = None
        current_lines = []
        current_atoms = []
        current_affinity = None

    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                mm = MODEL_RE.match(line)
                if mm:
                    saw_model = True
                    flush()
                    current_num = int(mm.group(1))
                    current_lines = [line]
                    continue

                if saw_model and current_num is None:
                    continue

                current_lines.append(line)
                vm = VINA_RESULT_RE.match(line)
                if vm:
                    try:
                        current_affinity = float(vm.group(1))
                    except Exception:
                        pass
                atom = parse_pdb_like_atom(line)
                if atom is not None:
                    current_atoms.append(atom)
                if line.startswith("ENDMDL"):
                    flush()

        flush()
    except OSError:
        return {}

    if not saw_model and models:
        # flush() will have created model 1 from the whole file.
        return models
    return models


def parse_receptor_atoms(path: Path) -> List[AtomRecord]:
    suffix = path.suffix.lower()
    atoms: List[AtomRecord] = []

    if suffix == ".mol2":
        in_atoms = False
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atoms = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atoms:
                    break
                if not in_atoms or not line.strip():
                    continue
                toks = line.split()
                if len(toks) < 6:
                    continue
                try:
                    serial = int(toks[0])
                    name = toks[1]
                    x, y, z = map(float, toks[2:5])
                    atom_type = toks[5]
                except Exception:
                    continue
                element = infer_element(name, atom_type=atom_type.split(".")[0])
                atoms.append(AtomRecord(serial, name, element, atom_type, (x, y, z)))
        return atoms

    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            atom = parse_pdb_like_atom(line)
            if atom is not None:
                atoms.append(atom)
    return atoms


# -----------------------------------------------------------------------------
# Manifest / score loading
# -----------------------------------------------------------------------------

def load_job_manifest(results_dir: Path) -> Dict[str, Dict[str, str]]:
    path = results_dir / "docking_jobs.csv"
    if not path.is_file():
        return {}

    mapping: Dict[str, Dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out_text = (row.get("OutputPDBQT") or "").strip()
            if not out_text:
                continue
            out_path = resolve_maybe_relative(out_text, results_dir)
            mapping[str(out_path)] = row
    return mapping


def config_receptor_path(outfile: Path) -> Optional[Path]:
    cfg = outfile.with_name("config.txt")
    if not cfg.is_file():
        return None
    try:
        for line in cfg.read_text(errors="ignore").splitlines():
            if line.strip().lower().startswith("receptor") and "=" in line:
                value = line.split("=", 1)[1].strip().strip('"').strip("'")
                p = Path(value).expanduser()
                if not p.is_absolute():
                    p = (cfg.parent / p).resolve()
                return p.resolve()
    except Exception:
        pass
    return None


def load_pose_records(
    scores_csv: Path,
    results_dir: Path,
    job_manifest: Dict[str, Dict[str, str]],
    requested_ligands: Optional[set[str]] = None,
) -> Tuple[List[PoseRecord], List[str]]:
    warnings: List[str] = []
    rows: List[Dict[str, str]] = []
    with open(scores_csv, newline="", encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh)
        required = {"Receptor", "LigandBase", "LigandVariant", "Pose", "Binding_Affinity", "OutFile"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"ERROR: Step-4 CSV is missing columns {sorted(missing)}. Found: {reader.fieldnames}"
            )
        rows = list(reader)

    # Parse every physical out.pdbqt only once.
    model_cache: Dict[str, Dict[int, ParsedModel]] = {}
    records: List[PoseRecord] = []
    missing_outfiles = Counter()
    missing_models = Counter()

    for row in rows:
        ligand_base = (row.get("LigandBase") or row.get("Ligand") or "").strip()
        if requested_ligands and ligand_base not in requested_ligands:
            continue

        out_text = (row.get("OutFile") or "").strip()
        if not out_text:
            continue
        outfile = Path(out_text).expanduser()
        if not outfile.is_absolute():
            outfile = (results_dir.parent / outfile).resolve()
        else:
            outfile = outfile.resolve()

        if not outfile.is_file():
            missing_outfiles[str(outfile)] += 1
            continue

        key = str(outfile)
        if key not in model_cache:
            model_cache[key] = parse_pdbqt_models(outfile)

        try:
            pose_num = int(float(row.get("Pose") or 0))
        except Exception:
            continue
        model = model_cache[key].get(pose_num)
        if model is None:
            missing_models[key] += 1
            continue
        if not model.atoms:
            missing_models[key] += 1
            continue

        try:
            affinity = float(row.get("Binding_Affinity") or model.affinity)
        except Exception:
            if model.affinity is None:
                continue
            affinity = float(model.affinity)

        manifest_row = job_manifest.get(key, {})
        receptor_file: Optional[Path] = None
        receptor_text = (manifest_row.get("ReceptorFile") or "").strip()
        if receptor_text:
            receptor_file = Path(receptor_text).expanduser().resolve()
        if receptor_file is None or not receptor_file.is_file():
            receptor_file = config_receptor_path(outfile)

        records.append(
            PoseRecord(
                receptor=(row.get("Receptor") or "").strip(),
                ligand_base=ligand_base,
                ligand_variant=(row.get("LigandVariant") or row.get("Ligand") or "").strip(),
                pose=pose_num,
                affinity=affinity,
                outfile=outfile,
                protomer_tag=(row.get("ProtomerTag") or "").strip(),
                tautomer_tag=(row.get("TautomerTag") or "").strip(),
                conformer_tag=(row.get("ConformerTag") or "").strip(),
                state_tag=(row.get("StateTag") or "").strip(),
                seed=(manifest_row.get("Seed") or "").strip(),
                receptor_file=receptor_file,
                model=model,
            )
        )

    if missing_outfiles:
        warnings.append(f"{sum(missing_outfiles.values())} Step-4 rows referenced missing out.pdbqt files.")
    if missing_models:
        warnings.append(f"{sum(missing_models.values())} Step-4 rows referenced MODEL blocks that could not be parsed.")

    return records, warnings


# -----------------------------------------------------------------------------
# Mapping signatures / RMSD / clustering
# -----------------------------------------------------------------------------

def atom_mapping_signature(pose: PoseRecord) -> Tuple[Tuple[str, str, str], ...]:
    # Conservative by design: if atom identity/order differs, do not silently mix.
    sig = []
    for atom in pose.heavy_atoms:
        sig.append((atom.name.upper(), atom.element.upper(), atom.atom_type.upper()))
    return tuple(sig)


def signature_hash(sig: Tuple[Tuple[str, str, str], ...]) -> str:
    payload = json.dumps(sig, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()[:10]


def condensed_rmsd(coords: np.ndarray) -> np.ndarray:
    """
    Lower-triangle receptor-frame RMSD vector in Butina's expected order:
      for i in range(n):
          for j in range(i):
              d(i,j)
    No fitting/superposition is performed.
    """
    n = int(coords.shape[0])
    total = n * (n - 1) // 2
    dists = np.empty(total, dtype=np.float32)
    offset = 0
    for i in range(1, n):
        diff = coords[:i] - coords[i]
        vals = np.sqrt(np.mean(np.sum(diff * diff, axis=2), axis=1, dtype=np.float64))
        dists[offset: offset + i] = vals.astype(np.float32, copy=False)
        offset += i
    return dists


def condensed_get(dists: np.ndarray, i: int, j: int) -> float:
    if i == j:
        return 0.0
    if i < j:
        i, j = j, i
    return float(dists[i * (i - 1) // 2 + j])


def medoid_for_members(dists: np.ndarray, members: Sequence[int]) -> Tuple[int, float, float]:
    if len(members) == 1:
        return int(members[0]), 0.0, 0.0
    best_idx = int(members[0])
    best_mean = float("inf")
    best_max = float("inf")
    for i in members:
        vals = [condensed_get(dists, int(i), int(j)) for j in members if int(j) != int(i)]
        mean_d = float(statistics.fmean(vals)) if vals else 0.0
        max_d = max(vals) if vals else 0.0
        if mean_d < best_mean - 1e-9 or (abs(mean_d - best_mean) <= 1e-9 and int(i) < best_idx):
            best_idx = int(i)
            best_mean = mean_d
            best_max = max_d
    return best_idx, best_mean, best_max


def pct(num: int, den: int) -> float:
    return (100.0 * num / den) if den else 0.0


def cluster_mapping_group(
    poses: List[PoseRecord],
    cluster_rmsd: float,
    ensemble_pose_count: int,
) -> Tuple[List[ClusterResult], np.ndarray]:
    if not poses:
        return [], np.zeros(0, dtype=np.float32)

    coords = np.stack([p.heavy_coords for p in poses], axis=0)
    n = len(poses)
    if n == 1:
        raw_clusters = ((0,),)
        dists = np.zeros(0, dtype=np.float32)
    else:
        dists = condensed_rmsd(coords)
        raw_clusters = Butina.ClusterData(
            dists,
            nPts=n,
            distThresh=float(cluster_rmsd),
            isDistData=True,
            reordering=True,
        )

    # Present clusters by occupancy, deterministically.
    clusters_sorted = sorted(
        [tuple(int(x) for x in c) for c in raw_clusters],
        key=lambda c: (-len(c), min(c)),
    )

    all_jobs = {p.job_key for p in poses}
    all_confs = {p.conformer_key for p in poses}
    all_variants = {p.ligand_variant for p in poses}

    # Best-scoring pose for each physical Vina job.
    top_pose_by_job: Dict[str, int] = {}
    for idx, pose in enumerate(poses):
        current = top_pose_by_job.get(pose.job_key)
        if current is None or pose.affinity < poses[current].affinity:
            top_pose_by_job[pose.job_key] = idx

    results: List[ClusterResult] = []
    for cid, members in enumerate(clusters_sorted, 1):
        medoid, mean_rmsd, max_rmsd = medoid_for_members(dists, members)
        member_set = set(members)
        member_poses = [poses[i] for i in members]
        member_jobs = {p.job_key for p in member_poses}
        member_confs = {p.conformer_key for p in member_poses}
        member_variants = {p.ligand_variant for p in member_poses}
        top_job_count = sum(1 for idx in top_pose_by_job.values() if idx in member_set)
        affinities = [p.affinity for p in member_poses]

        for i in members:
            poses[i].cluster_id = cid
            poses[i].rmsd_to_medoid = condensed_get(dists, i, medoid)

        results.append(
            ClusterResult(
                receptor=poses[0].receptor,
                ligand_base=poses[0].ligand_base,
                mapping_group=poses[0].mapping_group,
                cluster_id=cid,
                member_indices=members,
                medoid_index=medoid,
                pose_count=len(members),
                pose_support_ensemble_pct=pct(len(members), ensemble_pose_count),
                pose_support_mapping_pct=pct(len(members), n),
                job_count=len(member_jobs),
                job_support_pct=pct(len(member_jobs), len(all_jobs)),
                top_pose_job_count=top_job_count,
                top_pose_job_support_pct=pct(top_job_count, len(all_jobs)),
                conformer_count=len(member_confs),
                conformer_support_pct=pct(len(member_confs), len(all_confs)),
                variant_count=len(member_variants),
                variant_support_pct=pct(len(member_variants), len(all_variants)),
                best_affinity=min(affinities),
                median_affinity=float(statistics.median(affinities)),
                mean_affinity=float(statistics.fmean(affinities)),
                mean_rmsd_to_medoid=mean_rmsd,
                max_rmsd_to_medoid=max_rmsd,
            )
        )

    return results, dists


# -----------------------------------------------------------------------------
# Fast ligand-only Shrake-Rupley SASA against receptor environment
# -----------------------------------------------------------------------------

def element_radius(element: str) -> float:
    e = normalize_element(element)
    if e in ATOMIC_RADII:
        return float(ATOMIC_RADII[e])
    # Biopython keys are generally uppercase; fall back to conventional values.
    return float(FALLBACK_RADII.get(e, 1.70))


def golden_sphere(n_points: int) -> np.ndarray:
    n = max(20, int(n_points))
    dl = math.pi * (3.0 - math.sqrt(5.0))
    dz = 2.0 / n
    longitude = 0.0
    z = 1.0 - dz / 2.0
    coords = np.zeros((n, 3), dtype=np.float32)
    for k in range(n):
        r = math.sqrt(max(0.0, 1.0 - z * z))
        coords[k, 0] = math.cos(longitude) * r
        coords[k, 1] = math.sin(longitude) * r
        coords[k, 2] = z
        z -= dz
        longitude += dl
    return coords


def make_receptor_context(path: Path, probe_radius: float) -> ReceptorContext:
    atoms = parse_receptor_atoms(path)
    if not atoms:
        raise RuntimeError(f"No receptor atoms parsed from {path}")
    coords = np.asarray([a.coord for a in atoms], dtype=np.float64)
    elements = [a.element for a in atoms]
    expanded = np.asarray([element_radius(e) + probe_radius for e in elements], dtype=np.float64)
    return ReceptorContext(
        path=path,
        atoms=atoms,
        coords=coords,
        elements=elements,
        expanded_radii=expanded,
        kdtree=KDTree(coords, 10),
    )


def ligand_relative_sasa(
    ligand_atoms: Sequence[AtomRecord],
    receptor: ReceptorContext,
    probe_radius: float,
    sphere: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (free_sasa, bound_sasa, relative_percent) for every ligand atom.

    We calculate SASA only for ligand atoms. Ligand atoms occlude one another in
    both states; receptor atoms additionally occlude the bound state. This is
    equivalent to the relevant ligand portion of a Shrake-Rupley complex SASA,
    without wasting time calculating receptor SASA values that Step 8 does not use.
    """
    lig_coords = np.asarray([a.coord for a in ligand_atoms], dtype=np.float64)
    lig_elements = [a.element for a in ligand_atoms]
    lig_radii = np.asarray([element_radius(e) + probe_radius for e in lig_elements], dtype=np.float64)
    lig_kdt = KDTree(lig_coords, 10)

    n_points = int(sphere.shape[0])
    max_expanded = max(
        float(np.max(lig_radii)) if len(lig_radii) else 0.0,
        float(np.max(receptor.expanded_radii)) if len(receptor.expanded_radii) else 0.0,
    )
    neighbor_cutoff = 2.0 * max_expanded + 1e-6

    free = np.zeros(len(ligand_atoms), dtype=np.float64)
    bound = np.zeros(len(ligand_atoms), dtype=np.float64)
    point_ids = set(range(n_points))

    for i in range(len(ligand_atoms)):
        r_i = lig_radii[i]
        surface = sphere.astype(np.float64, copy=False) * r_i + lig_coords[i]
        sphere_kdt = KDTree(surface, 10)
        available_free = point_ids.copy()

        # Occlusion by the rest of the ligand.
        for neighbor in lig_kdt.search(lig_coords[i], neighbor_cutoff):
            j = int(neighbor.index)
            if j == i:
                continue
            if float(neighbor.radius) < (r_i + lig_radii[j]):
                available_free -= {
                    int(pt.index) for pt in sphere_kdt.search(lig_coords[j], lig_radii[j])
                }
                if not available_free:
                    break

        available_bound = available_free.copy()
        if available_bound:
            # Additional occlusion by receptor atoms.
            for neighbor in receptor.kdtree.search(lig_coords[i], neighbor_cutoff):
                j = int(neighbor.index)
                r_j = receptor.expanded_radii[j]
                if float(neighbor.radius) < (r_i + r_j):
                    available_bound -= {
                        int(pt.index) for pt in sphere_kdt.search(receptor.coords[j], r_j)
                    }
                    if not available_bound:
                        break

        area_factor = r_i * r_i * (4.0 * math.pi / n_points)
        free[i] = len(available_free) * area_factor
        bound[i] = len(available_bound) * area_factor

    relative = np.zeros(len(ligand_atoms), dtype=np.float64)
    mask = free > 1e-12
    relative[mask] = 100.0 * bound[mask] / free[mask]
    relative = np.clip(relative, 0.0, 100.0)
    return free, bound, relative


# -----------------------------------------------------------------------------
# Output writers
# -----------------------------------------------------------------------------

def extract_model_to_file(outfile: Path, pose_num: int, dest: Path) -> None:
    models = parse_pdbqt_models(outfile)
    model = models.get(int(pose_num))
    if model is None:
        raise RuntimeError(f"MODEL {pose_num} not found in {outfile}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(model.lines) + "\n", encoding="utf-8")


def write_pose_assignments(path: Path, poses: Iterable[PoseRecord]) -> None:
    fields = [
        "Receptor", "LigandBase", "LigandVariant", "MappingGroup", "ClusterID",
        "Pose", "Binding_Affinity", "RMSD_to_Medoid_A", "OutFile", "Seed",
        "ProtomerTag", "TautomerTag", "ConformerTag", "StateTag", "ConformerUnit",
        "ReceptorFile",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for p in poses:
            writer.writerow({
                "Receptor": p.receptor,
                "LigandBase": p.ligand_base,
                "LigandVariant": p.ligand_variant,
                "MappingGroup": p.mapping_group,
                "ClusterID": p.cluster_id,
                "Pose": p.pose,
                "Binding_Affinity": p.affinity,
                "RMSD_to_Medoid_A": "" if p.rmsd_to_medoid is None else round(p.rmsd_to_medoid, 6),
                "OutFile": str(p.outfile),
                "Seed": p.seed,
                "ProtomerTag": p.protomer_tag,
                "TautomerTag": p.tautomer_tag,
                "ConformerTag": p.conformer_tag,
                "StateTag": p.state_tag,
                "ConformerUnit": p.conformer_key,
                "ReceptorFile": str(p.receptor_file) if p.receptor_file else "",
            })


def write_binding_modes(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_atom_consensus(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = [
        "Receptor", "LigandBase", "MappingGroup", "ClusterID",
        "ClusterPoseSupportPct_Ensemble", "ClusterPoseSupportPct_MappingGroup",
        "ClusterJobSupportPct", "ClusterConformerSupportPct",
        "AtomIndex", "AtomName", "Element", "PDBQTType",
        "SASA_PoseCount", "MedianFreeSASA_A2", "MedianBoundSASA_A2",
        "MeanRelativeSASA_Pct", "MedianRelativeSASA_Pct",
        "RelativeSASA_Q1_Pct", "RelativeSASA_Q3_Pct",
        "ExposedAt25Pct_Poses", "ExposedAt25Pct_Pct",
        "ExposedAt50Pct_Poses", "ExposedAt50Pct_Pct",
        "ExposedAt75Pct_Poses", "ExposedAt75Pct_Pct",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


# -----------------------------------------------------------------------------
# Main analysis
# -----------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Consensus binding-mode + solvent-exposure analysis for Vina docking ensembles."
    )
    ap.add_argument("--dir", help="Docking_Results_* directory. Omit for interactive selection.")
    ap.add_argument("--scores-csv", help="Full Step-4 per-directory scores CSV. Auto-detected if omitted.")
    ap.add_argument("--ligand-base", action="append", default=[], help="Analyze only this LigandBase; repeatable.")
    ap.add_argument("--cluster-rmsd", type=float, default=2.0, help="Receptor-frame heavy-atom RMSD cutoff in A (default 2.0).")
    ap.add_argument("--probe-radius", type=float, default=1.4, help="SASA solvent probe radius in A (default 1.4).")
    ap.add_argument("--sasa-points", type=int, default=100, help="Shrake-Rupley sphere points per atom (default 100).")
    ap.add_argument(
        "--sasa-top-clusters", type=int, default=3,
        help="SASA-analyze top N clusters per MappingGroup; 0=all (default 3).",
    )
    ap.add_argument("--skip-sasa", action="store_true", help="Run clustering/medoids only; skip SASA for a fast validation run.")
    ap.add_argument("--output-dir", help="Output directory. Default: PROTACability_Results_<timestamp>.")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    start = datetime.now()
    cwd = Path.cwd().resolve()

    print("\n=== Step 8: Consensus PROTACability Analysis ===")

    if args.dir:
        results_dir = Path(args.dir).expanduser().resolve()
        if not results_dir.is_dir():
            raise SystemExit(f"ERROR: invalid --dir: {results_dir}")
    else:
        candidates = find_result_dirs(cwd)
        if not candidates:
            raise SystemExit("ERROR: no Docking_Results_* directories found in the current directory.")
        results_dir = choose_one(candidates, "Select Docking_Results_* directory:")

    if args.scores_csv:
        scores_csv = Path(args.scores_csv).expanduser().resolve()
    else:
        scores_csv = latest_step4_csv(results_dir)
        if scores_csv is None:
            raise SystemExit(
                "ERROR: no full Step-4 per-directory CSV was found. Run 4_ParseScores.py first, "
                "or provide --scores-csv. Do NOT use a Step-4C TOP-N CSV for consensus percentages."
            )

    if not scores_csv.is_file():
        raise SystemExit(f"ERROR: score CSV does not exist: {scores_csv}")

    if args.cluster_rmsd <= 0:
        raise SystemExit("ERROR: --cluster-rmsd must be > 0")
    if args.probe_radius <= 0:
        raise SystemExit("ERROR: --probe-radius must be > 0")
    if args.sasa_points < 20:
        raise SystemExit("ERROR: --sasa-points must be >= 20")

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    if args.output_dir:
        output_dir = Path(args.output_dir).expanduser().resolve()
    else:
        output_dir = (results_dir.parent / f"PROTACability_Results_{safe_component(results_dir.name)}_{timestamp}").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reps_dir = output_dir / "Representatives"
    reps_dir.mkdir(parents=True, exist_ok=True)

    provenance_path = results_dir / "docking_parameters.json"
    job_manifest = load_job_manifest(results_dir)
    requested = set(args.ligand_base) if args.ligand_base else None

    print(f"Docking results : {results_dir}")
    print(f"Step-4 scores   : {scores_csv}")
    print(f"Dock provenance : {provenance_path if provenance_path.is_file() else '(not found)'}")
    print(f"Job manifest    : {results_dir / 'docking_jobs.csv' if job_manifest else '(not found)'}")
    print(f"Cluster RMSD    : {args.cluster_rmsd:.3f} A (receptor-frame; NO fitting)")
    if args.skip_sasa:
        print("SASA            : SKIPPED (--skip-sasa)")
    else:
        label = "ALL clusters" if args.sasa_top_clusters == 0 else f"top {args.sasa_top_clusters} cluster(s) per MappingGroup"
        print(f"SASA            : {label}; probe={args.probe_radius:.2f} A; points={args.sasa_points}")

    print("\nReading complete Step-4 pose ensemble...")
    poses, warnings = load_pose_records(scores_csv, results_dir, job_manifest, requested_ligands=requested)
    if not poses:
        raise SystemExit("ERROR: no pose records could be loaded from the selected Step-4 CSV.")

    for w in warnings:
        print(f"WARNING: {w}")

    print(f"Loaded poses    : {len(poses):,}")
    print(f"Ligand ensembles: {len({(p.receptor, p.ligand_base) for p in poses}):,}")

    # Group by receptor + parent ligand.
    ensembles: Dict[Tuple[str, str], List[PoseRecord]] = defaultdict(list)
    for p in poses:
        ensembles[(p.receptor, p.ligand_base)].append(p)

    all_cluster_results: List[ClusterResult] = []
    binding_rows: List[Dict[str, object]] = []
    atom_rows: List[Dict[str, object]] = []
    analysis_warnings: List[str] = list(warnings)
    receptor_cache: Dict[str, ReceptorContext] = {}
    sphere = golden_sphere(args.sasa_points)

    ordered_ensembles = sorted(ensembles.items(), key=lambda kv: (kv[0][0].lower(), kv[0][1].lower()))

    for ensemble_num, ((receptor, ligand_base), ensemble_poses) in enumerate(ordered_ensembles, 1):
        print(
            f"\n[{ensemble_num}/{len(ordered_ensembles)}] {receptor} :: {ligand_base} "
            f"({len(ensemble_poses):,} poses)"
        )

        signature_groups: Dict[Tuple[Tuple[str, str, str], ...], List[PoseRecord]] = defaultdict(list)
        for p in ensemble_poses:
            sig = atom_mapping_signature(p)
            if sig:
                signature_groups[sig].append(p)

        sorted_sig_groups = sorted(signature_groups.items(), key=lambda kv: (-len(kv[1]), signature_hash(kv[0])))
        if len(sorted_sig_groups) > 1:
            msg = (
                f"{receptor}/{ligand_base}: {len(sorted_sig_groups)} incompatible PDBQT heavy-atom mappings detected; "
                "analyzed as separate MappingGroups rather than silently mixing coordinates."
            )
            analysis_warnings.append(msg)
            print(f"  WARNING: {msg}")

        for mg_num, (sig, group_poses) in enumerate(sorted_sig_groups, 1):
            mg = f"MG{mg_num:02d}_{signature_hash(sig)}"
            for p in group_poses:
                p.mapping_group = mg

            n_atoms = len(sig)
            n = len(group_poses)
            pair_count = n * (n - 1) // 2
            print(f"  {mg}: {n:,} poses; {n_atoms} heavy atoms; {pair_count:,} pairwise RMSDs")

            clusters, dists = cluster_mapping_group(
                group_poses,
                cluster_rmsd=args.cluster_rmsd,
                ensemble_pose_count=len(ensemble_poses),
            )
            all_cluster_results.extend(clusters)
            print(f"    -> {len(clusters)} binding-mode cluster(s)")

            # Which clusters receive SASA in V1?
            if args.skip_sasa:
                sasa_cluster_ids: set[int] = set()
            elif args.sasa_top_clusters == 0:
                sasa_cluster_ids = {c.cluster_id for c in clusters}
            else:
                sasa_cluster_ids = {c.cluster_id for c in clusters[:max(0, args.sasa_top_clusters)]}

            # Representative files + cluster rows.
            for cluster in clusters:
                medoid_pose = group_poses[cluster.medoid_index]
                rep_subdir = reps_dir / safe_component(receptor) / safe_component(ligand_base)
                rep_name = f"{mg}_Cluster{cluster.cluster_id:03d}_medoid.pdbqt"
                rep_path = rep_subdir / rep_name
                try:
                    extract_model_to_file(medoid_pose.outfile, medoid_pose.pose, rep_path)
                except Exception as exc:
                    analysis_warnings.append(f"Could not write representative {rep_path}: {exc}")

                binding_rows.append({
                    "Receptor": receptor,
                    "LigandBase": ligand_base,
                    "MappingGroup": mg,
                    "MappingHeavyAtoms": n_atoms,
                    "ClusterID": cluster.cluster_id,
                    "PoseCount": cluster.pose_count,
                    "PoseSupportPct_Ensemble": fmt_pct(cluster.pose_support_ensemble_pct),
                    "PoseSupportPct_MappingGroup": fmt_pct(cluster.pose_support_mapping_pct),
                    "JobCount": cluster.job_count,
                    "JobSupportPct": fmt_pct(cluster.job_support_pct),
                    "TopPoseJobCount": cluster.top_pose_job_count,
                    "TopPoseJobSupportPct": fmt_pct(cluster.top_pose_job_support_pct),
                    "ConformerCount": cluster.conformer_count,
                    "ConformerSupportPct": fmt_pct(cluster.conformer_support_pct),
                    "VariantCount": cluster.variant_count,
                    "VariantSupportPct": fmt_pct(cluster.variant_support_pct),
                    "BestAffinity": round(cluster.best_affinity, 6),
                    "MedianAffinity": round(cluster.median_affinity, 6),
                    "MeanAffinity": round(cluster.mean_affinity, 6),
                    "MeanRMSDToMedoid_A": round(cluster.mean_rmsd_to_medoid, 6),
                    "MaxRMSDToMedoid_A": round(cluster.max_rmsd_to_medoid, 6),
                    "RepresentativeLigandVariant": medoid_pose.ligand_variant,
                    "RepresentativePose": medoid_pose.pose,
                    "RepresentativeAffinity": round(medoid_pose.affinity, 6),
                    "RepresentativeOutFile": str(medoid_pose.outfile),
                    "RepresentativePDBQT": str(rep_path),
                    "RepresentativeSeed": medoid_pose.seed,
                    "SASASelected": "yes" if cluster.cluster_id in sasa_cluster_ids else "no",
                    "SASAAnalyzedPoseCount": 0,
                })

            # SASA consensus for selected clusters.
            if not sasa_cluster_ids:
                continue

            receptor_path = next(
                (p.receptor_file for p in group_poses if p.receptor_file is not None and p.receptor_file.is_file()),
                None,
            )
            if receptor_path is None:
                msg = f"{receptor}/{ligand_base}/{mg}: receptor file unavailable; SASA skipped."
                analysis_warnings.append(msg)
                print(f"    WARNING: {msg}")
                continue

            rec_key = str(receptor_path.resolve())
            if rec_key not in receptor_cache:
                print(f"    Loading receptor environment for SASA: {receptor_path.name}")
                try:
                    receptor_cache[rec_key] = make_receptor_context(receptor_path, args.probe_radius)
                except Exception as exc:
                    msg = f"Could not build receptor SASA context from {receptor_path}: {exc}"
                    analysis_warnings.append(msg)
                    print(f"    WARNING: {msg}")
                    continue
            rec_ctx = receptor_cache[rec_key]

            for cluster in clusters:
                if cluster.cluster_id not in sasa_cluster_ids:
                    continue
                member_poses = [group_poses[i] for i in cluster.member_indices]
                print(
                    f"    SASA Cluster {cluster.cluster_id}: {len(member_poses):,} poses "
                    f"({cluster.pose_support_mapping_pct:.1f}% of {mg})"
                )

                # Each entry: arrays for all heavy atoms in canonical group order.
                free_matrix: List[np.ndarray] = []
                bound_matrix: List[np.ndarray] = []
                relative_matrix: List[np.ndarray] = []
                failures = 0

                for pose_i, pose in enumerate(member_poses, 1):
                    try:
                        free_all, bound_all, rel_all = ligand_relative_sasa(
                            pose.model.atoms if pose.model else [],
                            rec_ctx,
                            probe_radius=args.probe_radius,
                            sphere=sphere,
                        )
                        heavy_indices = [
                            idx for idx, atom in enumerate(pose.model.atoms if pose.model else [])
                            if atom.element.upper() != "H"
                        ]
                        free_matrix.append(free_all[heavy_indices])
                        bound_matrix.append(bound_all[heavy_indices])
                        relative_matrix.append(rel_all[heavy_indices])
                    except Exception as exc:
                        failures += 1
                        if failures <= 3:
                            print(f"      WARNING: SASA failed for {pose.ligand_variant} pose {pose.pose}: {exc}")

                    if pose_i % 100 == 0 or pose_i == len(member_poses):
                        sys.stdout.write(f"\r      SASA progress: {pose_i}/{len(member_poses)}")
                        sys.stdout.flush()
                print()

                if not relative_matrix:
                    analysis_warnings.append(
                        f"{receptor}/{ligand_base}/{mg}/Cluster{cluster.cluster_id}: all SASA calculations failed."
                    )
                    continue

                free_arr = np.stack(free_matrix, axis=0)
                bound_arr = np.stack(bound_matrix, axis=0)
                rel_arr = np.stack(relative_matrix, axis=0)
                analyzed_n = int(rel_arr.shape[0])
                cluster.sasa_analyzed = True
                cluster.sasa_pose_count = analyzed_n

                # Update corresponding binding row.
                for br in reversed(binding_rows):
                    if (
                        br["Receptor"] == receptor
                        and br["LigandBase"] == ligand_base
                        and br["MappingGroup"] == mg
                        and br["ClusterID"] == cluster.cluster_id
                    ):
                        br["SASAAnalyzedPoseCount"] = analyzed_n
                        break

                ref_atoms = member_poses[0].heavy_atoms
                for atom_idx in range(rel_arr.shape[1]):
                    vals = rel_arr[:, atom_idx]
                    free_vals = free_arr[:, atom_idx]
                    bound_vals = bound_arr[:, atom_idx]
                    q1, q3 = np.percentile(vals, [25, 75])
                    n25 = int(np.sum(vals >= 25.0))
                    n50 = int(np.sum(vals >= 50.0))
                    n75 = int(np.sum(vals >= 75.0))
                    atom = ref_atoms[atom_idx]
                    atom_rows.append({
                        "Receptor": receptor,
                        "LigandBase": ligand_base,
                        "MappingGroup": mg,
                        "ClusterID": cluster.cluster_id,
                        "ClusterPoseSupportPct_Ensemble": fmt_pct(cluster.pose_support_ensemble_pct),
                        "ClusterPoseSupportPct_MappingGroup": fmt_pct(cluster.pose_support_mapping_pct),
                        "ClusterJobSupportPct": fmt_pct(cluster.job_support_pct),
                        "ClusterConformerSupportPct": fmt_pct(cluster.conformer_support_pct),
                        "AtomIndex": atom_idx + 1,
                        "AtomName": atom.name,
                        "Element": atom.element,
                        "PDBQTType": atom.atom_type,
                        "SASA_PoseCount": analyzed_n,
                        "MedianFreeSASA_A2": round(float(np.median(free_vals)), 6),
                        "MedianBoundSASA_A2": round(float(np.median(bound_vals)), 6),
                        "MeanRelativeSASA_Pct": round(float(np.mean(vals)), 6),
                        "MedianRelativeSASA_Pct": round(float(np.median(vals)), 6),
                        "RelativeSASA_Q1_Pct": round(float(q1), 6),
                        "RelativeSASA_Q3_Pct": round(float(q3), 6),
                        "ExposedAt25Pct_Poses": n25,
                        "ExposedAt25Pct_Pct": fmt_pct(pct(n25, analyzed_n)),
                        "ExposedAt50Pct_Poses": n50,
                        "ExposedAt50Pct_Pct": fmt_pct(pct(n50, analyzed_n)),
                        "ExposedAt75Pct_Poses": n75,
                        "ExposedAt75Pct_Pct": fmt_pct(pct(n75, analyzed_n)),
                    })

    # Stable output order.
    poses.sort(key=lambda p: (p.receptor.lower(), p.ligand_base.lower(), p.mapping_group, p.cluster_id, p.affinity, p.job_key, p.pose))
    binding_rows.sort(key=lambda r: (str(r["Receptor"]).lower(), str(r["LigandBase"]).lower(), str(r["MappingGroup"]), int(r["ClusterID"])))
    atom_rows.sort(key=lambda r: (str(r["Receptor"]).lower(), str(r["LigandBase"]).lower(), str(r["MappingGroup"]), int(r["ClusterID"]), int(r["AtomIndex"])))

    assignments_path = output_dir / "Protacability_PoseAssignments.csv"
    modes_path = output_dir / "Protacability_BindingModes.csv"
    atoms_path = output_dir / "Protacability_AtomConsensus.csv"
    summary_path = output_dir / "Protacability_Summary.json"

    write_pose_assignments(assignments_path, poses)
    write_binding_modes(modes_path, binding_rows)
    write_atom_consensus(atoms_path, atom_rows)

    docking_provenance = None
    if provenance_path.is_file():
        try:
            docking_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        except Exception as exc:
            analysis_warnings.append(f"Could not parse docking_parameters.json: {exc}")

    elapsed = (datetime.now() - start).total_seconds()
    summary = {
        "schema_version": 1,
        "analysis": "Consensus PROTACability V1",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "runtime_seconds": round(elapsed, 3),
        "inputs": {
            "docking_results_dir": str(results_dir),
            "step4_scores_csv": str(scores_csv),
            "docking_parameters_json": str(provenance_path) if provenance_path.is_file() else None,
            "docking_jobs_csv": str(results_dir / "docking_jobs.csv") if job_manifest else None,
        },
        "parameters": {
            "cluster_rmsd_angstrom": args.cluster_rmsd,
            "rmsd_coordinate_frame": "receptor-frame; no ligand superposition",
            "mapping_policy": "strict PDBQT heavy-atom name/element/type/order signature",
            "probe_radius_angstrom": args.probe_radius,
            "sasa_points_per_atom": args.sasa_points,
            "sasa_top_clusters_per_mapping_group": args.sasa_top_clusters,
            "sasa_skipped": bool(args.skip_sasa),
            "exposure_thresholds_relative_sasa_pct": [25, 50, 75],
        },
        "counts": {
            "poses_loaded": len(poses),
            "receptor_ligand_ensembles": len(ensembles),
            "binding_mode_clusters": len(binding_rows),
            "atom_consensus_rows": len(atom_rows),
        },
        "interpretation": {
            "cluster_support": "docking ensemble support, not thermodynamic probability",
            "sasa": "relative bound/free SASA from docking-coordinate atoms; heavy-atom rows reported",
        },
        "warnings": analysis_warnings,
        "docking_provenance": docking_provenance,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("\n=== Step 8 complete ===")
    print(f"Output directory : {output_dir}")
    print(f"Binding modes    : {modes_path}")
    print(f"Pose assignments : {assignments_path}")
    print(f"Atom consensus   : {atoms_path}")
    print(f"Summary JSON     : {summary_path}")
    print(f"Representatives  : {reps_dir}")
    print(f"Runtime          : {elapsed:.1f} seconds")

    if analysis_warnings:
        print(f"Warnings         : {len(analysis_warnings)} (see Protacability_Summary.json)")


if __name__ == "__main__":
    main()
