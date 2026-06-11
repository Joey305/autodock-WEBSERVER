#!/usr/bin/env python3
from __future__ import annotations

import os, sys, re, csv, argparse, threading
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Optional, Tuple

from ligand_manifest import (
    CHEMICAL_METADATA_COLUMNS,
    find_ligand_state_manifests,
    load_ligand_state_manifest,
    merge_ligand_metadata,
)
from ligand_naming import parse_ligand_variant

# ==========================================================
# Command-line argument parser
# ==========================================================
def parse_cli():
    ap = argparse.ArgumentParser(description="4_ParseScores: parse Vina results for one or more Docking_Results_* dirs.")
    ap.add_argument("--dir", help="Single Docking_Results_* directory to process (non-interactive).")
    ap.add_argument("--workers", type=int, default=8, help="Worker threads per job (default: 8).")
    ap.add_argument("--heartbeat", type=int, default=1000, help="Files per progress update.")
    ap.add_argument("--fallback-crawl", action="store_true", help="If log missing, crawl directory.")
    return ap.parse_args()

# ==========================================================
# Regex patterns
# ==========================================================
VINA_RE = re.compile(r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)")
OK_LINE = re.compile(r"^[\s✅\-\+]?\s*(.+?)\s*→\s*([A-Za-z0-9._-]+)\s*$")  # ligand  → receptor.pdbqt

# ==========================================================
# Helper functions
# ==========================================================
def parse_vina_pdbqt(p: Path):
    """Parse a PDBQT file for all VINA RESULT lines."""
    poses = []
    try:
        with open(p, "r", errors="ignore") as fh:
            pose_idx = 0
            for line in fh:
                m = VINA_RE.match(line)
                if m:
                    pose_idx += 1
                    poses.append((pose_idx, float(m.group(1))))
    except Exception:
        pass
    return poses

def find_candidates() -> List[Path]:
    """Return all Docking_Results_* directories in the current folder."""
    return [d for d in sorted(Path(".").glob("Docking_Results_*")) if d.is_dir()]

def match_run_log(results_dir: Path) -> Optional[Path]:
    """Find the corresponding run_log_*.txt file for a given results dir."""
    rest = results_dir.name[len("Docking_Results_"):]
    exact = results_dir.parent / f"run_log_{rest}.txt"
    if exact.is_file():
        return exact
    # fuzzy matching
    cands = sorted(results_dir.parent.glob(f"run_log_{rest}*.txt"))
    return cands[0] if cands else None

def parse_log_to_targets(results_dir: Path, log_path: Path) -> List[Tuple[str, str, Path]]:
    """Extract receptor-ligand-outfile mappings from a run_log_*.txt."""
    ok_pairs = []
    with open(log_path, "r", encoding="utf-8", errors="ignore") as lf:
        for line in lf:
            if "✅" not in line and "->" not in line and "→" not in line:
                continue
            line = line.replace("->", "→")
            m = OK_LINE.match(line.strip())
            if not m:
                continue
            ligand_dir = m.group(1).strip()
            receptor_file = m.group(2).strip()
            receptor_dir = receptor_file.replace(".pdbqt", "").replace(".pdb", "").replace(".mol2", "")
            p = results_dir / receptor_dir / ligand_dir / "out.pdbqt"
            ok_pairs.append((receptor_dir, ligand_dir, p))
    return ok_pairs

def choose_multi(items: List[Path]) -> List[Path]:
    """Interactive picker for Docking_Results_* directories."""
    print("\nSelect Docking_Results_* directories:")
    print(" [0] ALL")
    for i, d in enumerate(items, 1):
        print(f" [{i}] {d.name}")
    while True:
        s = input("Enter indexes (e.g., 1,2 or 3-7; 0=ALL): ").strip()
        if s == "0":
            return items
        try:
            idxs = []
            for tok in s.split(","):
                tok = tok.strip()
                if not tok:
                    continue
                if "-" in tok:
                    a, b = tok.split("-", 1)
                    a, b = int(a), int(b)
                    idxs.extend(range(a, b + 1))
                else:
                    idxs.append(int(tok))
            idxs = [i - 1 for i in idxs]
            sel = [items[i] for i in idxs if 0 <= i < len(items)]
            if sel:
                return sel
        except Exception:
            pass
        print("Invalid input. Try again.")

# ==========================================================
# Core processing
# ==========================================================
OUTPUT_COLUMNS = [
    "Receptor",
    "Ligand",
    "LigandBase",
    "LigandVariant",
    "Pose",
    "Binding_Affinity",
    "OutFile",
    "ProtomerTag",
    "TautomerTag",
    "ConformerTag",
    "LegacyPoseTag",
    "StateTag",
    "ProtomerIndex",
    "TautomerIndex",
    "ConformerIndex",
] + CHEMICAL_METADATA_COLUMNS


def load_manifest_map(search_roots: List[Path]) -> Dict[str, Dict[str, str]]:
    merged: Dict[str, Dict[str, str]] = {}
    for manifest_path in find_ligand_state_manifests(search_roots):
        merged.update(load_ligand_state_manifest(manifest_path))
    return merged


def build_score_rows(
    receptor: str,
    ligand_variant: str,
    outfile_path: Path,
    poses: List[Tuple[int, float]],
    manifest_row: Optional[Dict[str, str]] = None,
) -> List[Dict[str, object]]:
    parsed = parse_ligand_variant(ligand_variant)
    metadata = merge_ligand_metadata(ligand_variant, manifest_row=manifest_row)
    rows = []
    for pose_idx, aff in poses:
        row = {
            "Receptor": receptor,
            "Ligand": metadata["LigandBase"] or parsed["LigandBase"],
            "LigandBase": metadata["LigandBase"] or parsed["LigandBase"],
            "LigandVariant": metadata["LigandVariant"] or parsed["LigandVariant"],
            "Pose": pose_idx,
            "Binding_Affinity": float(aff),
            "OutFile": str(outfile_path.resolve()),
            "ProtomerTag": metadata["ProtomerTag"],
            "TautomerTag": metadata["TautomerTag"],
            "ConformerTag": metadata["ConformerTag"],
            "LegacyPoseTag": metadata["LegacyPoseTag"],
            "StateTag": metadata["StateTag"],
            "ProtomerIndex": metadata["ProtomerIndex"],
            "TautomerIndex": metadata["TautomerIndex"],
            "ConformerIndex": metadata["ConformerIndex"],
        }
        for key in CHEMICAL_METADATA_COLUMNS:
            row[key] = metadata.get(key, "")
        rows.append(row)
    return rows


def process_one_dir(results_dir: Path, workers: int, heartbeat: int, fallback_crawl: bool) -> List[Dict[str, object]]:
    """Parse all out.pdbqt files for one Docking_Results_* directory."""
    log_path = match_run_log(results_dir)
    targets: List[Tuple[str, str, Path]] = []
    manifest_map = load_manifest_map([results_dir.parent, results_dir, Path.cwd()])

    if log_path and log_path.is_file():
        print(f"✅ {results_dir.name} -> {log_path.name}")
        targets = parse_log_to_targets(results_dir, log_path)
    elif fallback_crawl:
        print(f"🐢 {results_dir.name}: no log; crawling …")
        for root, dirs, files in os.walk(results_dir):
            dirs[:] = [x for x in dirs if not x.startswith(".")]
            if "out.pdbqt" in files:
                p = Path(root) / "out.pdbqt"
                receptor = p.parent.parent.name
                ligdir = p.parent.name
                targets.append((receptor, ligdir, p))
    else:
        print(f"⏭️  {results_dir.name}: no matching run_log_; skipping (use --fallback-crawl to include).")
        return []

    if not targets:
        print(f"⚠️ {results_dir.name}: 0 candidate out.pdbqt")
        return []

    rows_buffer = []
    processed = 0
    hits = 0
    lock = threading.Lock()

    def work(tup):
        receptor, ligand_dir, p = tup
        if not Path(p).is_file():
            return []
        poses = parse_vina_pdbqt(p)
        if not poses:
            return []
        return build_score_rows(receptor, ligand_dir, Path(p), poses, manifest_map.get(ligand_dir))

    print(f"🚀 {results_dir.name}: parsing {len(targets)} files with {workers} threads …")
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(work, t) for t in targets]
        for fut in as_completed(futs):
            try:
                part = fut.result()
            except Exception as e:
                print(f"❌ Worker error: {e}")
                continue
            processed += 1
            if part:
                hits += 1
                with lock:
                    rows_buffer.extend(part)
            if processed % heartbeat == 0:
                pct = int(processed * 100 / max(1, len(targets)))
                sys.stdout.write(f"\r🔎 {results_dir.name}: {processed}/{len(targets)} ({pct}%)  files_with_results={hits}")
                sys.stdout.flush()

    print(f"\n✅ {results_dir.name}: parsed {processed} files; {hits} had results.")
    return rows_buffer

# ==========================================================
# Main entry point
# ==========================================================
def main():
    args = parse_cli()

    # Non-interactive mode for LSF or automation
    if args.dir:
        results_dir = Path(args.dir).resolve()
        if not results_dir.is_dir():
            print(f"❌ Invalid --dir: {results_dir}")
            sys.exit(1)

        ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        rows = process_one_dir(results_dir, args.workers, args.heartbeat, args.fallback_crawl)
        if not rows:
            print(f"⚠️ No results parsed in {results_dir.name}")
            sys.exit(0)

        rows.sort(key=lambda r: (str(r["Receptor"]), float(r["Binding_Affinity"])))
        per_dir_csv = results_dir.parent / f"{results_dir.name}_{ts}_vina_docking_scores_sorted.csv"
        with open(per_dir_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS)
            w.writeheader()
            w.writerows(rows)
        print(f"📊 Per-dir CSV → {per_dir_csv}")
        sys.exit(0)

    # Interactive multi-directory mode
    candidates = find_candidates()
    if not candidates:
        print("❌ No Docking_Results_* directories found.")
        sys.exit(1)

    selected = choose_multi(candidates)
    selected = [Path(p).resolve() for p in selected]
    print("\n▶️ Selected dirs:")
    for p in selected:
        print("  -", p)

    default_workers = min(64, (os.cpu_count() or 8) * 2)
    workers = args.workers or default_workers
    heartbeat = max(50, args.heartbeat)

    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    all_rows: List[Dict[str, object]] = []

    for results_dir in selected:
        rows = process_one_dir(results_dir, workers, heartbeat, args.fallback_crawl)
        if not rows:
            continue
        rows.sort(key=lambda r: (str(r["Receptor"]), float(r["Binding_Affinity"])))
        per_dir_csv = results_dir.parent / f"{results_dir.name}_{ts}_vina_docking_scores_sorted.csv"
        with open(per_dir_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS)
            w.writeheader()
            w.writerows(rows)
        print(f"📊 Per-dir CSV → {per_dir_csv}")
        all_rows.extend(rows)

    if not all_rows:
        print("⚠️ No rows parsed from selected directories.")
        sys.exit(2)

    # Optional: combined CSV (disabled by default for memory safety)
    # all_rows.sort(key=lambda r: (r[0], float(r[3])))
    # all_csv = Path(".").resolve() / f"ALL_Docking_Results_{ts}_vina_docking_scores_ALL_sorted.csv"
    # with open(all_csv, "w", newline="", encoding="utf-8") as fh:
    #     w = csv.writer(fh)
    #     w.writerow(["Receptor","Ligand","Pose","Binding_Affinity","OutFile"])
    #     w.writerows(all_rows)
    # print(f"🏁 ALL-IN-ONE CSV → {all_csv}")

if __name__ == "__main__":
    main()
