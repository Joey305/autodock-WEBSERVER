#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import csv
import re
import argparse
import shlex
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Tuple

from ligand_manifest import CHEMICAL_METADATA_COLUMNS, merge_ligand_metadata
from ligand_naming import parse_ligand_variant

REQUIRED = ["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile"]

TIER_TO_TOP = {
    0: 0,
    1: 100,
    2: 1000,
    3: 10000,
}


def find_result_dirs(base: Path) -> List[Path]:
    return [d for d in sorted(base.glob("Docking_Results_*")) if d.is_dir()]


def latest_perdir_csv(results_dir: Path) -> Optional[Path]:
    # Matches per-dir sorted output you already generate
    pat = f"{results_dir.name}*_vina_docking_scores_sorted.csv"
    cands = sorted(results_dir.parent.glob(pat))
    return cands[-1] if cands else None


def choose_multi(items: List[Path]) -> List[Path]:
    print("\nSelect Docking_Results_* directories to combine:")
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
                    idxs.extend(range(int(a), int(b) + 1))
                else:
                    idxs.append(int(tok))

            idxs = [i - 1 for i in idxs]
            sel = [items[i] for i in idxs if 0 <= i < len(items)]
            if sel:
                return sel
        except Exception:
            pass

        print("Invalid input. Try again.")


def ask_top_tier() -> int:
    print("\nHow many TOP rows per directory?")
    print(" [0] ALL")
    print(" [1] 100")
    print(" [2] 1000")
    print(" [3] 10000")

    while True:
        s = input("Choose 0/1/2/3 [default 1=100]: ").strip() or "1"
        try:
            v = int(s)
            if v in (0, 1, 2, 3):
                return TIER_TO_TOP[v]
        except ValueError:
            pass
        print("Invalid choice. Try again.")


def derive_provenance(outfile_path: str):
    """
    From:
      /.../Docking_Results_XYZ/<ReceptorDir>/<LigandDir>/out.pdbqt

    Return:
      (results_root, receptor_dir, ligand_dir, parsed_metadata)

    IMPORTANT:
      ligand metadata is derived from the full ligand folder name using the shared parser.
      We do NOT split on internal '__', '____', underscores, or hyphens.
    """
    try:
        p = Path(outfile_path)
        ligand_dir = p.parent.name
        receptor_dir = p.parent.parent.name
        results_root = p.parent.parent.parent.name
        parsed = parse_ligand_variant(ligand_dir)
        return results_root, receptor_dir, ligand_dir, parsed
    except Exception:
        return "", "", "", parse_ligand_variant("")


def pick_variant_name(row: dict, ligand_dir: str) -> str:
    for key in ("LigandVariant", "LigandDir"):
        value = (row.get(key) or "").strip()
        if value:
            return value
    if ligand_dir:
        return ligand_dir
    value = (row.get("Ligand") or "").strip()
    if value:
        return value
    return ligand_dir


def enrich_row(row: dict) -> dict:
    results_root, receptor_dir, ligand_dir, parsed_from_outfile = derive_provenance(row["OutFile"])
    variant = pick_variant_name(row, ligand_dir)
    merged = merge_ligand_metadata(variant, base_row=row)
    ligand_base = (merged.get("LigandBase") or row.get("Ligand") or variant).strip()
    ligand_variant = (merged.get("LigandVariant") or variant).strip()
    receptor = (row.get("Receptor") or receptor_dir).strip()
    enriched = {
        "Receptor": receptor,
        "Ligand": ligand_base,
        "LigandBase": ligand_base,
        "LigandVariant": ligand_variant,
        "Pose": row.get("Pose", ""),
        "Binding_Affinity": row.get("Binding_Affinity", ""),
        "OutFile": row.get("OutFile", ""),
        "ResultsRoot": results_root,
        "ReceptorDir": receptor_dir or receptor,
        "LigandDir": ligand_dir or ligand_variant,
        "ProtomerTag": merged.get("ProtomerTag") or parsed_from_outfile["ProtomerTag"],
        "TautomerTag": merged.get("TautomerTag") or parsed_from_outfile["TautomerTag"],
        "ConformerTag": merged.get("ConformerTag") or parsed_from_outfile["ConformerTag"],
        "LegacyPoseTag": merged.get("LegacyPoseTag") or parsed_from_outfile["LegacyPoseTag"],
        "StateTag": merged.get("StateTag") or parsed_from_outfile["StateTag"],
        "ProtomerIndex": merged.get("ProtomerIndex") or parsed_from_outfile["ProtomerIndex"],
        "TautomerIndex": merged.get("TautomerIndex") or parsed_from_outfile["TautomerIndex"],
        "ConformerIndex": merged.get("ConformerIndex") or parsed_from_outfile["ConformerIndex"],
    }
    for key in CHEMICAL_METADATA_COLUMNS:
        enriched[key] = merged.get(key, "")
    return enriched


def external_sort_csv(in_path: Path, out_path: Path, header_cols: List[str]):
    sort_exe = shutil.which("sort")
    header = ",".join(header_cols)

    if not sort_exe:
        # Fallback (RAM sort) only if 'sort' is unavailable
        rows = []
        with open(in_path, newline="", encoding="utf-8") as f:
            r = csv.reader(f)
            next(r, None)
            for row in r:
                rows.append(row)

        # Sort by Receptor (col 1) then Binding_Affinity (col 6 numeric ascending)
        rows.sort(
            key=lambda r: (
                r[0],
                float(r[5]) if r[5] not in ("", None) else 1e9
            )
        )

        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(header_cols)
            w.writerows(rows)
        return

    tmp_sorted = in_path.with_suffix(".nosort.tmp")
    cmd = (
        f"tail -n +2 {shlex.quote(str(in_path))} | "
        f"{sort_exe} -t, -k1,1 -k6,6n > {shlex.quote(str(tmp_sorted))}"
    )
    subprocess.run(["bash", "-lc", cmd], check=True)

    with open(out_path, "w", encoding="utf-8", newline="") as out_f:
        out_f.write(header + "\n")
        with open(tmp_sorted, "r", encoding="utf-8", errors="ignore") as src:
            shutil.copyfileobj(src, out_f)

    try:
        tmp_sorted.unlink()
        in_path.unlink()
    except Exception:
        pass


def main():
    ap = argparse.ArgumentParser(
        description=(
            "4C_ConcatenateScores: concatenate per-dir Vina CSVs with provenance columns, "
            "with optional TOP-N-per-dir limiter."
        )
    )
    ap.add_argument(
        "--dirs",
        nargs="*",
        help="Specific Docking_Results_* dirs. Omit for interactive pick."
    )
    ap.add_argument(
        "--outfile",
        help="Output CSV path (sorted). Default includes timestamp and TOP value."
    )
    ap.add_argument(
        "--tier",
        type=int,
        choices=[0, 1, 2, 3],
        help="0=ALL, 1=100, 2=1000, 3=10000"
    )
    ap.add_argument(
        "--top",
        type=int,
        default=None,
        help="Override: exact number of top rows per directory (0=ALL)."
    )
    ap.add_argument(
        "--no-sort",
        action="store_true",
        help="Write unsorted CSV (faster, but not ranked globally)."
    )
    args = ap.parse_args()

    base = Path(".").resolve()
    all_dirs = find_result_dirs(base)

    if not all_dirs and not args.dirs:
        print("❌ No Docking_Results_* directories found.")
        sys.exit(1)

    selected = [Path(d).resolve() for d in args.dirs] if args.dirs else choose_multi(all_dirs)

    # Determine TOP per-dir
    if args.top is not None:
        top = max(0, int(args.top))
    elif args.tier is not None:
        top = TIER_TO_TOP[args.tier]
    else:
        top = ask_top_tier()  # interactive default 100

    # Find per-dir CSVs
    csvs = []
    for d in selected:
        c = latest_perdir_csv(d)
        if not c:
            print(f"⚠️ No per-dir CSV found for {d.name} — skipping")
            continue
        csvs.append(c)

    if not csvs:
        print("❌ No per-dir CSVs to combine.")
        sys.exit(2)

    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    top_tag = f"TOP{top}" if top > 0 else "TOPALL"

    header = [
        "Receptor",
        "Ligand",
        "LigandBase",
        "LigandVariant",
        "Pose",
        "Binding_Affinity",
        "OutFile",
        "ResultsRoot",
        "ReceptorDir",
        "LigandDir",
        "ProtomerTag",
        "TautomerTag",
        "ConformerTag",
        "LegacyPoseTag",
        "StateTag",
        "ProtomerIndex",
        "TautomerIndex",
        "ConformerIndex",
    ] + CHEMICAL_METADATA_COLUMNS

    tmp = base / f"ALL_Docking_Results_with_provenance_{top_tag}_{ts}_UNSORTED.tmp.csv"
    out = Path(args.outfile) if args.outfile else (
        base / f"ALL_Docking_Results_with_provenance_{top_tag}_{ts}.csv"
    )

    total_rows = 0

    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()

        for c in csvs:
            taken = 0
            with open(c, newline="", encoding="utf-8") as f:
                r = csv.DictReader(f)

                # sanity check
                for k in REQUIRED:
                    if k not in r.fieldnames:
                        raise SystemExit(f"❌ {c.name} missing column {k}; has: {r.fieldnames}")

                for row in r:
                    w.writerow(enrich_row(row))

                    total_rows += 1
                    taken += 1

                    if top > 0 and taken >= top:
                        break  # respect TOP per-dir

    if args.no_sort:
        tmp.rename(out)
        print(f"✅ Wrote UNSORTED combined CSV with provenance → {out}")
        print(
            f"   Rows: {total_rows}  |  From per-dir files: {len(csvs)}  |  "
            f"TOP per-dir: {top if top > 0 else 'ALL'}"
        )
        return

    # External sort by Receptor (col1) then Binding_Affinity (col6 numeric)
    external_sort_csv(tmp, out, header_cols=header)

    print(f"✅ Wrote combined (sorted) CSV with provenance → {out}")
    print(
        f"   Rows: {total_rows}  |  From per-dir files: {len(csvs)}  |  "
        f"TOP per-dir: {top if top > 0 else 'ALL'}"
    )


if __name__ == "__main__":
    main()
