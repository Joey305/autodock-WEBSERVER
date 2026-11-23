#!/usr/bin/env python3
"""Report live, truthful progress for a generated docking campaign."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


def _count_rows(path: Path) -> int:
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            return sum(1 for _ in csv.DictReader(handle))
    except OSError:
        return 0


def _unique_jobs(path: Path) -> set[tuple[str, str, str]]:
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            return {(r.get("ResolvedReceptor", ""), r.get("Ligand", ""), r.get("Seed", "")) for r in csv.DictReader(handle) if r.get("Ligand")}
    except OSError:
        return set()


def _job_state(job_name: str) -> str:
    try:
        run = subprocess.run(["bjobs", "-noheader", "-o", "stat", "-J", job_name], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=8, check=False)
        states = [line.split()[0] for line in run.stdout.splitlines() if line.strip()]
        return states[-1] if states else "DONE"
    except (OSError, subprocess.TimeoutExpired):
        return "UNKNOWN"


def _pdbqt_count(batch: Path) -> int:
    return sum(1 for path in batch.glob("*_Ligands_PDBQT_*/*.pdbqt") if path.is_file())


def _input_count(batch: Path) -> int:
    inputs = [path for path in batch.iterdir() if path.is_file() and path.suffix.lower() in {".csv", ".sdf", ".smi", ".smiles"}]
    if not inputs:
        return 0
    source = inputs[0]
    if source.suffix.lower() == ".csv":
        with source.open(newline="", encoding="utf-8-sig", errors="replace") as handle:
            return sum(1 for row in csv.reader(handle) if any(cell.strip() for cell in row)) - 1
    if source.suffix.lower() == ".sdf":
        return sum(1 for line in source.open(encoding="utf-8", errors="replace") if line.strip() == "$$$$")
    return sum(1 for line in source.open(encoding="utf-8", errors="replace") if line.strip())


def _center_count(root: Path) -> int:
    centers = sorted(root.rglob("vina_centers*.csv"))
    if not centers:
        return 0
    with centers[0].open(newline="", encoding="utf-8-sig", errors="replace") as handle:
        return sum(1 for row in csv.DictReader(handle) if any((value or "").strip() for value in row.values()))


def _legacy_maximums(settings: dict[str, Any], input_ligands: int, docking_targets: int) -> tuple[int, int]:
    flags = list(settings.get("state_flags", []))
    enumerate_states = bool(settings.get("enumerate_protomers")) or "--enumerate-protomers" in flags
    def flag_value(flag: str, default: int) -> int:
        if flag in flags:
            try:
                return int(flags[flags.index(flag) + 1])
            except (IndexError, ValueError):
                pass
        return int(settings.get(flag.removeprefix("--").replace("-", "_"), default) or default)
    states = flag_value("--max-protomers", 4) * flag_value("--max-tautomers", 4) if enumerate_states else 1
    pdbqt_maximum = input_ligands * int(settings.get("confgen_poses", 64) or 64) * states
    return pdbqt_maximum, pdbqt_maximum * docking_targets


def _pct(done: int, total: int | None) -> float | None:
    return round(100 * done / total, 1) if total else None


def campaign_progress(root: Path) -> dict[str, Any]:
    path = root / "campaign.json"
    if not path.is_file():
        raise ValueError(f"campaign.json was not found in {root}")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    batches = manifest.get("batches", [])
    settings = manifest.get("settings", {})
    docking_targets = _center_count(root)
    if batches and isinstance(batches[0], str):  # early campaign compatibility
        batches = [{"name": name} for name in batches]
    result, pdbqt_done, pdbqt_max, vina_done, vina_total, vina_exact = [], 0, 0, 0, 0, True
    for index, item in enumerate(batches, start=1):
        name = str(item.get("name") or f"batch_{index:03d}")
        batch, results = root / "batches" / name, root / "results" / name
        input_ligands = int(item.get("input_ligands", 0) or 0) or _input_count(batch)
        maximum = int(item.get("pdbqt_maximum", 0) or 0)
        planned = int(item.get("vina_maximum", 0) or 0)
        if not maximum:
            maximum, planned = _legacy_maximums(settings, input_ligands, docking_targets)
        produced = _pdbqt_count(batch)
        actual = _count_rows(results / "docking_jobs.csv")
        total, exact = actual or planned, bool(actual)
        success, failed = _unique_jobs(results / "compact_results.csv"), _unique_jobs(results / "compact_failures.csv")
        completed = success | failed
        result.append({
            "name": name,
            "input_ligands": input_ligands,
            "confgen_status": _job_state(f"{manifest.get('campaign', '')}-conf-{index:03d}"),
            "vina_status": _job_state(f"{manifest.get('campaign', '')}-vina-{index:03d}"),
            "pdbqt": {"produced": produced, "maximum": maximum, "percent_of_maximum": _pct(produced, maximum)},
            "vina": {"completed": len(completed), "successful": len(success), "failed": len(failed), "total": total, "total_is_exact": exact, "percent": _pct(len(completed), total)},
        })
        pdbqt_done += produced; pdbqt_max += maximum; vina_done += len(completed); vina_total += total; vina_exact = vina_exact and exact
    return {"campaign": manifest.get("campaign", root.name), "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "notes": {"pdbqt_target": "Maximum planned PDBQT files; valid molecules can yield fewer states or conformers.", "vina_target": "Exact after docking_jobs.csv exists; before then it is the maximum planned count."}, "pdbqt": {"produced": pdbqt_done, "maximum": pdbqt_max, "percent_of_maximum": _pct(pdbqt_done, pdbqt_max)}, "vina": {"completed": vina_done, "total": vina_total, "total_is_exact": vina_exact, "percent": _pct(vina_done, vina_total)}, "batches": result}


def _show(report: dict[str, Any]) -> None:
    p, v = report["pdbqt"], report["vina"]
    print(f"Campaign: {report['campaign']}")
    print(f"PDBQT generation: {p['produced']:,}/{p['maximum']:,} maximum planned ({p['percent_of_maximum'] or 0:.1f}%)")
    print(f"Vina jobs: {v['completed']:,}/{v['total']:,} {'exact total' if v['total_is_exact'] else 'maximum planned total'} ({v['percent'] or 0:.1f}%)")
    print("\nBatch       ConfGen  Vina     PDBQT produced / max        Docking complete / total")
    for b in report["batches"]:
        p, v = b["pdbqt"], b["vina"]
        print(f"{b['name']:<11} {b['confgen_status']:<8} {b['vina_status']:<8} {p['produced']:>7,}/{p['maximum']:<7,} ({p['percent_of_maximum'] or 0:>5.1f}%)  {v['completed']:>7,}/{v['total']:<7,} ({v['percent'] or 0:>5.1f}%)")


def main() -> int:
    parser = argparse.ArgumentParser(description="Show PDBQT and Vina progress for a docking campaign.")
    parser.add_argument("--root", default=Path(__file__).resolve().parent, type=Path)
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--write", action="store_true", help="Update campaign_status.json atomically.")
    parser.add_argument("--watch", type=float, metavar="SECONDS")
    args = parser.parse_args()
    while True:
        try:
            report = campaign_progress(args.root.resolve())
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"campaign status error: {exc}", file=sys.stderr); return 2
        if args.write:
            target = args.root / "campaign_status.json"; temp = target.with_suffix(".json.tmp")
            temp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8"); temp.replace(target)
        print(json.dumps(report, indent=2) if args.json else "") if args.json else _show(report)
        if not args.watch: return 0
        time.sleep(max(args.watch, 1.0))


if __name__ == "__main__":
    raise SystemExit(main())
