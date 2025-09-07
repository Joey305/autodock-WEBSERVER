#!/usr/bin/env python3
import os, sys, re, csv, argparse, time, threading
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

VINA_RE = re.compile(r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)")

# ---------- helpers ----------
def parse_vina_pdbqt(p: Path):
    """Return list of (pose_idx, affinity) for this PDBQT. Empty if none."""
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
        return []
    return poses

def guess_ligand_name(p: Path):
    stem = p.stem
    for suf in ("_out", "_OUT", "-out"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
    m = re.match(r"(.+?)_pose\d+$", stem)
    if m:
        stem = m.group(1)
    return stem

def choose(prompt: str, items: list[Path]) -> Path:
    print(f"\n{prompt}")
    for i, item in enumerate(items):
        print(f" [{i}] {item.name}")
    while True:
        s = input("Enter index: ").strip()
        try:
            i = int(s)
            if 0 <= i < len(items):
                return items[i]
        except ValueError:
            pass
        print("Invalid choice. Try again.")

def choose_scope(results_dir: Path) -> list[Path]:
    subs = [d for d in sorted(results_dir.iterdir()) if d.is_dir()]
    if not subs:
        print(f"⚠️ No subfolders under {results_dir.name}. Using it directly.")
        return [results_dir]
    print(f"\nSelect subfolder(s) under {results_dir.name}:")
    print(" [0] All subfolders")
    for i, d in enumerate(subs, 1): print(f" [{i}] {d.name}")
    while True:
        s = input("Enter index or 0 for ALL: ").strip()
        if s == "0": return subs
        try:
            i = int(s)
            if 1 <= i <= len(subs): return [subs[i-1]]
        except ValueError:
            pass
        print("Invalid choice. Try again.")

def find_candidates() -> list[Path]:
    return [d for d in sorted(Path(".").glob("Docking_Results_*")) if d.is_dir()]

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Parse Vina results into CSV (with live streaming & logging).")
    ap.add_argument("--fast", action="store_true",
                    help="Minimal extras; omit timestamp/mtime columns.")
    ap.add_argument("--workers", type=int, default=min(32, (os.cpu_count() or 8)*2),
                    help="Worker threads (default ≈ 2× CPU cores).")
    ap.add_argument("--live", action="store_true",
                    help="Stream rows to a LIVE CSV as they are found (tail -f friendly).")
    ap.add_argument("--log", action="store_true",
                    help="Write a progress log that includes found results and periodic heartbeats.")
    ap.add_argument("--heartbeat", type=int, default=1000,
                    help="Log a heartbeat every N files processed (default: 1000).")
    ap.add_argument("--name", default="out.pdbqt",
                    help="Only parse files with this exact name (default: out.pdbqt). Use '*.pdbqt' to scan all, slower.")
    ap.add_argument("--no-sort", action="store_true",
                    help="Skip final sort (useful if you only want the LIVE CSV).")
    args = ap.parse_args()

    candidates = find_candidates()
    if not candidates:
        print("❌ No Docking_Results_* directories in current folder.")
        sys.exit(1)

    results_dir = choose("Select a docking results folder:", candidates)
    scope_dirs = choose_scope(results_dir)

    # outputs
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    suffix = "_ALL" if len(scope_dirs) > 1 else f"_{scope_dirs[0].name}"
    out_csv = results_dir.parent / f"{results_dir.name}{suffix}_{ts}_vina_docking_scores_sorted.csv"
    live_csv = results_dir.parent / f"{results_dir.name}{suffix}_{ts}_LIVE.csv"
    log_path = results_dir.parent / f"{results_dir.name}{suffix}_{ts}_parse_progress.log"

    # collect targets
    targets = []
    if args.name == "*.pdbqt":
        for d in scope_dirs:
            targets.extend(d.rglob("*.pdbqt"))
    else:
        # exact filename match is MUCH faster on large trees
        for d in scope_dirs:
            for root, dirs, files in os.walk(d):
                # skip hidden dirs a bit
                dirs[:] = [x for x in dirs if not x.startswith(".")]
                if args.name in files:
                    targets.append(Path(root) / args.name)

    if not targets:
        print("⚠️ No target files found in selection.")
        sys.exit(2)

    print(f"\n📦 Found {len(targets)} files matching '{args.name}'. Parsing with {args.workers} workers…")

    # writers / locks
    csv_lock = threading.Lock()
    log_lock = threading.Lock()
    rows_buffer = []   # kept for final sort/write even if --live
    processed = 0
    hits = 0

    # init live CSV and log
    if args.live:
        with open(live_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            if args.fast:
                w.writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile"])
            else:
                w.writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile", "Run_Timestamp"])
        print(f"📝 LIVE CSV → {live_csv} (tail -f it)")

    if args.log:
        with open(log_path, "w", encoding="utf-8") as lf:
            lf.write(f"Parse started {ts}\nTotal targets: {len(targets)}\n\n")
        print(f"📄 Progress log → {log_path}")

    run_ts_str = ts

    def work(p: Path):
        # skip any receptor files if they sneak in
        if "receptor" in p.name.lower():
            return [], (str(p), False)
        receptor = p.parent.name
        ligand = p.parent.name if p.name == "out.pdbqt" else guess_ligand_name(p)
        # If the directory layout is .../<RECEPTOR>/<LIGAND>/out.pdbqt
        # receptor should be parent of ligand folder:
        try:
            receptor = p.parent.parent.name
            ligand = p.parent.name
        except Exception:
            pass

        poses = parse_vina_pdbqt(p)
        out_rows = []
        if poses:
            if args.fast:
                for pose_idx, aff in poses:
                    out_rows.append((receptor, ligand, pose_idx, aff, str(p.relative_to(results_dir))))
            else:
                for pose_idx, aff in poses:
                    out_rows.append((receptor, ligand, pose_idx, aff, str(p.relative_to(results_dir)), run_ts_str))
        return out_rows, (str(p), bool(poses))

    # parse concurrently and stream results/log heartbeats
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(work, p) for p in targets]
        for fut in as_completed(futs):
            out_rows, (ppath, has_results) = fut.result()

            # accumulate for final sort
            if out_rows:
                rows_buffer.extend(out_rows)

            # LIVE streaming
            if args.live and out_rows:
                with csv_lock:
                    with open(live_csv, "a", newline="", encoding="utf-8") as fh:
                        w = csv.writer(fh)
                        w.writerows(out_rows)

            # progress counters & logging
            processed += 1
            if has_results:
                hits += 1
                if args.log:
                    with log_lock:
                        with open(log_path, "a", encoding="utf-8") as lf:
                            lf.write(f"[FOUND] {ppath}  poses={len(out_rows)}\n")

            if args.log and processed % args.heartbeat == 0:
                with log_lock:
                    with open(log_path, "a", encoding="utf-8") as lf:
                        lf.write(f"[HEARTBEAT] processed={processed}  with_results={hits}\n")
                # also reflect basic progress to stdout
                pct = int(processed * 100 / max(1, len(targets)))
                sys.stdout.write(f"\r🔎 {processed}/{len(targets)} ({pct}%)  hits={hits}")
                sys.stdout.flush()

    # final heartbeat
    if args.log:
        with open(log_path, "a", encoding="utf-8") as lf:
            lf.write(f"[DONE] processed={processed} with_results={hits}\n")

    sys.stdout.write(f"\n✅ Parsed {processed} files; {hits} had results.\n")

    if args.no-sort:
        if args.live:
            print("ℹ️ Skipped final sorting (use the LIVE CSV).")
        else:
            print("ℹ️ Skipped final sorting and no LIVE CSV requested (nothing written).")
        return

    # sort & write final CSV
    def key(r):
        aff = r[3]
        try:
            return (r[0], float(aff))
        except Exception:
            return (r[0], 1e9)

    rows_buffer.sort(key=key)

    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if args.fast:
            w.writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile"])
        else:
            w.writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile", "Run_Timestamp"])
        w.writerows(rows_buffer)

    print(f"📊 Final sorted CSV → {out_csv}")
    if args.live:
        print("Tip: LIVE CSV remains for incremental inspection.")

if __name__ == "__main__":
    main()
