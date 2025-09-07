#!/usr/bin/env python3
import os, sys, re, csv, argparse, threading
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# --------- patterns ----------
VINA_RE = re.compile(r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)")
OK_LINE = re.compile(r"^[\s✅\-\+]?\s*(.+?)\s*→\s*([A-Za-z0-9._-]+)\s*$")  # ligand  → receptor.pdbqt

# --------- helpers ----------
def parse_vina_pdbqt(p: Path):
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

def choose(prompt: str, items):
    print(f"\n{prompt}")
    for i, item in enumerate(items):
        name = item.name if hasattr(item, "name") else str(item)
        print(f" [{i}] {name}")
    while True:
        s = input("Enter index: ").strip()
        try:
            i = int(s)
            if 0 <= i < len(items):
                return items[i]
        except ValueError:
            pass
        print("Invalid choice. Try again.")

def prompt_bool(question: str, default: bool) -> bool:
    d = "Y/n" if default else "y/N"
    while True:
        s = input(f"{question} [{d}]: ").strip().lower()
        if s == "" : return default
        if s in ("y","yes"): return True
        if s in ("n","no"):  return False
        print("Please answer y or n.")

def prompt_int(question: str, default: int) -> int:
    while True:
        s = input(f"{question} [{default}]: ").strip()
        if s == "": return default
        try:
            v = int(s)
            if v > 0: return v
        except ValueError:
            pass
        print("Please enter a positive integer.")

def choose_scope(results_dir: Path):
    subs = [d for d in sorted(results_dir.iterdir()) if d.is_dir()]
    if not subs:
        print(f"⚠️ No subfolders under {results_dir.name}. Using it directly.")
        return [results_dir]
    print(f"\nSelect subfolder(s) under {results_dir.name}:")
    print(" [0] All subfolders")
    for i, d in enumerate(subs, 1):
        print(f" [{i}] {d.name}")
    while True:
        s = input("Enter index or 0 for ALL: ").strip()
        if s == "0": return subs
        try:
            i = int(s)
            if 1 <= i <= len(subs):
                return [subs[i-1]]
        except ValueError:
            pass
        print("Invalid choice. Try again.")

def find_candidates():
    return [d for d in sorted(Path(".").glob("Docking_Results_*")) if d.is_dir()]

# --------- main ----------
def main():
    # flags remain available for non-interactive use
    ap = argparse.ArgumentParser(
        description="Parse Vina results into CSV. Interactive prompts when flags are omitted."
    )
    ap.add_argument("--from-log", dest="from_log", help="Path to run_log_*.txt (skips directory crawl).")
    ap.add_argument("--live", action="store_true", help="Stream rows to LIVE CSV (tail -f).")
    ap.add_argument("--log", action="store_true", help="Write parse progress log.")
    ap.add_argument("--no-sort", dest="no_sort", action="store_true", help="Skip final sorted CSV.")
    ap.add_argument("--workers", type=int, default=None, help="Worker threads (default ≈ 2× CPU cores, capped at 32).")
    ap.add_argument("--heartbeat", type=int, default=None, help="Heartbeat interval for progress log.")
    args = ap.parse_args()

    # Choose results dir
    candidates = find_candidates()
    if not candidates:
        print("❌ No Docking_Results_* directories in current folder."); sys.exit(1)
    results_dir = choose("Select a docking results folder:", candidates)

    # Discover logs & prompt if --from-log not given
    if not args.from_log:
        # look for likely logs
        logs = sorted(Path(".").glob("run_log_*")) + sorted(Path(".").glob("*.log"))
        use_log = False
        if logs:
            use_log = prompt_bool("Use a run log to drive parsing (faster, no directory crawl)?", True)
        else:
            print("ℹ️ No run logs detected in this folder.")
        if use_log and logs:
            # add manual path option
            choice_items = logs + ["<Enter a different path>"]
            sel = choose("Select a run log (or choose the last option to type a path):", choice_items)
            if isinstance(sel, str):
                manual = input("Enter path to run log: ").strip()
                args.from_log = manual
            else:
                args.from_log = str(sel.resolve())
        elif use_log:
            manual = input("Enter path to run log: ").strip()
            args.from_log = manual

    # Prompt for live/log/no-sort/workers/heartbeat if not set
    default_workers = min(32, (os.cpu_count() or 8) * 2)
    default_heartbeat = 5000

    live = args.live if args.live else prompt_bool("Write LIVE CSV while parsing?", True)
    write_parse_log = args.log if args.log else prompt_bool("Write a parse progress log?", True)
    no_sort = args.no_sort if args.no_sort else (not prompt_bool("Also create a final sorted CSV when done?", True))
    workers = args.workers if args.workers else prompt_int("How many worker threads?", default_workers)
    heartbeat = args.heartbeat if args.heartbeat else prompt_int("Heartbeat interval (files per log heartbeat)?", default_heartbeat)

    # Output names
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    suffix = "_ALL"  # log-driven spans receptors; crawl mode may too
    out_base = f"{results_dir.name}{suffix}_{ts}"
    out_csv  = results_dir.parent / f"{out_base}_vina_docking_scores_sorted.csv"
    live_csv = results_dir.parent / f"{out_base}_LIVE.csv"
    parse_log = results_dir.parent / f"{out_base}_parse_progress.log"

    # Build target list
    targets = []

    if args.from_log:
        log_path = Path(args.from_log).expanduser().resolve()
        if not log_path.is_file():
            print(f"❌ Log not found: {log_path}"); sys.exit(1)

        ok_pairs = []
        with open(log_path, "r", encoding="utf-8", errors="ignore") as lf:
            for line in lf:
                if "✅" not in line and "->" not in line and "→" not in line:
                    continue
                # Normalize arrows
                line = line.replace("->", "→")
                m = OK_LINE.match(line.strip())
                if not m:
                    continue
                ligand_dir = m.group(1).strip()
                receptor_file = m.group(2).strip()
                receptor_dir = receptor_file.replace(".pdbqt", "").replace(".pdb", "").replace(".mol2", "")
                p = results_dir / receptor_dir / ligand_dir / "out.pdbqt"
                ok_pairs.append((receptor_dir, ligand_dir, p))

        if not ok_pairs:
            print("⚠️ No successful entries found in the log (no '✅ … → receptor').")
            sys.exit(2)

        targets = ok_pairs
        print(f"\n🗒️  Log-driven mode: {len(targets)} completed pairs in {log_path.name}")

    else:
        # Crawl mode (slower) — still keep it around as a fallback
        scope_dirs = choose_scope(results_dir)
        for d in scope_dirs:
            for root, dirs, files in os.walk(d):
                dirs[:] = [x for x in dirs if not x.startswith(".")]
                if "out.pdbqt" in files:
                    p = Path(root) / "out.pdbqt"
                    receptor = p.parent.parent.name
                    ligand   = p.parent.name
                    targets.append((receptor, ligand, p))
        if not targets:
            print("⚠️ No out.pdbqt files found."); sys.exit(2)
        print(f"\n📦 Found {len(targets)} out.pdbqt files (crawl mode).")

    # Init LIVE CSV / log
    if live:
        with open(live_csv, "w", newline="", encoding="utf-8") as fh:
            csv.writer(fh).writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity"])
        print(f"📝 LIVE CSV → {live_csv}  (tail -f)")

    if write_parse_log:
        with open(parse_log, "w", encoding="utf-8") as lf:
            lf.write(f"Parse started {ts}\nTotal targets: {len(targets)}\n\n")
        print(f"📄 Progress log → {parse_log}")

    csv_lock = threading.Lock()
    log_lock = threading.Lock()
    processed = 0
    hits = 0
    rows_buffer = []

    def work(rec_lig_path):
        receptor, ligand, p = rec_lig_path
        if not p.is_file():
            return [], (str(p), False)
        poses = parse_vina_pdbqt(p)
        rows = []
        for pose_idx, aff in poses:
            rows.append((receptor, ligand, pose_idx, aff))
        return rows, (str(p), bool(poses))


    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(work, item) for item in targets]
        for fut in as_completed(futs):
            out_rows, (ppath, has_results) = fut.result()
            processed += 1
            if has_results:
                hits += 1

            if out_rows:
                rows_buffer.extend(out_rows)
                if live:
                    with csv_lock:
                        with open(live_csv, "a", newline="", encoding="utf-8") as fh:
                            csv.writer(fh).writerows(out_rows)
                if write_parse_log:
                    with log_lock:
                        with open(parse_log, "a", encoding="utf-8") as lf:
                            lf.write(f"[FOUND] {ppath}  poses={len(out_rows)}\n")

            if write_parse_log and processed % max(1, heartbeat) == 0:
                pct = int(processed * 100 / max(1, len(targets)))
                with log_lock:
                    with open(parse_log, "a", encoding="utf-8") as lf:
                        lf.write(f"[HEARTBEAT] processed={processed}/{len(targets)} ({pct}%)  hits={hits}\n")
                sys.stdout.write(f"\r🔎 {processed}/{len(targets)} ({pct}%)  hits={hits}")
                sys.stdout.flush()

    sys.stdout.write(f"\n✅ Parsed {processed} completed pairs; {hits} files had results.\n")

    if no_sort:
        if live:
            print("ℹ️ Skipped final sorting (use the LIVE CSV).")
        else:
            print("ℹ️ Skipped final sorting and no LIVE CSV requested (nothing written).")
        return

    # Final sort & write (by receptor, then affinity asc)
    rows_buffer.sort(key=lambda r: (r[0], float(r[3]) if isinstance(r[3], (int,float)) else 1e9))
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Receptor", "Ligand", "Pose", "Binding_Affinity", "OutFile"])
        w.writerows(rows_buffer)
    print(f"📊 Final sorted CSV → {out_csv}")
    if live:
        print("Tip: LIVE CSV remains for incremental inspection.")

if __name__ == "__main__":
    main()
