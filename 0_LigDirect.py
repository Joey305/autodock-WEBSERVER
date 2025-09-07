#!/usr/bin/env python3
import os
import sys
import csv
import argparse
from pathlib import Path
from rdkit import Chem

def sanitize(s: str) -> str:
    return "".join(c if (c.isalnum() or c in "-_.") else "_" for c in s).strip("_") or "mol"

def next_unique_path(base: Path) -> Path:
    if not base.exists():
        return base
    k = 1
    while True:
        cand = base.with_name(f"{base.name}_{k}")
        if not cand.exists():
            return cand
        k += 1

def list_sdf(cwd: Path) -> list[Path]:
    return sorted([p for p in cwd.glob("*.sdf") if p.is_file()], key=lambda x: x.name.lower())

def choose_index(items, title):
    print(f"\n{title}")
    for i, p in enumerate(items, start=1):
        print(f" [{i}] {p.name}")
    s = input("Enter index: ").strip()
    if not s.isdigit() or not (1 <= int(s) <= len(items)):
        print("❌ Invalid selection.", file=sys.stderr); sys.exit(2)
    return items[int(s)-1]

def split_sdf(input_sdf: Path, outdir: Path, prefix: str, pad: int, keep_names: bool) -> int:
    suppl = Chem.SDMolSupplier(str(input_sdf))
    if suppl is None:
        print(f"❌ Could not read SDF: {input_sdf}", file=sys.stderr)
        return 0

    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = outdir / "manifest.csv"
    n = 0

    with open(manifest_path, "w", newline="") as mf:
        w = csv.writer(mf)
        w.writerow(["index", "orig_name", "output_file"])
        for i, mol in enumerate(suppl, start=1):
            if mol is None:
                continue
            n += 1

            orig_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i:0{pad}d}"
            base_name = sanitize(orig_name) if keep_names else f"mol_{i:0{pad}d}"

            stem = f"{prefix}{i:0{pad}d}_{base_name}"
            out_path = outdir / f"{stem}.sdf"

            # ensure uniqueness if a collision somehow occurs
            if out_path.exists():
                out_path = next_unique_path(out_path)

            writer = Chem.SDWriter(str(out_path))
            writer.write(mol)
            writer.close()

            w.writerow([i, orig_name, out_path.name])

    return n

def build_parser():
    ap = argparse.ArgumentParser(
        description="Split a multi-record SDF into individual SDF files. Output dir named like 'Ligand_<TAG>'."
    )
    ap.add_argument("-i", "--input", help="Input multi-record SDF path (no default; pick interactively if omitted)")
    ap.add_argument("-n", "--name", help="Molecule tag used to form output dir 'Ligand_<TAG>' (e.g., CPD31)")
    ap.add_argument("-o", "--outdir", help="Explicit output directory (overrides --name). No extension.")
    ap.add_argument("--prefix", default="0_", help="Filename prefix (default: 0_)")
    ap.add_argument("--pad", type=int, default=4, help="Zero-padding width (default: 4 -> 0001)")
    ap.add_argument("--keep-names", action="store_true",
                    help="Use original _Name in filenames; else use mol_<index>.")
    ap.add_argument("--unique", action="store_true",
                    help="If output dir exists, append suffix (_1, _2, ...) to make a new dir.")
    return ap

def main():
    args = build_parser().parse_args()
    cwd = Path(".").resolve()

    # 1) Resolve input SDF
    if not args.input:
        candidates = list_sdf(cwd)
        if not candidates:
            print("❌ No .sdf files found in current directory. Provide --input.", file=sys.stderr)
            sys.exit(2)
        input_sdf = choose_index(candidates, "Select input SDF:")
    else:
        input_sdf = Path(args.input)
        if not input_sdf.is_file():
            print(f"❌ Input SDF not found: {input_sdf}", file=sys.stderr)
            sys.exit(2)

    # 2) Resolve output directory name
    if args.outdir:
        outdir_name = args.outdir
    else:
        mol_tag = args.name
        if not mol_tag:
            mol_tag = input("Enter molecule tag for output dir name 'Ligand_<TAG>': ").strip()
            if not mol_tag:
                print("❌ No tag provided.", file=sys.stderr); sys.exit(2)
        outdir_name = f"Ligand_{mol_tag}"

    # enforce "no extension" directory name
    outdir = cwd / Path(outdir_name).name  # strip any accidental path/extension
    if outdir.suffix:
        outdir = outdir.with_suffix("")  # remove an accidental extension

    # 3) Handle existing directory policy
    if outdir.exists():
        if args.unique:
            outdir = next_unique_path(outdir)
            print(f"ℹ️ Output dir exists; using unique path: {outdir.name}")
        else:
            print(f"⚠️ Output dir already exists: {outdir}\n"
                  f"    (Pass --unique to auto-make a new one, or use --outdir to specify another.)", file=sys.stderr)
            sys.exit(1)

    # 4) Split
    count = split_sdf(input_sdf, outdir, args.prefix, args.pad, args.keep_names)
    if count == 0:
        print("⚠️ No molecules were written (input may be empty/unreadable).", file=sys.stderr)
        sys.exit(1)

    print(f"✅ Split {count} molecules into '{outdir}/'")
    print(f"🗒️  Wrote manifest: {outdir/'manifest.csv'}")

if __name__ == "__main__":
    main()
