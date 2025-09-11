#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
autodock-proteinPREP — receptor prep for AutoDock Vina (interactive + headless)

FEATURES
- Reads PDB / mmCIF
- Lists only true HET groups (ligands/ions/waters/sugars)
- Optional chain removal
- AltLoc collapse (highest occupancy, tie: ' ' > 'A' > lexicographic)
- Converts to PDBQT using first available backend:
    1) Meeko (mk_prepare_receptor.py)
    2) AutoDockTools_py3 (installed module)
    3) AutoDockTools_py3 (local checkout ./AutoDockTools_py3)
    4) Legacy MGLTools (prepare_receptor4.py on PATH)
- FINAL OUTPUT: PDBQT-only (sweeps non-.pdbqt files from output dir by default)
- NEW: Interactive option "0" to open a GUI folder picker (Tk) that comes to the front

USAGE (INTERACTIVE)
    python 3a_PDB2PDBQTbatch.py
    # Follow the prompts (same heads-up display as before)
    # At the folder prompt, enter 0 to use a GUI folder browser (requires a desktop/display & Tk)

USAGE (HEADLESS / AUTOMATION)
    # Batch mode, remove all HETs, keep all chains, PDBQT-only (defaults)
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --headless

    # Batch mode, remove specific HET codes, remove chains B,C
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --headless \
        --mode batch \
        --remove-het HOH,EDO,SO4 \
        --remove-chains B,C

    # Batch mode, pick HETs by aggregated INDICES from a scan (after scanning)
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --headless \
        --mode batch \
        --remove-het-indices 1,3,5

    # Per-file headless via JSON config (filenames relative to --folder)
    # config.json:
    # {
    #   "7hg9.cif": {"remove_het": "all", "remove_chains": ["A","C"]},
    #   "8BB2_cleaned.pdb": {"remove_het": ["HOH","EDO"], "remove_chains": []}
    # }
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --headless \
        --mode per-file \
        --per-file-config config.json

    # Limit to a subset of files
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --files 7hg9.cif 8BB2_cleaned.pdb \
        --headless

    # Choose backend explicitly (auto/meeko/adt/mgl), change altloc policy,
    # and allow non-PDBQT artifacts (summary / clean PDBs)
    python 3a_PDB2PDBQTbatch.py \
        --folder Receptors \
        --headless \
        --backend meeko \
        --altloc collapse \
        --no-pdbqt-only \
        --write-summary \
        --keep-clean-pdb

ARGUMENTS (most useful)
  --folder FOLDER               Input folder containing .pdb/.cif
  --files F1 [F2 ...]           Optional subset of files to process (relative to folder)
  --headless                    Run without prompts (see defaults below)
  --mode {batch,per-file}       Headless mode selection (default: batch)
  --remove-het ALL|CODES        Comma list (e.g. HOH,EDO,SO4) or 'all' (default in headless: all)
  --remove-het-indices IDS      Comma list of indices from aggregated HET scan (batch only)
  --remove-chains IDS           Comma list of chain IDs (e.g. A,B). Default headless: none
  --per-file-config JSON        JSON mapping filename -> {remove_het, remove_chains}
  --backend {auto,meeko,adt,mgl}Backend preference (default: auto)
  --altloc {collapse,all}       AltLoc policy (default: collapse)
  --output-dir PATH             Output directory (default: <folder>_PDBQT_Converted)
  --pdbqt-only / --no-pdbqt-only  Keep only .pdbqt files (default: on)
  --write-summary               Write HET summary (only if not PDBQT-only)
  --keep-clean-pdb              Save cleaned PDBs (only if not PDBQT-only)
  --browse                      Force opening a GUI folder picker at start (interactive only)

HEADLESS DEFAULTS (when --headless is used)
  mode=batch, remove_het=all, remove_chains=none, pdbqt_only=True

REQUIREMENTS
  pip install gemmi
  (optional) pip install meeko
  (optional) pip install -e ./AutoDockTools_py3
  (optional for GUI) conda install tk  # or system Tk libraries; not available over pure SSH without a display
"""

import os
import sys
import json
import shutil
import subprocess
import argparse
import importlib.util
from collections import Counter, defaultdict

# =========================
# Dependencies
# =========================
try:
    import gemmi
except ImportError:
    print("❌ Requires 'gemmi'. Install with: pip install gemmi")
    sys.exit(1)

# =========================
# Arg parsing
# =========================
def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Prep PDB/mmCIF receptors and convert to PDBQT (interactive or headless)."
    )
    p.add_argument("--folder", type=str, help="Folder containing receptor structures")
    p.add_argument("--files", nargs="+", help="Subset of files to process (relative to --folder)")

    p.add_argument("--headless", action="store_true", help="Run without prompts")
    p.add_argument("--mode", choices=["batch", "per-file"], default="batch", help="Headless selection mode")

    het = p.add_mutually_exclusive_group()
    het.add_argument("--remove-het", type=str, help="Comma-separated 3-letter codes or 'all'")
    het.add_argument("--remove-het-indices", type=str, help="Comma-separated indices from aggregated HET scan (batch only)")

    p.add_argument("--remove-chains", type=str, help="Comma-separated chain IDs (e.g. A,B)")
    p.add_argument("--per-file-config", type=str, help="JSON config: filename->{remove_het, remove_chains}")

    p.add_argument("--backend", choices=["auto", "meeko", "adt", "mgl"], default="auto", help="Backend preference")
    p.add_argument("--altloc", choices=["collapse", "all"], default="collapse", help="AltLoc handling")

    p.add_argument("--output-dir", type=str, help="Output directory name")

    pdbqt = p.add_mutually_exclusive_group()
    pdbqt.add_argument("--pdbqt-only", dest="pdbqt_only", action="store_true", help="Keep only .pdbqt in final folder")
    pdbqt.add_argument("--no-pdbqt-only", dest="pdbqt_only", action="store_false", help="Allow non-.pdbqt artifacts")
    p.set_defaults(pdbqt_only=True)

    p.add_argument("--write-summary", action="store_true", help="Write HET summary (suppressed if --pdbqt-only)")
    p.add_argument("--keep-clean-pdb", action="store_true", help="Keep cleaned PDBs (suppressed if --pdbqt-only)")

    # NEW: force GUI folder picker
    p.add_argument("--browse", action="store_true", help="Open a GUI folder picker at start (interactive only)")

    return p.parse_args()

# =========================
# Backend resolution
# =========================
def has_module(modname: str) -> bool:
    return importlib.util.find_spec(modname) is not None

def resolve_receptor_preparer(prefer="auto"):
    """
    Returns (kind, base_cmd_list, extra_env_or_None)
    kind ∈ {'meeko','adt_module','adt_local','mgltools',None}
    """
    def _meeko():
        path = shutil.which("mk_prepare_receptor.py")
        return ("meeko", [path], None) if path else (None, None, None)

    def _adt_module():
        if has_module("AutoDockTools.Utilities24.prepare_receptor4"):
            return ("adt_module", [sys.executable, "-m", "AutoDockTools.Utilities24.prepare_receptor4"], None)
        return (None, None, None)

    def _adt_local():
        local_script = os.path.join("AutoDockTools_py3", "AutoDockTools", "Utilities24", "prepare_receptor4.py")
        if os.path.exists(local_script):
            env = os.environ.copy()
            adt_root = os.path.abspath("AutoDockTools_py3")
            env["PYTHONPATH"] = adt_root + (os.pathsep + env.get("PYTHONPATH",""))
            return ("adt_local", [sys.executable, local_script], env)
        return (None, None, None)

    def _mgl():
        path = shutil.which("prepare_receptor4.py")
        return ("mgltools", [path], None) if path else (None, None, None)

    order_map = {
        "auto":   [ _adt_module, _meeko, _adt_local, _mgl],
        "meeko":  [_meeko],
        "adt":    [_adt_module, _adt_local],
        "mgl":    [_mgl],
    }
    for fn in order_map.get(prefer, order_map["auto"]):
        kind, base, env = fn()
        if kind:
            return kind, base, env
    return None, None, None

# =========================
# HET / altloc / PDB writer
# =========================
STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "SEC","PYL","MSE"
}
STD_NT = {"A","C","G","T","U","I","DA","DC","DG","DT","DI","RA","RC","RG","RU"}
WATER_NAMES = {"HOH","WAT","H2O"}
COMMON_SUGARS = {"NAG","BMA","MAN","GAL","FUC","NDG"}

def is_true_het(res: 'gemmi.Residue', chain=None) -> bool:
    name = res.name.strip().upper()
    if name in WATER_NAMES:
        return True
    et = res.entity_type
    if et in (gemmi.EntityType.NonPolymer, gemmi.EntityType.Water):
        return True
    if name in STD_AA or name in STD_NT:
        return False
    if name in COMMON_SUGARS:
        return True
    if res.het_flag != ' ':
        return True
    return False

def iter_atoms_with_altloc_policy(res, policy="collapse"):
    if policy == "all":
        for a in res:
            yield a
        return
    groups = {}
    for a in res:
        key = a.name.strip()
        alt = getattr(a, "altloc", "") or ""
        occ = getattr(a, "occ", 1.0)
        prev = groups.get(key)
        if prev is None:
            groups[key] = (a, alt, occ)
        else:
            _, p_alt, p_occ = prev
            better = (occ > p_occ) or (occ == p_occ and (alt or " ") < (p_alt or " "))
            if better:
                groups[key] = (a, alt, occ)
    for a, _, _ in groups.values():
        yield a

def format_atom_name(atom_name: str) -> str:
    n = atom_name.strip()
    return n[:4] if len(n) >= 4 else f"{n:>4}"

def derive_element(atom) -> str:
    try:
        elem = atom.element.name.strip()
        if elem and elem != "X":
            return elem
    except Exception:
        pass
    letters = "".join([c for c in atom.name if c.isalpha()])
    if not letters:
        return ""
    if len(letters) >= 2 and letters[1].islower():
        return letters[:2].title()
    return letters[0].upper()

def pdb_write_manual(structure: 'gemmi.Structure', out_pdb_path: str, altloc_policy="collapse"):
    serial = 1
    lines = []
    wrote_any = False
    for model in structure:
        for chain in model:
            chain_id = (chain.name[:1] if chain.name else "A")
            chain_wrote = False
            last_resname, last_resseq, last_icode = "UNK", 0, " "
            for res in chain:
                resname = res.name.strip().upper()
                try:
                    resseq = res.seqid.num
                    icode = res.seqid.icode if res.seqid.icode and res.seqid.icode != "\x00" else " "
                except Exception:
                    resseq, icode = 0, " "
                record = "HETATM" if is_true_het(res, chain) else "ATOM  "
                for atom in iter_atoms_with_altloc_policy(res, policy=altloc_policy):
                    x, y, z = atom.pos.x, atom.pos.y, atom.pos.z
                    occ = getattr(atom, "occ", 1.00)
                    bfac = getattr(atom, "b_iso", 0.00)
                    altloc = atom.altloc if getattr(atom, "altloc", "") else " "
                    atname = format_atom_name(atom.name)
                    element = derive_element(atom)
                    line = (
                        f"{record:<6}{serial:>5} "
                        f"{atname:<4}{altloc:1}"
                        f"{resname:>3} "
                        f"{chain_id:1}{resseq:>4}{icode:1}"
                        f"   "
                        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                        f"{occ:>6.2f}{bfac:>6.2f}"
                        f"          {element:>2}  "
                    )
                    lines.append(line)
                    serial += 1
                    wrote_any = True
                    chain_wrote = True
                    last_resname, last_resseq, last_icode = resname, resseq, icode
            if chain_wrote:
                ter = f"{'TER':<6}{serial:>5}      {last_resname:>3} {chain_id:1}{last_resseq:>4}{last_icode:1}"
                lines.append(ter)
                serial += 1
    lines.append("END")
    with open(out_pdb_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    if not wrote_any:
        raise RuntimeError("No atoms written — structure empty after filtering.")

def write_pdb_compat(structure: 'gemmi.Structure', out_pdb_path: str, altloc_policy="collapse"):
    if hasattr(gemmi, "write_minimal_pdb") and altloc_policy == "collapse":
        try:
            gemmi.write_minimal_pdb(structure, out_pdb_path)
            return
        except Exception:
            pass
    pdb_write_manual(structure, out_pdb_path, altloc_policy=altloc_policy)

# =========================
# IO / scanning helpers
# =========================
def read_structure(path: str) -> 'gemmi.Structure':
    return gemmi.read_structure(path)

def is_cif_name(filename: str) -> bool:
    return os.path.splitext(filename.lower())[1] in {".cif", ".mmcif"}

def scan_file_for_hets_and_chains(path: str):
    hets = Counter()
    chains = set()
    st = read_structure(path)
    for model in st:
        for chain in model:
            chains.add(chain.name)
            for res in chain:
                if is_true_het(res, chain):
                    hets[res.name.strip().upper()] += 1
    return hets, chains

def scan_all(files, folder):
    het_counter = Counter()
    het_presence = defaultdict(set)
    chains_per_file = defaultdict(set)
    for fname in files:
        fpath = os.path.join(folder, fname)
        try:
            h, chs = scan_file_for_hets_and_chains(fpath)
            for k, v in h.items():
                het_counter[k] += v
                het_presence[k].add(fname)
            chains_per_file[fname] = chs
        except Exception as e:
            print(f"⚠️ Skipping {fname}: {e}")
    return het_counter, het_presence, chains_per_file

def clean_structure_to_pdb(in_path, out_pdb_path, residues_to_drop, chains_to_drop, altloc_policy="collapse"):
    st = gemmi.read_structure(in_path)
    new_st = gemmi.Structure()
    new_st.cell = st.cell
    new_st.spacegroup_hm = st.spacegroup_hm
    new_st.name = st.name

    res_drop = {r.upper() for r in residues_to_drop}
    chain_drop = set(chains_to_drop)

    for mi, model in enumerate(st):
        new_model = gemmi.Model(str(mi + 1))
        for chain in model:
            if chain.name in chain_drop:
                continue
            new_chain = gemmi.Chain(chain.name)
            for res in chain:
                if is_true_het(res, chain) and res.name.strip().upper() in res_drop:
                    continue
                new_chain.add_residue(res)
            if len(new_chain) > 0:
                new_model.add_chain(new_chain)
        if len(new_model) > 0:
            new_st.add_model(new_model)

    write_pdb_compat(new_st, out_pdb_path, altloc_policy=altloc_policy)

# =========================
# GUI folder picker (Tk)
# =========================
def choose_folder_via_tk():
    """Open a Tk folder picker, bring it to front, return selected path or ''."""
    try:
        import tkinter as tk
        from tkinter import filedialog
    except Exception as e:
        print(f"⚠️ GUI picker unavailable (tkinter not found): {e}")
        return ""
    try:
        root = tk.Tk()
        root.withdraw()
        # Bring dialog to front on all platforms
        try:
            root.attributes("-topmost", True)
            root.update()
        except Exception:
            pass
        # Show
        folder = filedialog.askdirectory(title="Select receptor folder (contains .pdb / .cif)")
        try:
            root.destroy()
        except Exception:
            pass
        return folder or ""
    except Exception as e:
        print(f"⚠️ GUI picker failed (no display or Tk error): {e}")
        return ""

# =========================
# Utilities
# =========================
def parse_csv_list(s):
    return [x.strip() for x in s.split(",") if x.strip()]

def load_per_file_config(path):
    with open(path, "r") as f:
        cfg = json.load(f)
    # Normalize: codes to upper, allow 'all'
    norm = {}
    for k, v in cfg.items():
        rm_het = v.get("remove_het", [])
        if isinstance(rm_het, str):
            rm_het_norm = rm_het.strip().lower()
            if rm_het_norm == "all":
                rm_het = "all"
            else:
                rm_het = [x.strip().upper() for x in parse_csv_list(rm_het)]
        else:
            rm_het = [str(x).strip().upper() for x in rm_het]
        rm_ch = [str(x).strip() for x in v.get("remove_chains", [])]
        norm[k] = {"remove_het": rm_het, "remove_chains": rm_ch}
    return norm

# =========================
# Main
# =========================
def main():
    args = parse_args()

    # Choose folder
    if args.headless:
        if not args.folder:
            print("❌ In headless mode you must pass --folder.")
            sys.exit(1)
        folder = args.folder
    else:
        # Interactive path
        folder = args.folder
        if args.browse and not folder:
            picked = choose_folder_via_tk()
            if picked:
                folder = picked

        if not folder:
            print("\n📁 Available folders in current directory:")
            folders = [f for f in os.listdir() if os.path.isdir(f)]
            for i, folder_name in enumerate(folders, 1):
                print(f"{i}. {folder_name}")
            print("0. 🖱️  Browse for a folder (GUI)")

            choice = input("🔍 Enter the number of the folder containing receptor structures (.pdb / .cif): ").strip()
            if choice == "0":
                picked = choose_folder_via_tk()
                if picked:
                    folder = picked
                else:
                    print("⚠️ GUI folder selection unavailable. Please choose from the list.")
            if not folder:
                try:
                    folder_index = int(choice) - 1
                    folder = folders[folder_index]
                except (ValueError, IndexError):
                    raise ValueError("❌ Invalid folder selection.")

    # Collect files
    valid_exts = {".pdb", ".ent", ".cif", ".mmcif"}
    try:
        file_candidates = os.listdir(folder)
    except FileNotFoundError:
        print(f"❌ Folder not found: {folder}")
        sys.exit(1)

    all_files = sorted(
        f for f in file_candidates
        if os.path.isfile(os.path.join(folder, f))
        and os.path.splitext(f.lower())[1] in valid_exts
    )
    if not all_files:
        print("❌ No .pdb / .cif / .mmcif files found in the selected folder.")
        sys.exit(1)

    if args.files:
        # Validate subset
        missing = [f for f in args.files if f not in all_files]
        if missing:
            print(f"❌ These files are not in {folder}: {', '.join(missing)}")
            sys.exit(1)
        files = sorted(args.files)
    else:
        files = all_files

    # Interactive vs headless mode
    interactive = not args.headless

    # Head-up display: batch vs per-file (interactive path)
    if interactive and len(files) > 1:
        ans = input("\n🗂️ Multiple structures detected. Batch process with the SAME choices for all? [Y/n]: ").strip().lower()
        mode = "batch" if ans in ("", "y", "yes") else "per-file"
    else:
        mode = args.mode  # headless default 'batch'

    # Backend resolve
    prefer = args.backend if args.backend else "auto"
    prep_kind, prep_base, prep_env = resolve_receptor_preparer(prefer=prefer)
    if not prep_kind:
        raise RuntimeError(
            "No receptor preparer found. Install Meeko (`pip install meeko`) or ADT_py3 "
            "(`pip install -e ./AutoDockTools_py3` or `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3`) "
            "or ensure MGLTools' prepare_receptor4.py is on PATH."
        )

    # Output dirs
    output_dir = args.output_dir if args.output_dir else (folder + "_PDBQT_Converted")
    os.makedirs(output_dir, exist_ok=True)
    temp_cleaned_dir = ".temp_cleaned_pdbs"
    os.makedirs(temp_cleaned_dir, exist_ok=True)

    # Banner
    if args.pdbqt_only:
        print("\n🧾 PDBQT-only mode: final output folder will contain ONLY .pdbqt files.\n")

    print(f"🔧 Using receptor preparer: {prep_kind}")

    # Scan for HETs/chains
    if mode == "batch":
        print(f"\n📦 Found {len(files)} structure file(s). Scanning for HET ligands/waters and chains...")
        het_counter, het_presence, chains_per_file = scan_all(files, folder)
        het_names_sorted = sorted(het_counter.items(), key=lambda x: (-x[1], x[0]))

        # Display HUD in both interactive and headless (nice logs)
        print("\n🔎 HET (ligands/ions/waters) detected across all files:")
        if not het_names_sorted:
            print("  (none)")
        else:
            for idx, (name, count) in enumerate(het_names_sorted, 1):
                files_here = ", ".join(sorted(het_presence[name]))
                print(f"  {idx:>3}. {name:<6} | count: {count:<4} | files: {files_here}")

        print("\n🔗 Chains detected (per file):")
        for fname in files:
            chs = sorted(chains_per_file.get(fname, []))
            print(f"  • {fname}: {', '.join(chs) if chs else '(none)'}")

        residues_to_remove = set()
        if interactive and het_names_sorted:
            sel = input("\n🧽 Select HET indices to remove (e.g. 1,2,5) or 'all' to remove all shown. Leave blank to skip: ").strip().lower()
            if sel == "all":
                residues_to_remove = {name for name, _ in het_names_sorted}
            elif sel:
                try:
                    chosen = {int(x.strip()) for x in sel.split(",") if x.strip()}
                    for i in chosen:
                        if 1 <= i <= len(het_names_sorted):
                            residues_to_remove.add(het_names_sorted[i-1][0])
                        else:
                            print(f"  ⚠️ Ignoring out-of-range index: {i}")
                except ValueError:
                    print("  ⚠️ Invalid index list, skipping automatic HET removal.")
        else:
            # Headless defaults & flags
            if args.remove_het_indices:
                # Map indices to names from aggregated list
                indices = []
                try:
                    indices = [int(x.strip()) for x in parse_csv_list(args.remove_het_indices)]
                except ValueError:
                    print("❌ --remove-het-indices must be integers")
                    sys.exit(1)
                for i in indices:
                    if 1 <= i <= len(het_names_sorted):
                        residues_to_remove.add(het_names_sorted[i-1][0])
                    else:
                        print(f"  ⚠️ Ignoring out-of-range index: {i}")
            elif args.remove_het:
                s = args.remove_het.strip().lower()
                if s == "all":
                    residues_to_remove = {name for name, _ in het_names_sorted} if het_names_sorted else set()
                else:
                    residues_to_remove = {x.upper() for x in parse_csv_list(args.remove_het)}
            else:
                # headless default: remove ALL HETs
                residues_to_remove = {name for name, _ in het_names_sorted} if het_names_sorted else set()

        if interactive:
            extra_res = input("➕ Enter any EXTRA 3-letter residue names to remove (comma-separated, optional): ").strip()
            if extra_res:
                residues_to_remove |= {r.strip().upper() for r in extra_res.split(",") if r.strip()}
            if residues_to_remove:
                print("🧼 Will remove HET residues (for all files):", ", ".join(sorted(residues_to_remove)))
            else:
                print("ℹ️ No HET residue names selected for removal.")
            chain_resp = input("🧩 Remove specific chains for ALL files? Enter chain IDs (comma-separated, leave blank to skip): ").strip()
            chains_to_remove_global = {c.strip() for c in chain_resp.split(",") if c.strip()} if chain_resp else set()
            if chains_to_remove_global:
                print("✂️ Will remove chains (for all files):", ", ".join(sorted(chains_to_remove_global)))
            else:
                print("ℹ️ No chains selected for removal.")
        else:
            chains_to_remove_global = set(parse_csv_list(args.remove_chains)) if args.remove_chains else set()
            if residues_to_remove:
                print("🧼 Will remove HET residues (headless / all files):", ", ".join(sorted(residues_to_remove)))
            if chains_to_remove_global:
                print("✂️ Will remove chains (headless / all files):", ", ".join(sorted(chains_to_remove_global)))

        per_file_choices = {fname: (set(residues_to_remove), set(chains_to_remove_global)) for fname in files}

        # Summary file (only if allowed)
        if (not args.pdbqt_only) and args.write_summary:
            summary_path = os.path.join(output_dir, "HETATM_summary.txt")
            with open(summary_path, "w") as sf:
                sf.write(f"Folder: {folder}\n")
                sf.write("HET (ligand/ion/water) residue counts (index | name | count | files)\n")
                if het_names_sorted:
                    for idx, (name, count) in enumerate(het_names_sorted, 1):
                        files_str = ", ".join(sorted(het_presence[name]))
                        sf.write(f"{idx:>3} | {name:<6} | {count:<4} | {files_str}\n")
                else:
                    sf.write("(none)\n")
            print(f"\n🧾 Wrote HET summary → {summary_path}")
        elif args.pdbqt_only and args.write_summary:
            print("ℹ️ --write-summary requested, but --pdbqt-only is active; summary will be swept. Use --no-pdbqt-only to keep it.")

    else:
        # per-file mode
        if interactive:
            print(f"\n📦 Found {len(files)} structure file(s). Per-file selection enabled.\n")
        # Load per-file config (headless)
        cfg = load_per_file_config(args.per_file_config) if (args.headless and args.per_file_config) else {}

        per_file_choices = {}
        for fname in files:
            fpath = os.path.join(folder, fname)
            if interactive:
                print(f"── {fname} ─────────────────────────────────────────────")
            try:
                hets, chains = scan_file_for_hets_and_chains(fpath)
            except Exception as e:
                print(f"⚠️ Skipping (scan failed): {e}")
                continue

            het_sorted = sorted(hets.items(), key=lambda x: (-x[1], x[0]))
            if interactive:
                print("🔎 HET (ligands/ions/waters) in this file:")
                if not het_sorted:
                    print("  (none)")
                else:
                    for idx, (name, count) in enumerate(het_sorted, 1):
                        print(f"  {idx:>3}. {name:<6} | count: {count:<4}")
                print("🔗 Chains in this file:", ", ".join(sorted(chains)) if chains else "(none)")

            if interactive:
                res_remove = set()
                if het_sorted:
                    sel = input("🧽 Select HET indices to remove for THIS file (1,2,...) or 'all' or blank to skip: ").strip().lower()
                    if sel == "all":
                        res_remove = {name for name, _ in het_sorted}
                    elif sel:
                        try:
                            chosen = {int(x.strip()) for x in sel.split(",") if x.strip()}
                            for i in chosen:
                                if 1 <= i <= len(het_sorted):
                                    res_remove.add(het_sorted[i-1][0])
                                else:
                                    print(f"  ⚠️ Ignoring out-of-range index: {i}")
                        except ValueError:
                            print("  ⚠️ Invalid index list, skipping automatic HET removal.")
                extra_res = input("➕ Extra 3-letter residue names to remove (optional): ").strip()
                if extra_res:
                    res_remove |= {r.strip().upper() for r in extra_res.split(",") if r.strip()}
                if res_remove:
                    print("🧼 Will remove HET residues (this file):", ", ".join(sorted(res_remove)))
                else:
                    print("ℹ️ No HET residue names selected (this file).")
                chain_resp = input("🧩 Chains to remove for THIS file (comma-separated, blank to keep all): ").strip()
                chains_to_remove = {c.strip() for c in chain_resp.split(",") if c.strip()} if chain_resp else set()
                if chains_to_remove:
                    print("✂️ Will remove chains (this file):", ", ".join(sorted(chains_to_remove)))
                else:
                    print("ℹ️ Keeping all chains (this file).")
            else:
                # headless per-file: use config or defaults (all HETs, no chains)
                entry = cfg.get(fname, {})
                rm_het = entry.get("remove_het", "all")
                if isinstance(rm_het, str) and rm_het.strip().lower() == "all":
                    res_remove = {name for name, _ in het_sorted}
                else:
                    res_remove = {str(x).strip().upper() for x in rm_het}
                chains_to_remove = {str(x).strip() for x in entry.get("remove_chains", [])}

                # Log choice
                if res_remove:
                    print(f"🧼 [{fname}] remove HET: {', '.join(sorted(res_remove))}")
                if chains_to_remove:
                    print(f"✂️ [{fname}] remove chains: {', '.join(sorted(chains_to_remove))}")

            per_file_choices[fname] = (res_remove, chains_to_remove)
            if interactive:
                print()

    # PROCESS
    success = 0
    fail = 0
    for fname in files:
        in_path = os.path.join(folder, fname)
        base = os.path.splitext(fname)[0]
        cleaned_pdb = os.path.join(".temp_cleaned_pdbs", base + ".clean.pdb")
        residues_to_remove, chains_to_remove = per_file_choices.get(fname, (set(), set()))

        try:
            clean_structure_to_pdb(in_path, cleaned_pdb, residues_to_remove, chains_to_remove, altloc_policy=args.altloc)
            print(f"\n🧹 Cleaning {fname} → {os.path.basename(cleaned_pdb)}")
        except Exception as e:
            print(f"❌ Failed cleaning {fname}: {e}")
            fail += 1
            continue

        out_pdbqt = os.path.join(output_dir, base + ".converted.pdbqt")
        # Build conversion command by backend
        if prep_kind == "meeko":
            # Meeko CLI expects --read_pdb/--write_pdbqt (or -i/--read_with_prody)
            cmd = list(prep_base) + [
                "--read_pdb", cleaned_pdb,
                "--write_pdbqt", out_pdbqt
            ]
            # optional: you can add box or other params here if desired
        elif prep_kind in ("adt_module", "adt_local", "mgltools"):
            # ADT/MGL takes -r/-o and usually benefits from -A checkhydrogens
            cmd = list(prep_base) + [
                "-r", cleaned_pdb,
                "-o", out_pdbqt,
                "-A", "checkhydrogens"
            ]
        else:
            raise RuntimeError(f"Unsupported preparer kind: {prep_kind}")

        env = prep_env if prep_env is not None else os.environ.copy()

        try:
            print(f"⚙️ Converting {fname} → {os.path.basename(out_pdbqt)}")
            subprocess.run(cmd, check=True, env=env)
            print(f"✅ Saved: {out_pdbqt}")
            # Optional artifacts if allowed
            if (not args.pdbqt_only) and args.keep_clean_pdb and (not is_cif_name(fname)):
                shutil.copy2(cleaned_pdb, os.path.join(output_dir, base + ".clean.pdb"))
            success += 1
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed conversion on {fname}: {e}")
            fail += 1

    # Cleanup temp
    try:
        shutil.rmtree(".temp_cleaned_pdbs")
    except Exception:
        pass

    # Final sweep (PDBQT-only)
    if args.pdbqt_only:
        for g in os.listdir(output_dir):
            if not g.lower().endswith(".pdbqt"):
                try:
                    os.remove(os.path.join(output_dir, g))
                except Exception:
                    pass

    print("\n🎉 Done.")
    print(f"   ✅ Converted: {success}")
    print(f"   ❌ Failed:    {fail}")
    print(f"📦 Output folder: {output_dir}")

if __name__ == "__main__":
    main()