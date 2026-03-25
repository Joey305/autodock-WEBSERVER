#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, csv, shutil, subprocess
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ---- Headless PyMOL launch ----
import pymol as _pymol
try:
    _pymol.pymol_argv = ['pymol', '-cq']  # -c no GUI, -q quiet
    _pymol.finish_launching()
except Exception:
    pass
from pymol import cmd

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Geometry import Point3D

# ==================== regex & name helpers ====================
POSE_PAT   = re.compile(r"__pose(\d+)$")
BUCKET_PAT = re.compile(r"(CPD(?:XL|XH)?\d+)", re.IGNORECASE)
POSE_DIR_PAT = re.compile(r"(.+?)_pose\d+$", re.IGNORECASE)

def base_ligand_from_outfile(outfile: str) -> str | None:
    """
    Extract base ligand ID from OutFile path.
    Example:
      kif20A_Allo/UMF-604_pose17/out.pdbqt → UMF-604
    """
    if not outfile:
        return None
    try:
        lig_dir = Path(outfile).parent.name
        m = POSE_DIR_PAT.match(lig_dir)
        return m.group(1) if m else lig_dir
    except Exception:
        return None

def strip_pose_suffix(s: str) -> str:
    return POSE_PAT.sub("", s or "")

def extract_pose_num(s: str) -> str | None:
    m = POSE_PAT.search(s or "")
    return m.group(1) if m else None

def list_and_select(prompt, options):
    print(f"\n📁 {prompt}")
    for i, opt in enumerate(options, 1):
        print(f"{i}. {opt}")
    idx = int(input("🔢 Select: ")) - 1
    return options[idx]

def safe_float(x, default=1e9):
    try:
        return float(x)
    except Exception:
        return default

def sanitize_pymol_name(name: str) -> str:
    cleaned = re.sub(r'[^A-Za-z0-9_]', '_', name or "")
    if not cleaned or not cleaned[0].isalpha():
        cleaned = 'obj_' + cleaned
    return cleaned[:80]

def norm_variants(s: str):
    s = (s or "").strip()
    xs = {s, s.replace("-", "_"), s.replace("_", "-")}
    xs |= {x.lower() for x in list(xs)}
    return list(xs)

def norm_base(s: str):
    return re.sub(r'[^a-z0-9]', '', (s or '').lower())

# -------------------- CSV inspection & grouping --------------------

def detect_csv_style(rows):
    if not rows:
        return "simple"
    usable = 0
    for r in rows:
        out = (r.get("OutFile") or "").strip()
        rr  = (r.get("ResultsRoot") or "").strip()
        rd  = (r.get("ReceptorDir") or "").strip()
        ld  = (r.get("LigandDir") or "").strip()
        if out or (rr and rd and ld):
            usable += 1
    return "provenance" if usable >= max(1, int(0.1 * len(rows))) else "simple"

def load_csv_grouped(csv_path: Path, top_n: int, mode: str):
    print(f"\n📄 Reading top poses from: {csv_path}")

    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        raise SystemExit("❌ CSV is empty.")

    # --------------------------------------------------
    # Required columns
    # --------------------------------------------------
    required = ["Receptor", "Ligand", "Pose"]
    missing = [k for k in required if k not in rows[0]]
    if missing:
        raise SystemExit(f"❌ CSV missing required columns: {missing}")

    # --------------------------------------------------
    # Detect binding affinity column
    # --------------------------------------------------
    bind_col = (
        "Binding_Affinity"
        if "Binding_Affinity" in rows[0]
        else (
            "Binding_Affinity_kcal_per_mol"
            if "Binding_Affinity_kcal_per_mol" in rows[0]
            else None
        )
    )

    if bind_col:
        rows.sort(key=lambda r: safe_float(r.get(bind_col, "1e9")))

    # --------------------------------------------------
    # Detect CSV style (simple vs provenance)
    # --------------------------------------------------
    style = detect_csv_style(rows)

    # --------------------------------------------------
    # 🔥 DERIVE BaseLigand ONCE, BEFORE GROUPING 🔥
    # --------------------------------------------------
    for row in rows:
        base = base_ligand_from_outfile(row.get("OutFile", ""))
        if base:
            row["BaseLigand"] = base
        else:
            # fallback for legacy/simple CSVs
            row["BaseLigand"] = strip_pose_suffix(row["Ligand"])

    # --------------------------------------------------
    # Group rows
    # --------------------------------------------------
    grouped = defaultdict(list)

    if mode == "per_ligand":
        # Top N poses PER ligand (multi-state objects)
        per_combo_counts = defaultdict(lambda: defaultdict(int))  # receptor → base_ligand → count

        for row in rows:
            rec  = row["Receptor"]
            base = row["BaseLigand"]

            if per_combo_counts[rec][base] < top_n:
                grouped[rec].append(row)
                per_combo_counts[rec][base] += 1

    elif mode == "per_receptor":
        # Top N ligands PER receptor (best pose only)
        seen_bases = defaultdict(set)

        for row in rows:
            rec  = row["Receptor"]
            base = row["BaseLigand"]

            if base not in seen_bases[rec] and len(seen_bases[rec]) < top_n:
                grouped[rec].append(row)
                seen_bases[rec].add(base)

    else:
        raise SystemExit("❌ Unknown mode: choose 'per_ligand' or 'per_receptor'.")

    return grouped, style

# -------------------- locating docked out.pdbqt --------------------

def resolve_outfile_provenance(row: dict, cwd: Path) -> Path | None:
    out = (row.get("OutFile") or "").strip()
    if out:
        p = Path(out)
        if p.is_file():
            return p
        if p.is_absolute() and not p.exists():
            tail = Path(*p.parts[-4:]) if len(p.parts) >= 4 else p
            cand = cwd / tail
            if cand.is_file():
                return cand
    rr = (row.get("ResultsRoot") or "").strip()
    rd = (row.get("ReceptorDir") or "").strip()
    ld = (row.get("LigandDir") or "").strip()
    if rr and rd and ld:
        guess = cwd / rr / rd / ld / "out.pdbqt"
        if guess.is_file():
            return guess
    return None

def find_docked_pose_simple(cwd: Path, receptor: str, ligand_field: str) -> Path | None:
    rec_vars = norm_variants(receptor)
    lig_vars = norm_variants(ligand_field)

    for root in cwd.glob("Docking_Results*"):
        for rv in rec_vars:
            rec_dir = root / rv
            if rec_dir.is_dir():
                for lv in lig_vars:
                    cand = rec_dir / lv / "out.pdbqt"
                    if cand.is_file():
                        return cand

    for root in cwd.glob("Docking_Results*"):
        for rv in rec_vars:
            rec_root = root / rv
            if rec_root.is_dir():
                for p in rec_root.rglob("out.pdbqt"):
                    parent_str = str(p.parent).lower()
                    if any(lv in parent_str for lv in lig_vars):
                        return p

    for p in cwd.rglob("*.pdbqt"):
        sp = str(p).lower()
        if any(lv in sp for lv in lig_vars):
            return p
    return None

def resolve_outfile_unified(row: dict, cwd: Path, style: str) -> Path | None:
    if style == "provenance":
        p = resolve_outfile_provenance(row, cwd)
        if p and p.exists():
            return p
    receptor = row["Receptor"]
    ligand_field = row["Ligand"]
    return find_docked_pose_simple(cwd, receptor, ligand_field)

# -------------------- reference SDF/MOL2 lookup --------------------

def infer_bucket_from_row(row: dict) -> str | None:
    for key in ("ResultsRoot","LigandDir","OutFile"):
        s = (row.get(key) or "")
        m = BUCKET_PAT.search(s)
        if m:
            return m.group(1).upper()
    return None

def get_mol_from_sdf_dir(sdf_dir: Path, base_id: str):
    want = norm_base(strip_pose_suffix(base_id))
    sdf_files = list(sdf_dir.rglob("*.sdf"))
    if not sdf_files:
        return None, None

    for sdf in sdf_files:
        try:
            supp = Chem.SDMolSupplier(str(sdf), sanitize=False)
        except Exception:
            continue
        for m in supp:
            if m is None:
                continue
            title = m.GetProp("_Name") if m.HasProp("_Name") else ""
            if norm_base(title) == want:
                return m, sdf

    for sdf in sdf_files:
        if norm_base(sdf.stem) == want:
            try:
                supp = Chem.SDMolSupplier(str(sdf), sanitize=False)
                m = next((x for x in supp if x is not None), None)
                if m:
                    return m, sdf
            except Exception:
                pass
    return None, None

def get_ref_from_simple_folder(ligands_root: Path, base_id_or_lig_field: str):
    want = norm_base(strip_pose_suffix(base_id_or_lig_field))
    cand = None
    for ext in (".sdf", ".mol2"):
        for p in ligands_root.rglob(f"*{ext}"):
            if norm_base(p.stem) == want:
                cand = p
                break
        if cand:
            break
    if not cand:
        return None, None
    if cand.suffix.lower() == ".sdf":
        sup = Chem.SDMolSupplier(str(cand), sanitize=False)
        ref = next((m for m in sup if m is not None), None)
    else:
        ref = Chem.MolFromMol2File(str(cand), sanitize=False)
    return ref, cand

# -------------------- chemistry: PDBQT -> coords, rigid MCS fit --------------------

def pdbqt_to_mol(pdbqt_file, tmp_mol_file="/tmp/tmp_ligand_coords.mol"):
    subprocess.run(
        ["obabel", pdbqt_file, "-O", tmp_mol_file],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    return Chem.MolFromMolFile(tmp_mol_file, sanitize=False, removeHs=False)

def _noH_with_map(m: Chem.Mol):
    idx_map = []
    for i, a in enumerate(m.GetAtoms()):
        if a.GetAtomicNum() > 1:
            idx_map.append(i)
    noH = Chem.RemoveHs(Chem.Mol(m), sanitize=False)
    return noH, idx_map

def rigid_fit_by_mcs(coords_full: Chem.Mol, ref_full: Chem.Mol) -> Chem.Mol | None:
    """
    Rigidly align the ENTIRE reference (ref_full) onto the docked coordinates (coords_full)
    using a heavy-atom MCS atom map. This avoids per-atom coordinate overwrites and
    prevents 'flying atoms'.
    """
    if coords_full.GetNumConformers() == 0 or ref_full.GetNumConformers() == 0:
        print("⚠️ Missing conformers; cannot align.")
        return None

    # Work on heavy atoms for matching; keep maps to full indices
    coords_noH, coords_map = _noH_with_map(coords_full)
    ref_noH,   ref_map     = _noH_with_map(ref_full)

    mcs = rdFMCS.FindMCS([ref_noH, coords_noH])
    if mcs.canceled or not mcs.smartsString:
        print("⚠️ MCS alignment canceled/failed")
        return None
    patt = Chem.MolFromSmarts(mcs.smartsString)
    ref_match    = list(ref_noH.GetSubstructMatch(patt))
    coords_match = list(coords_noH.GetSubstructMatch(patt))
    if len(ref_match) != len(coords_match) or len(ref_match) == 0:
        print("⚠️ MCS mismatch in atom counts")
        return None

    # Build atomMap in FULL indices (pairs of (refIdx, coordsIdx))
    atomMap = []
    for i_ref_noH, i_coords_noH in zip(ref_match, coords_match):
        ref_idx_full    = ref_map[i_ref_noH]
        coords_idx_full = coords_map[i_coords_noH]
        atomMap.append((ref_idx_full, coords_idx_full))

    # Make a copy; AlignMol moves the probe (first arg) to match the reference (second arg)
    ref_aln = Chem.Mol(ref_full)
    try:
        rdMolAlign.AlignMol(ref_aln, coords_full, atomMap=atomMap)
    except Exception as e:
        print(f"⚠️ AlignMol failed: {e}")
        return None

    # (Optional) remove hydrogens to match your SDF style
    ref_aln = Chem.RemoveHs(Chem.Mol(ref_aln), sanitize=False)
    try:
        Chem.SanitizeMol(
            ref_aln,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES
                        | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
                        | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
    except Exception:
        pass

    return ref_aln

# -------------------- filesystem helpers --------------------

def find_receptor_file(receptors_root: Path, rec_name: str) -> Path | None:
    stems = {rec_name, rec_name.replace(".pdbqt","").replace(".converted","")}
    exts  = [".pdb", ".pdbqt", ".mol2", ".converted.pdbqt", ""]
    for root, _, files in os.walk(receptors_root):
        for f in files:
            p = Path(root) / f
            s = p.stem
            if s in stems or p.name.startswith(rec_name):
                if p.suffix in exts or any(str(p).endswith(ext) for ext in exts):
                    return p
    return None

def strip_conect_lines(pdb_path: Path):
    try:
        txt = pdb_path.read_text()
        lines = [ln for ln in txt.splitlines() if not ln.startswith("CONECT")]
        pdb_path.write_text("\n".join(lines) + "\n")
    except Exception as e:
        print(f"⚠️ Could not strip CONECT from {pdb_path}: {e}")

# -------------------- main --------------------

def main():
    cwd = Path(".").resolve()
    out_dir_pdb   = cwd / "11_AA_PDBOUTPUT"
    out_dir_sdf   = cwd / "11_AA_Ligands"
    out_dir_pdbqt = cwd / "11_AA_PDBQT"
    out_dir_pdb.mkdir(exist_ok=True)
    out_dir_sdf.mkdir(exist_ok=True)
    out_dir_pdbqt.mkdir(exist_ok=True)

    cmd.set('pdb_conect_all', 0)
    cmd.set('pdb_conect_nodup', 1)

    # 1) Select CSV
    csv_files = [f for f in os.listdir(cwd) if f.endswith(".csv")]
    if not csv_files:
        raise SystemExit("❌ No CSV files found in current directory.")
    csv_file = list_and_select("Available CSV files:", csv_files)
    csv_path = cwd / csv_file

    # 2) Mode selection
    print("\n⚙️ Selection mode:")
    print("1. Top N poses per ligand (per receptor–ligand pair)  ← multi-state per base ligand")
    print("2. Top N ligands total (per receptor, best pose only)")
    mode_choice = input("🔢 Choose [1/2]: ").strip()
    mode = "per_ligand" if mode_choice == "1" else ("per_receptor" if mode_choice == "2" else None)
    if not mode:
        raise SystemExit("❌ Invalid choice.")

    # 3) N
    top_n = int(input("🔝 How many top entries? (e.g., 5): "))

    # 4) Pick Receptors folder
    rec_dirs = [d for d in os.listdir(cwd) if os.path.isdir(d) and "Receptor" in d]
    if not rec_dirs:
        raise SystemExit("❌ No 'Receptor*' folder found in current dir.")
    receptors_root = cwd / list_and_select("Receptor folders:", rec_dirs)

    # 5) Load + group rows & detect style
    top_poses, style = load_csv_grouped(csv_path, top_n, mode=mode)
    print(f"🔎 Detected CSV style: {style}")

   # 6) If simple style, ask for Ligands folder once (reference SDF/MOL2 by base name)
    ligands_root = None
    if style == "simple":
        lig_dirs = [d for d in os.listdir(cwd) if os.path.isdir(d) and "Ligand" in d]
        if not lig_dirs:
            print("ℹ️ Simple mode: no 'Ligand*' folder found; will proceed without ref SDF (coords only).")
        else:
            ligands_root = cwd / list_and_select("Ligand folders (SDF/MOL2, base-name match):", lig_dirs)

    # 7) Build PyMOL session + exports
    print("\n🔬 Building PyMOL session...")
    cmd.reinitialize()
    try:
        cmd.feedback("disable", "all", "warnings")  # quiet name-cleanup spam
    except Exception:
        pass

    per_receptor_rank = defaultdict(int)
    audit_rows = []
    ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

    loaded_objects = set()

    for receptor, rows in top_poses.items():
        rec_name = receptor.replace(".pdbqt", "").replace(".converted", "")
        rec_obj  = sanitize_pymol_name(rec_name)

        rec_file = find_receptor_file(receptors_root, rec_name)
        if rec_file is None:
            print(f"❌ Receptor file not found for {rec_name}. Skipping receptor.")
            continue

        try:
            cmd.load(str(rec_file), rec_obj)
        except Exception as e:
            print(f"⚠️ Failed to load receptor {rec_file}: {e}")
            continue

        # regroup by base ligand to build 1 object per base with multiple states
        by_base = defaultdict(list)
        for row in rows:
            base = row["BaseLigand"]
            by_base[base].append(row)

        for base_ligand, base_rows in by_base.items():
            base_tag = sanitize_pymol_name(f"{rec_name}_{base_ligand}")
            lig_tag  = sanitize_pymol_name(base_ligand)
            state_index = 0

            for row in base_rows:
                ligand_field = row["Ligand"]             # with __pose##
                pose_num     = (row.get("Pose") or "").strip() or extract_pose_num(ligand_field) or ""

                dock_path = resolve_outfile_unified(row, cwd, style)
                if not dock_path or not dock_path.exists():
                    print(f"❌ Missing docked pose for {rec_name}/{ligand_field} (base {base_ligand})")
                    continue

                try:
                    coords = pdbqt_to_mol(str(dock_path))
                except subprocess.CalledProcessError:
                    print(f"⚠️ obabel failed on {dock_path}")
                    continue
                if coords is None:
                    print(f"⚠️ Failed to parse coords from {dock_path}")
                    continue

                # Locate reference SDF/MOL2
                ref = ref_src = None
                if style == "provenance":
                    bucket = infer_bucket_from_row(row)
                    if bucket:
                        sdf_root = None
                        for d in cwd.iterdir():
                            if d.is_dir() and d.name.startswith(f"Ligands_{bucket}") and "_PDBQT" not in d.name:
                                sdf_root = d; break
                        if not sdf_root:
                            for rootd, dirs, _ in os.walk(cwd):
                                for name in dirs:
                                    if name.startswith(f"Ligands_{bucket}") and "_PDBQT" not in name:
                                        sdf_root = Path(rootd)/name; break
                                if sdf_root: break
                        if sdf_root:
                            ref, ref_src = get_mol_from_sdf_dir(sdf_root, base_ligand or ligand_field)
                else:
                    if ligands_root:
                        ref, ref_src = get_ref_from_simple_folder(ligands_root, base_ligand or ligand_field)

                # Build rigidly fitted molecule (or coords-only) and append as a NEW STATE
                state_index += 1
                if ref is None:
                    tmp_coords = Path("/tmp") / f"{base_tag}_coords_only_state{state_index}.mol"
                    Chem.MolToMolFile(coords, str(tmp_coords))
                    cmd.load(str(tmp_coords), base_tag, state=state_index, discrete=0)
                else:
                    ref_orig = Chem.Mol(ref)
                    fitted = rigid_fit_by_mcs(coords, Chem.Mol(ref))
                    if fitted is None:
                        # fallback: coords-only to avoid “flying atoms”
                        tmp_coords = Path("/tmp") / f"{base_tag}_coords_only_state{state_index}.mol"
                        Chem.MolToMolFile(coords, str(tmp_coords))
                        cmd.load(str(tmp_coords), base_tag, state=state_index, discrete=0)
                    else:
                        tmp_path = Path("/tmp") / f"{base_tag}_state{state_index}.mol"
                        with open(tmp_path, "w") as fp:
                            fp.write(Chem.MolToMolBlock(fitted, kekulize=False))
                        cmd.load(str(tmp_path), base_tag, state=state_index, discrete=0)

                        # export SDFs for reference and this fitted state
                        ref_sdf  = out_dir_sdf / f"{base_ligand}.sdf"
                        dock_sdf = out_dir_sdf / f"{base_ligand}_docked_state{state_index}.sdf"
                        try:
                            w = SDWriter(str(ref_sdf));  w.SetKekulize(False);  w.write(ref_orig);  w.close()
                        except Exception as e:
                            print(f"⚠️ Failed to write reference SDF {ref_sdf}: {e}")
                        try:
                            w = SDWriter(str(dock_sdf)); w.SetKekulize(False); w.write(fitted); w.close()
                        except Exception as e:
                            print(f"⚠️ Failed to write docked SDF {dock_sdf}: {e}")

                # Export complex PDB (receptor + current state)
                # per-receptor running rank (Top01, Top02, … across loaded states)
                per_receptor_rank[rec_obj] += 1
                rank = per_receptor_rank[rec_obj]
                rank2 = f"{rank:02d}"
                pdb_path   = out_dir_pdb   / f"{rec_name}_Top{rank2}_{lig_tag}_state{state_index}.pdb"
                pdbqt_name = out_dir_pdbqt / f"{rec_name}_Top{rank2}_{lig_tag}_state{state_index}.pdbqt"

                sel = f"({rec_obj}) or ({base_tag} and state {state_index})"
                try:
                    cmd.save(str(pdb_path), selection=sel)
                    strip_conect_lines(pdb_path)
                except Exception as e:
                    print(f"⚠️ Failed to save PDB {pdb_path}: {e}")

                # copy raw docked out.pdbqt with a unique name
                try:
                    shutil.copy2(dock_path, pdbqt_name)
                except Exception as e:
                    print(f"⚠️ Failed to copy PDBQT to {pdbqt_name}: {e}")

                # Group once objects exist; ignore error if already grouped
                try:
                    cmd.group(rec_obj, base_tag)
                    cmd.group(f"{rec_obj}/{lig_tag}", base_tag)
                except Exception:
                    pass

                # Audit
                audit_rows.append({
                    "timestamp": ts,
                    "csv_style": style,
                    "mode": mode_choice,
                    "receptor": rec_name,
                    "rank_for_receptor": rank,
                    "base_ligand": base_ligand,
                    "csv_pose": pose_num,
                    "dock_pdbqt": str(dock_path),
                    "ref_source": str(ref_src) if ref_src else "",
                    "saved_complex_pdb": str(pdb_path),
                    "saved_dock_pdbqt": str(pdbqt_name),
                    "session_object": base_tag,
                    "state_index": state_index,
                })

    # Save session + audit
    ts2 = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    session_name = out_dir_pdb / f"Top{top_n}_{'per_ligand' if mode=='per_ligand' else 'per_receptor'}_Docking_Results_{ts2}.pse"
    cmd.save(str(session_name))

    audit_csv = out_dir_pdb / f"top_hits_audit_{ts2}.csv"
    if audit_rows:
        with open(audit_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()))
            writer.writeheader()
            writer.writerows(audit_rows)

    print(f"\n🎉 Done!")
    print(f"🧪 PyMOL session: {session_name}")
    print(f"📦 PDB complexes: {out_dir_pdb}")
    print(f"📦 Ligand SDFs:   {out_dir_sdf}")
    print(f"📦 Raw PDBQT:     {out_dir_pdbqt}")
    if audit_rows:
        print(f"🧾 Audit CSV:     {audit_csv}")

if __name__ == "__main__":
    main()



