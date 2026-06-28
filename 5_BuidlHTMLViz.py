#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import re
import shutil
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

UPSTREAM_REPO_URL = "https://github.com/muntisa/py-VinaScope-Docking-Viewer"
UPSTREAM_VIEWER_URL = "https://muntisa.github.io/VinaDock-Viz/VinaDock_Viz.html"
UPSTREAM_AUTHOR = "Cristian R. Munteanu, PhD"
UPSTREAM_AUTHOR_AFFILIATION = "Professor of Computer Science, University of A Coruna, RNASA-IMEDIR"

UPSTREAM_LICENSE_TEXT = """MIT License

Copyright (c) 2026 Cristian R Munteanu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

POSE_RE = re.compile(r"^REMARK\s+VINA\s+RESULT:\s+(-?\d+(?:\.\d+)?)", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a portable HTML docking-visualization project from a sorted Vina score CSV."
    )
    parser.add_argument("--csv", help="Sorted docking CSV produced by 4_ParseScores.py or 4C_ConcatenateScores.py")
    parser.add_argument("--outdir", help="Output directory for the generated visualization project")
    parser.add_argument(
        "--receptor-roots",
        nargs="+",
        default=["Receptors", "Receptors_PDBQT", "."],
        help="Folders searched for receptor PDBQT files",
    )
    parser.add_argument("--top-ligands", type=int, default=25, help="Top ligand hits per receptor to include")
    parser.add_argument("--top-poses", type=int, default=3, help="Top poses retained inside each viewer")
    parser.add_argument("--project-name", default="Docking_HTML_Viz_Project", help="Output project folder name")
    parser.add_argument("--page-title", default="Docking Visualization Project", help="Display title embedded in the HTML pages")
    return parser.parse_args()


def choose_from_list(prompt: str, items: List[Path]) -> Path:
    print(f"\n{prompt}")
    for index, item in enumerate(items, start=1):
        print(f" [{index}] {item.name}")
    while True:
        raw = input("Enter index: ").strip()
        try:
            selected = int(raw)
        except ValueError:
            print("Please enter a number from the list.")
            continue
        if 1 <= selected <= len(items):
            return items[selected - 1]
        print("Selection out of range. Try again.")


def prompt_text(prompt: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{prompt}{suffix}: ").strip()
    return value or default


def prompt_int(prompt: str, default: int) -> int:
    while True:
        raw = input(f"{prompt} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = int(raw)
        except ValueError:
            print("Please enter an integer.")
            continue
        if value < 1:
            print("Please enter a value greater than 0.")
            continue
        return value


def prompt_csv_path() -> Path:
    candidates = sorted(Path.cwd().glob("*vina_docking_scores_sorted.csv"))
    if not candidates:
        candidates = sorted(Path.cwd().glob("*.csv"))
    if candidates:
        return choose_from_list("Select a docking score CSV:", candidates)

    while True:
        raw = input("Path to docking score CSV: ").strip()
        if not raw:
            print("A CSV path is required.")
            continue
        path = Path(raw).expanduser()
        if path.exists() and path.is_file():
            return path
        print(f"File not found: {path}")


def default_receptor_roots_for(csv_path: Path) -> List[str]:
    base = csv_path.resolve().parent
    defaults = ["Receptors", "Receptors_PDBQT", "."]
    existing: List[str] = []
    for candidate in defaults:
        path = base / candidate if candidate != "." else base
        if path.exists():
            existing.append(candidate)
    return existing or defaults


def resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.csv:
        return args

    print("\n5_BuidlHTMLViz interactive mode")
    print("No --csv flag was provided, so we'll gather the inputs step by step.")

    csv_path = prompt_csv_path()
    args.csv = str(csv_path)

    suggested_outdir = str(csv_path.resolve().parent)
    args.outdir = prompt_text("Output directory", suggested_outdir)

    roots_default = " ".join(default_receptor_roots_for(csv_path))
    receptor_roots_raw = prompt_text("Receptor search roots (space-separated)", roots_default)
    args.receptor_roots = receptor_roots_raw.split()

    args.top_ligands = prompt_int("Top ligands per receptor", args.top_ligands)
    args.top_poses = prompt_int("Top poses per ligand", args.top_poses)
    args.project_name = prompt_text("Project folder name", args.project_name)
    args.page_title = prompt_text("Viewer page title", args.page_title)
    return args


def safe_slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (value or "").strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("._")
    return cleaned or "item"


def safe_float(value: Any, default: float = 1e9) -> float:
    try:
        return float(value)
    except Exception:
        return default


def load_rows(csv_path: Path) -> List[Dict[str, str]]:
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [{str(key): str(value or "") for key, value in row.items()} for row in reader]
    rows.sort(
        key=lambda row: (
            row.get("Receptor", ""),
            safe_float(row.get("Binding_Affinity", "")),
            safe_float(row.get("Pose", ""), default=999999),
        )
    )
    return rows


def group_top_hits(rows: Iterable[Dict[str, str]], top_ligands: int) -> List[Dict[str, str]]:
    per_receptor: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
    for row in rows:
        receptor = row.get("Receptor", "").strip()
        ligand_key = row.get("LigandBase") or row.get("Ligand") or row.get("LigandVariant") or "ligand"
        variant_key = row.get("LigandVariant") or row.get("Ligand") or "variant"
        bucket = per_receptor[receptor]
        existing = bucket.get(variant_key)
        current_score = safe_float(row.get("Binding_Affinity", ""))
        if existing is None or current_score < safe_float(existing.get("Binding_Affinity", "")):
            row_copy = dict(row)
            row_copy["_LigandBaseKey"] = ligand_key
            row_copy["_LigandVariantKey"] = variant_key
            bucket[variant_key] = row_copy

    selected: List[Dict[str, str]] = []
    for receptor in sorted(per_receptor):
        variants = list(per_receptor[receptor].values())
        base_best: Dict[str, float] = {}
        for row in variants:
            base_key = row.get("_LigandBaseKey", "")
            score = safe_float(row.get("Binding_Affinity", ""))
            base_best[base_key] = min(score, base_best.get(base_key, score))

        ranked_bases = sorted(base_best.items(), key=lambda item: (item[1], item[0]))
        allowed_bases = {base for base, _score in ranked_bases[: max(0, top_ligands)]}

        ranked_variants = sorted(
            [row for row in variants if row.get("_LigandBaseKey", "") in allowed_bases],
            key=lambda row: (
                base_best.get(row.get("_LigandBaseKey", ""), 1e9),
                safe_float(row.get("Binding_Affinity", "")),
                row.get("_LigandBaseKey", ""),
                row.get("_LigandVariantKey", ""),
            ),
        )
        selected.extend(ranked_variants)
    return selected


def find_receptor_file(receptor_name: str, search_roots: Iterable[Path]) -> Optional[Path]:
    stem = Path(receptor_name).stem
    candidates = [
        receptor_name,
        f"{stem}.pdbqt",
        f"{stem}.pdb",
        f"{stem}.mol2",
    ]
    for root in search_roots:
        for candidate in candidates:
            path = (root / candidate).resolve()
            if path.exists():
                return path
    return None


def split_pose_blocks(pdbqt_text: str) -> List[str]:
    blocks: List[str] = []
    current: List[str] = []
    in_model = False
    for line in pdbqt_text.splitlines(keepends=True):
        if line.startswith("MODEL"):
            if current:
                blocks.append("".join(current))
            current = [line]
            in_model = True
            continue
        if line.startswith("ENDMDL"):
            current.append(line)
            blocks.append("".join(current))
            current = []
            in_model = False
            continue
        if in_model:
            current.append(line)
    if current:
        blocks.append("".join(current))
    return blocks


def trim_pose_file_text(pdbqt_text: str, top_poses: int) -> str:
    blocks = split_pose_blocks(pdbqt_text)
    if not blocks:
        return pdbqt_text
    keep = blocks[: max(1, top_poses)]
    return "".join(keep)


def pose_affinities_from_text(pdbqt_text: str) -> List[float]:
    return [float(match.group(1)) for match in POSE_RE.finditer(pdbqt_text)]


def b64encode_text(text: str) -> str:
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


def load_vinascope_generator_template() -> str:
    template_path = Path(__file__).resolve().parent / "viewer_templates" / "generate_vina_docking_viewer_upstream.py"
    script_text = template_path.read_text(encoding="utf-8")
    match = re.search(r'HTML_TEMPLATE\s*=\s*r"""(.*)"""', script_text, re.DOTALL)
    if not match:
        raise RuntimeError("Could not extract HTML_TEMPLATE from bundled VinaScope generator.")
    return match.group(1)


def build_viewer_html(
    page_title: str,
    receptor_label: str,
    ligand_label: str,
    receptor_text: str,
    ligand_text: str,
    best_affinity: str,
    source_outfile: str,
) -> str:
    title = f"{page_title} | {ligand_label}"
    receptor_b64 = b64encode_text(receptor_text)
    ligand_b64 = b64encode_text(ligand_text)
    template = load_vinascope_generator_template()
    template = "<!-- Adapted from py-VinaScope-Docking-Viewer | " + html.escape(source_outfile) + " -->\n" + template
    template = template.replace("%%RECEPTOR_B64%%", receptor_b64)
    template = template.replace("%%POSE_B64%%", ligand_b64)
    template = template.replace(
        "<title>VinaScope · Molecular Docking Viewer</title>",
        f"<title>{html.escape(title)}</title>",
        1,
    )
    template = template.replace(
        'return m ? m[1] : "Ligand";',
        f"return m ? m[1] : {json.dumps(ligand_label)};",
        1,
    )
    template = template.replace(
        'document.getElementById("ligand-label").querySelector("span").textContent = "Ligand: "+ligandName;',
        'document.getElementById("ligand-label").querySelector("span").textContent = "Ligand: "+ligandName;\n'
        f'  document.title = {json.dumps(title)};\n',
        1,
    )
    template = template.replace(
        'toast("Interactions CSV saved automatically");',
        f'toast("Loaded {ligand_label} ({best_affinity} kcal/mol)");',
        1,
    )

    binding_pocket_controls = """        <label class="toggle-row">
          <input type="checkbox" id="surface-toggle">
          <span>VDW surface overlay</span>
        </label>
      </div>

      <div class="divider"></div>

      <!-- BINDING POCKET STYLE -->
      <div class="section">
        <div class="sec-label">Binding Pocket</div>
        <div class="ctrl-group">
          <div class="ctrl-label">Style</div>
          <div class="btn-row" id="pocket-style-btns">
            <button class="btn active" data-v="stick">Sticks</button>
            <button class="btn" data-v="sphere">Spheres</button>
            <button class="btn" data-v="line">Lines</button>
          </div>
        </div>
        <label class="toggle-row">
          <input type="checkbox" id="pocket-toggle" checked>
          <span>Show atoms within 5 Å of active ligand</span>
        </label>
      </div>

      <div class="divider"></div>"""
    template = template.replace(
        """        <label class="toggle-row">
          <input type="checkbox" id="surface-toggle">
          <span>VDW surface overlay</span>
        </label>
      </div>

      <div class="divider"></div>""",
        binding_pocket_controls,
        1,
    )

    template = template.replace(
        "const COORD_RX = /^\\s*(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)/;",
        """const COORD_RX = /^\\s*(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)\\s+(-?\\d+\\.\\d+)/;
const POCKET_DISTANCE_CUTOFF = 5.0;
const POCKET_HIGHLIGHT = "#f97316";""",
        1,
    )
    template = template.replace(
        """    const name    = line.substring(12,16).trim();
    const resname = line.substring(17,20).trim();
    const chain   = line.substring(21,22).trim() || "A";
    const resi    = parseInt(line.substring(22,26)) || 0;
""",
        """    const serial  = parseInt(line.substring(6,11).trim()) || 0;
    const name    = line.substring(12,16).trim();
    const resname = line.substring(17,20).trim();
    const chain   = line.substring(21,22).trim() || "A";
    const resi    = parseInt(line.substring(22,26)) || 0;
""",
        1,
    )
    template = template.replace(
        "    atoms.push({ name, resname, chain, resi, x, y, z, charge, adtype });",
        "    atoms.push({ serial, name, resname, chain, resi, x, y, z, charge, adtype });",
        1,
    )
    template = template.replace(
        """function groupByResidue(rows) {
  const map = new Map();
  for (const r of rows) {
    const key = `${r.receptor_resname} ${r.receptor_resi}:${r.receptor_chain}`;
    if (!map.has(key)) map.set(key, { key, types:new Set(), minDist:Infinity });
    const g = map.get(key);
    g.types.add(r.type);
    g.minDist = Math.min(g.minDist, r.distance);
  }
  return [...map.values()].sort((a,b)=>a.minDist-b.minDist);
}
""",
        """function groupByResidue(rows) {
  const map = new Map();
  for (const r of rows) {
    const key = `${r.receptor_resname} ${r.receptor_resi}:${r.receptor_chain}`;
    if (!map.has(key)) map.set(key, { key, types:new Set(), minDist:Infinity });
    const g = map.get(key);
    g.types.add(r.type);
    g.minDist = Math.min(g.minDist, r.distance);
  }
  return [...map.values()].sort((a,b)=>a.minDist-b.minDist);
}

function buildPocketSelection(ligAtoms, recAtoms, cutoff) {
  const serials = new Set();
  for (const la of ligAtoms) {
    for (const ra of recAtoms) {
      if (dist3(la, ra) <= cutoff) serials.add(ra.serial);
    }
  }
  return [...serials];
}
""",
        1,
    )
    template = template.replace(
        """  const recModel = viewer.addModel(receptorPDB,"pdb");

  let recStyle="cartoon", recColor="spectrum", recSurfObj=null, surfObj=null, surfOn=false;
""",
        """  const recModel = viewer.addModel(receptorPDB,"pdb");

  let recStyle="cartoon", recColor="spectrum", recSurfObj=null, surfObj=null, surfOn=false;
  let pocketOn=true, pocketStyle="stick";
  let pocketSelection = [];
""",
        1,
    )
    template = template.replace(
        """  function applyRecStyle() {
    viewer.setStyle({model:recModel.getID()},{});
    const cs = recColor==="chain"?"chain": recColor==="ss"?"ssJmol": undefined;
    if (recSurfObj) { viewer.removeSurface(recSurfObj); recSurfObj = null; }
    if (recStyle==="cartoon")
      viewer.setStyle({model:recModel.getID()},{cartoon:{color:recColor==="spectrum"?"spectrum":undefined,colorscheme:cs}});
    else if (recStyle==="surface") {
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:0.3}});
      recSurfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.88,colorscheme:$3Dmol.elementColors.rasmol},{model:recModel.getID()});
    }
    else if (recStyle==="line")
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:1.5}});
    else
      viewer.setStyle({model:recModel.getID()},{stick:{colorscheme:cs??"element",radius:0.12}});
    if (surfObj) { viewer.removeSurface(surfObj); surfObj=null; }
    if (surfOn)
      surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.35,color:"#0a2040"},{model:recModel.getID()});
  }
""",
        """  function applyRecStyle() {
    viewer.setStyle({model:recModel.getID()},{});
    const cs = recColor==="chain"?"chain": recColor==="ss"?"ssJmol": undefined;
    if (recSurfObj) { viewer.removeSurface(recSurfObj); recSurfObj = null; }
    if (recStyle==="cartoon")
      viewer.setStyle({model:recModel.getID()},{cartoon:{color:recColor==="spectrum"?"spectrum":undefined,colorscheme:cs}});
    else if (recStyle==="surface") {
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:0.3}});
      recSurfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.88,colorscheme:$3Dmol.elementColors.rasmol},{model:recModel.getID()});
    }
    else if (recStyle==="line")
      viewer.setStyle({model:recModel.getID()},{line:{colorscheme:cs??"element",linewidth:1.5}});
    else
      viewer.setStyle({model:recModel.getID()},{stick:{colorscheme:cs??"element",radius:0.12}});
    if (surfObj) { viewer.removeSurface(surfObj); surfObj=null; }
    if (surfOn)
      surfObj = viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.35,color:"#0a2040"},{model:recModel.getID()});

    if (pocketOn && pocketSelection.length) {
      const pocketSel = {model:recModel.getID(), serial:pocketSelection};
      if (pocketStyle==="sphere")
        viewer.setStyle(pocketSel,{sphere:{color:POCKET_HIGHLIGHT,radius:0.45,opacity:0.95}});
      else if (pocketStyle==="line")
        viewer.setStyle(pocketSel,{line:{color:POCKET_HIGHLIGHT,linewidth:3}});
      else
        viewer.setStyle(pocketSel,{stick:{color:POCKET_HIGHLIGHT,radius:0.2,opacity:0.98}});
    }
  }
""",
        1,
    )
    template = template.replace(
        """  const recAtoms      = parseAtoms(receptorPDB);
  const allInteractions = poses.map((_,i)=>
    computeInteractions(parseAtoms(poseToPDB(poses[i])), recAtoms)
  );
""",
        """  const recAtoms      = parseAtoms(receptorPDB);
  const ligandAtomsByPose = poses.map((_,i)=>parseAtoms(poseToPDB(poses[i])));
  const allInteractions = poses.map((_,i)=>
    computeInteractions(ligandAtomsByPose[i], recAtoms)
  );
  const pocketSelections = poses.map((_,i)=>
    buildPocketSelection(ligandAtomsByPose[i], recAtoms, POCKET_DISTANCE_CUTOFF)
  );
""",
        1,
    )
    template = template.replace(
        """  function selectPose(i) {
    curPose=i;
    document.querySelectorAll(".pose-item").forEach((el,j)=>el.classList.toggle("active",j===i));
    const hex=poseHex(i);
    document.getElementById("hud-score").textContent = poses[i].score.toFixed(1);
    document.getElementById("hud-score").style.color = hex;
    document.getElementById("hud-pose").textContent  = `Pose ${i+1} of ${poses.length}`;
    applyLigStyles();
    renderInteractions(allInteractions[i], i);
    viewer.render();
  }
""",
        """  function selectPose(i) {
    curPose=i;
    pocketSelection = pocketSelections[i] || [];
    document.querySelectorAll(".pose-item").forEach((el,j)=>el.classList.toggle("active",j===i));
    const hex=poseHex(i);
    document.getElementById("hud-score").textContent = poses[i].score.toFixed(1);
    document.getElementById("hud-score").style.color = hex;
    document.getElementById("hud-pose").textContent  = `Pose ${i+1} of ${poses.length}`;
    applyRecStyle();
    applyLigStyles();
    renderInteractions(allInteractions[i], i);
    viewer.render();
  }
""",
        1,
    )
    template = template.replace(
        """  wireGroup("#rec-style-btns .btn","v", v=>{ recStyle=v; applyRecStyle(); viewer.render(); });
  wireGroup("#rec-color-btns .btn","v", v=>{ recColor=v; applyRecStyle(); viewer.render(); });
  wireGroup("#lig-style-btns .btn","v", v=>{ ligStyle=v; applyLigStyles(); viewer.render(); });
  wireGroup("#lig-color-btns .btn","v", v=>{ ligColor=v; applyLigStyles(); viewer.render(); });
""",
        """  wireGroup("#rec-style-btns .btn","v", v=>{ recStyle=v; applyRecStyle(); viewer.render(); });
  wireGroup("#rec-color-btns .btn","v", v=>{ recColor=v; applyRecStyle(); viewer.render(); });
  wireGroup("#pocket-style-btns .btn","v", v=>{ pocketStyle=v; applyRecStyle(); viewer.render(); });
  wireGroup("#lig-style-btns .btn","v", v=>{ ligStyle=v; applyLigStyles(); viewer.render(); });
  wireGroup("#lig-color-btns .btn","v", v=>{ ligColor=v; applyLigStyles(); viewer.render(); });
""",
        1,
    )
    template = template.replace(
        """  document.getElementById("surface-toggle").addEventListener("change",e=>{
    surfOn=e.target.checked; applyRecStyle(); viewer.render();
  });
  document.getElementById("all-poses-toggle").addEventListener("change",e=>{
    showAll=e.target.checked; applyLigStyles(); viewer.render();
  });
""",
        """  document.getElementById("surface-toggle").addEventListener("change",e=>{
    surfOn=e.target.checked; applyRecStyle(); viewer.render();
  });
  document.getElementById("pocket-toggle").addEventListener("change",e=>{
    pocketOn=e.target.checked; applyRecStyle(); viewer.render();
  });
  document.getElementById("all-poses-toggle").addEventListener("change",e=>{
    showAll=e.target.checked; applyLigStyles(); viewer.render();
  });
""",
        1,
    )
    return template


def build_index_html(page_title: str, entries: List[Dict[str, Any]]) -> str:
    cards = []
    for entry in entries:
        cards.append(
            f"""
<article class="card">
  <div class="eyebrow">{html.escape(entry["receptor"])}</div>
  <h2>{html.escape(entry["ligand"])}</h2>
  <p>Best affinity: <strong>{html.escape(entry["best_affinity"])}</strong> kcal/mol</p>
  <a href="{html.escape(entry["viewer_file"])}">Open viewer</a>
</article>
"""
        )
    joined_cards = "\n".join(cards) or "<p>No viewer entries were generated.</p>"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{html.escape(page_title)}</title>
<style>
body {{
  margin: 0;
  font-family: "IBM Plex Sans", system-ui, sans-serif;
  background: linear-gradient(145deg, #f5efe4 0%, #edf5f2 100%);
  color: #182230;
}}
main {{
  max-width: 1100px;
  margin: 0 auto;
  padding: 42px 22px 60px;
}}
h1 {{ font-size: 42px; margin: 0 0 10px; }}
.lede {{ max-width: 760px; color: #4d5b69; line-height: 1.6; }}
.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
  gap: 16px;
  margin-top: 28px;
}}
.card {{
  background: rgba(255,255,255,0.84);
  border: 1px solid rgba(24,34,48,0.08);
  border-radius: 18px;
  padding: 18px;
  box-shadow: 0 18px 50px rgba(24, 34, 48, 0.08);
}}
.card h2 {{ margin: 0 0 8px; font-size: 22px; }}
.card p {{ margin: 0 0 14px; color: #4d5b69; }}
.card a {{
  display: inline-block;
  padding: 10px 14px;
  border-radius: 999px;
  background: #006d77;
  color: white;
  text-decoration: none;
}}
.eyebrow {{
  text-transform: uppercase;
  letter-spacing: 0.14em;
  color: #6a7684;
  font-size: 11px;
  margin-bottom: 10px;
}}
.credit {{
  margin-top: 34px;
  padding-top: 18px;
  border-top: 1px solid rgba(24,34,48,0.12);
  color: #5f6b7a;
  line-height: 1.6;
}}
.credit a {{ color: #006d77; }}
</style>
</head>
<body>
<main>
  <div class="eyebrow">Docking HTML Visualization</div>
  <h1>{html.escape(page_title)}</h1>
  <p class="lede">This portable project packages top docking hits into standalone HTML viewers that can be opened locally or embedded inside your docking server.</p>
  <section class="grid">
    {joined_cards}
  </section>
  <div class="credit">
    Adapted from <a href="{html.escape(UPSTREAM_REPO_URL)}" target="_blank" rel="noreferrer">py-VinaScope-Docking-Viewer</a> and
    <a href="{html.escape(UPSTREAM_VIEWER_URL)}" target="_blank" rel="noreferrer">VinaDock-Viz</a>.<br>
    Original source acknowledged to {html.escape(UPSTREAM_AUTHOR)}<br>
    {html.escape(UPSTREAM_AUTHOR_AFFILIATION)}
  </div>
</main>
</body>
</html>
"""


def write_zip(project_dir: Path) -> Path:
    zip_path = project_dir.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(project_dir.rglob("*")):
            archive.write(path, path.relative_to(project_dir.parent))
    return zip_path


def build_project(
    csv_path: Path,
    outdir: Optional[Path] = None,
    receptor_roots: Optional[Iterable[Path]] = None,
    top_ligands: int = 25,
    top_poses: int = 3,
    project_name: str = "Docking_HTML_Viz_Project",
    page_title: str = "Docking Visualization Project",
) -> Dict[str, Any]:
    csv_path = csv_path.resolve()
    base_dir = csv_path.parent
    output_root = (outdir or base_dir).resolve()
    project_dir = output_root / safe_slug(project_name)
    load_vinascope_generator_template()
    if project_dir.exists():
        shutil.rmtree(project_dir)
    viewers_dir = project_dir / "viewers"
    inputs_dir = project_dir / "inputs"
    project_dir.mkdir(parents=True, exist_ok=True)
    viewers_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)

    rows = load_rows(csv_path)
    selected_rows = group_top_hits(rows, top_ligands=top_ligands)
    search_roots: List[Path] = []
    for root in receptor_roots or []:
        root_path = Path(root)
        if root_path.is_absolute():
            search_roots.append(root_path.resolve())
        else:
            search_roots.append((Path.cwd() / root_path).resolve())
            search_roots.append((base_dir / root_path).resolve())
    if not search_roots:
        search_roots = [base_dir / "Receptors", base_dir / "Receptors_PDBQT", base_dir]

    manifest_entries: List[Dict[str, Any]] = []
    missing_entries: List[Dict[str, Any]] = []

    for row in selected_rows:
        receptor_name = row.get("Receptor", "").strip()
        ligand_name = row.get("LigandBase") or row.get("Ligand") or row.get("LigandVariant") or "ligand"
        ligand_variant = row.get("LigandVariant") or str(ligand_name)
        outfile = Path(row.get("OutFile", "")).resolve()
        receptor_file = find_receptor_file(receptor_name, search_roots)

        if not outfile.exists() or receptor_file is None:
            missing_entries.append(
                {
                    "receptor": receptor_name,
                    "ligand": ligand_name,
                    "outfile": str(outfile),
                    "receptor_file": str(receptor_file) if receptor_file else "",
                }
            )
            continue

        receptor_text = receptor_file.read_text(encoding="utf-8", errors="ignore")
        ligand_text = trim_pose_file_text(outfile.read_text(encoding="utf-8", errors="ignore"), top_poses=top_poses)
        affinities = pose_affinities_from_text(ligand_text)
        best_affinity = f"{min(affinities):.2f}" if affinities else row.get("Binding_Affinity", "")
        slug = safe_slug(f"{Path(receptor_name).stem}__{ligand_variant}")

        receptor_copy = inputs_dir / f"{slug}_receptor.pdbqt"
        ligand_copy = inputs_dir / f"{slug}_poses.pdbqt"
        viewer_file = viewers_dir / f"{slug}.html"

        receptor_copy.write_text(receptor_text, encoding="utf-8")
        ligand_copy.write_text(ligand_text, encoding="utf-8")
        try:
            viewer_html = build_viewer_html(
                page_title=page_title,
                receptor_label=Path(receptor_file).name,
                ligand_label=str(ligand_variant),
                receptor_text=receptor_text,
                ligand_text=ligand_text,
                best_affinity=best_affinity,
                source_outfile=str(outfile),
            )
        except Exception:
            if project_dir.exists():
                shutil.rmtree(project_dir, ignore_errors=True)
            raise
        viewer_file.write_text(viewer_html, encoding="utf-8")

        manifest_entries.append(
            {
                "receptor": receptor_name,
                "ligand": str(ligand_name),
                "ligand_variant": row.get("LigandVariant", ""),
                "best_affinity": best_affinity,
                "source_outfile": str(outfile),
                "receptor_file": str(receptor_file),
                "viewer_file": str(viewer_file.relative_to(project_dir)),
                "receptor_copy": str(receptor_copy.relative_to(project_dir)),
                "pose_copy": str(ligand_copy.relative_to(project_dir)),
            }
        )

    manifest = {
        "project_name": project_name,
        "page_title": page_title,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_csv": str(csv_path),
        "top_ligands": int(top_ligands),
        "top_poses": int(top_poses),
        "entry_count": len(manifest_entries),
        "entries": manifest_entries,
        "missing_entries": missing_entries,
        "attribution": {
            "repo": UPSTREAM_REPO_URL,
            "viewer": UPSTREAM_VIEWER_URL,
            "author": UPSTREAM_AUTHOR,
            "affiliation": UPSTREAM_AUTHOR_AFFILIATION,
        },
    }

    (project_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (project_dir / "UPSTREAM_LICENSE.txt").write_text(UPSTREAM_LICENSE_TEXT, encoding="utf-8")
    (project_dir / "index.html").write_text(build_index_html(page_title, manifest_entries), encoding="utf-8")
    zip_path = write_zip(project_dir)
    manifest["project_dir"] = str(project_dir)
    manifest["zip_path"] = str(zip_path)
    return manifest


def main() -> None:
    args = resolve_args(parse_args())
    manifest = build_project(
        csv_path=Path(args.csv),
        outdir=Path(args.outdir).resolve() if args.outdir else None,
        receptor_roots=[Path(path) for path in args.receptor_roots],
        top_ligands=args.top_ligands,
        top_poses=args.top_poses,
        project_name=args.project_name,
        page_title=args.page_title,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
