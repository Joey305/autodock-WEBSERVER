#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from ligand_naming import parse_ligand_variant


def _load_base_builder_module():
    module_path = Path(__file__).resolve().with_name("5_BuidlHTMLViz.py")
    spec = importlib.util.spec_from_file_location("html_viz_base_builder", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader is not None
    spec.loader.exec_module(module)
    return module


BASE = _load_base_builder_module()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compacted HTML docking-visualization project grouped by base ligand."
    )
    parser.add_argument("--csv", help="Sorted docking CSV produced by 4_ParseScores.py or 4C_ConcatenateScores.py")
    parser.add_argument("--outdir", help="Output directory for the generated visualization project")
    parser.add_argument(
        "--receptor-roots",
        nargs="+",
        default=["Receptors", "Receptors_PDBQT", "."],
        help="Folders searched for receptor PDBQT files",
    )
    parser.add_argument(
        "--top-ligands",
        type=int,
        default=25,
        help="Top base ligands per receptor to include",
    )
    parser.add_argument(
        "--top-poses",
        type=int,
        default=25,
        help="Top docking poses retained per base ligand across all conformers/states",
    )
    parser.add_argument(
        "--project-name",
        default="Docking_HTML_Viz_Project_Compacted",
        help="Output project folder name",
    )
    parser.add_argument(
        "--page-title",
        default="Compacted Docking Visualization Project",
        help="Display title embedded in the HTML pages",
    )
    return parser.parse_args()


def resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.csv:
        return args

    print("\n5_CompactedHTMLViz interactive mode")
    print("No --csv flag was provided, so we'll gather the inputs step by step.")

    csv_path = BASE.prompt_csv_path()
    args.csv = str(csv_path)

    suggested_outdir = str(csv_path.resolve().parent)
    args.outdir = BASE.prompt_text("Output directory", suggested_outdir)

    roots_default = " ".join(BASE.default_receptor_roots_for(csv_path))
    receptor_roots_raw = BASE.prompt_text("Receptor search roots (space-separated)", roots_default)
    args.receptor_roots = receptor_roots_raw.split()

    args.top_ligands = BASE.prompt_int("Top base ligands per receptor", args.top_ligands)
    args.top_poses = BASE.prompt_int("Top poses per base ligand", args.top_poses)
    args.project_name = BASE.prompt_text("Project folder name", args.project_name)
    args.page_title = BASE.prompt_text("Viewer page title", args.page_title)
    return args


def safe_float(value: Any, default: float = 1e9) -> float:
    try:
        return float(value)
    except Exception:
        return default


def unique_output_path(directory: Path, filename: str, used_names: set[str]) -> Path:
    candidate = filename
    stem = Path(filename).stem
    suffix = Path(filename).suffix
    counter = 2
    while candidate in used_names:
        candidate = f"{stem}_{counter}{suffix}"
        counter += 1
    used_names.add(candidate)
    return directory / candidate


def row_ligand_variant(row: Dict[str, str]) -> str:
    for key in ("LigandVariant", "Ligand", "LigandDir"):
        value = (row.get(key) or "").strip()
        if value:
            return Path(value).name
    out = (row.get("OutFile") or "").strip()
    if out:
        return Path(out).parent.name
    return ""


def row_ligand_base(row: Dict[str, str]) -> str:
    explicit = (row.get("LigandBase") or "").strip()
    if explicit:
        return explicit
    variant = row_ligand_variant(row)
    if variant:
        return str(parse_ligand_variant(variant).get("LigandBase", variant))
    ligand = (row.get("Ligand") or "").strip()
    if ligand:
        return str(parse_ligand_variant(ligand).get("LigandBase", ligand))
    return "ligand"


def normalize_row(row: Dict[str, str]) -> Dict[str, str]:
    cooked = dict(row)
    variant = row_ligand_variant(cooked)
    parsed = parse_ligand_variant(variant or cooked.get("Ligand", "") or "")
    cooked["LigandVariant"] = variant or str(parsed.get("LigandVariant", "")) or cooked.get("Ligand", "") or "ligand"
    cooked["LigandBase"] = row_ligand_base(cooked)
    cooked["StateTag"] = str(parsed.get("StateTag", "") or "")
    cooked["ProtomerTag"] = str(parsed.get("ProtomerTag", "") or "")
    cooked["TautomerTag"] = str(parsed.get("TautomerTag", "") or "")
    cooked["ConformerTag"] = str(parsed.get("ConformerTag", "") or "")
    return cooked


def load_rows(csv_path: Path) -> List[Dict[str, str]]:
    rows = [normalize_row(row) for row in BASE.load_rows(csv_path)]
    rows.sort(
        key=lambda row: (
            row.get("Receptor", ""),
            safe_float(row.get("Binding_Affinity", "")),
            safe_float(row.get("Pose", ""), default=999999),
            row.get("LigandBase", ""),
            row.get("LigandVariant", ""),
        )
    )
    return rows


def select_compacted_groups(
    rows: Sequence[Dict[str, str]],
    top_ligands: int,
    top_poses: int,
) -> List[Dict[str, Any]]:
    grouped_rows: Dict[str, Dict[str, List[Dict[str, str]]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        receptor = (row.get("Receptor") or "").strip()
        base = row_ligand_base(row)
        if not receptor or not base:
            continue
        grouped_rows[receptor][base].append(row)

    groups: List[Dict[str, Any]] = []
    for receptor in sorted(grouped_rows):
        base_groups = grouped_rows[receptor]
        ranked_bases = sorted(
            base_groups.items(),
            key=lambda item: (
                min(safe_float(row.get("Binding_Affinity", "")) for row in item[1]),
                item[0],
            ),
        )
        for base, base_rows in ranked_bases[: max(0, top_ligands)]:
            selected_rows = sorted(
                base_rows,
                key=lambda row: (
                    safe_float(row.get("Binding_Affinity", "")),
                    safe_float(row.get("Pose", ""), default=999999),
                    row.get("LigandVariant", ""),
                ),
            )[: max(1, top_poses)]
            groups.append(
                {
                    "receptor": receptor,
                    "ligand_base": base,
                    "selected_rows": selected_rows,
                }
            )
    return groups


def print_selection_summary(groups: Sequence[Dict[str, Any]], top_ligands: int, top_poses: int) -> None:
    per_receptor_counts: Dict[str, Counter[str]] = defaultdict(Counter)
    for group in groups:
        receptor = str(group["receptor"])
        ligand_base = str(group["ligand_base"])
        per_receptor_counts[receptor][ligand_base] = len(group["selected_rows"])

    print("\n✅ Compacted Top-N selection summary")
    print("   Mode: per_ligand_compacted")
    print(f"   Top base ligands per receptor: {top_ligands}")
    print(f"   Top poses per base ligand: {top_poses}")
    for receptor, counts in per_receptor_counts.items():
        total_rows = sum(counts.values())
        print(f"   - {receptor}: {total_rows} selected poses across {len(counts)} base ligand(s)")
        for ligand_base, count in sorted(counts.items()):
            print(f"       {ligand_base}: {count}")


def extract_selected_pose_text(pdbqt_text: str, pose_index: int) -> tuple[str, str]:
    lines = pdbqt_text.splitlines()
    models: List[Tuple[int, List[str]]] = []
    current_num: Optional[int] = None
    current_lines: List[str] = []

    for line in lines:
        if line.startswith("MODEL"):
            if current_lines:
                models.append((current_num or len(models) + 1, current_lines[:]))
            current_lines = [line]
            match = re.search(r"MODEL\s+(\d+)", line)
            current_num = int(match.group(1)) if match else len(models) + 1
            continue
        if current_lines:
            current_lines.append(line)
            if line.startswith("ENDMDL"):
                models.append((current_num or len(models) + 1, current_lines[:]))
                current_lines = []
                current_num = None

    if current_lines:
        models.append((current_num or len(models) + 1, current_lines[:]))

    if not models:
        return pdbqt_text if pdbqt_text.endswith("\n") else pdbqt_text + "\n", "No MODEL blocks found; used whole PDBQT."

    selected = next((block for number, block in models if number == pose_index), None)
    warning = ""
    if selected is None:
        number, selected = models[0]
        warning = f"Requested Vina pose {pose_index} not found; used MODEL {number}."
    return "\n".join(selected) + "\n", warning


def compose_pose_block(row: Dict[str, str], block_text: str, model_index: int) -> str:
    variant = row.get("LigandVariant", "")
    base = row.get("LigandBase", "")
    state = row.get("StateTag", "")
    source_pose = str(row.get("Pose", "")).strip() or str(model_index)
    source_outfile = str(row.get("OutFile", "")).strip()
    affinity = str(row.get("Binding_Affinity", "")).strip()

    block_lines = block_text.splitlines()
    payload_lines = [line for line in block_lines if not line.startswith("MODEL") and not line.startswith("ENDMDL")]

    header_lines = [
        f"MODEL {model_index}",
        f"REMARK COMPACTED_BASE_LIGAND: {base}",
        f"REMARK COMPACTED_VARIANT: {variant}",
        f"REMARK COMPACTED_STATE: {state}",
        f"REMARK COMPACTED_SOURCE_POSE: {source_pose}",
        f"REMARK COMPACTED_SOURCE_OUTFILE: {source_outfile}",
    ]
    if affinity:
        header_lines.append(f"REMARK COMPACTED_SOURCE_AFFINITY: {affinity}")

    return "\n".join(header_lines + payload_lines + ["ENDMDL"]) + "\n"


def build_combined_pose_text(selected_rows: Sequence[Dict[str, str]]) -> tuple[str, List[str], List[str]]:
    warnings: List[str] = []
    variants_included: List[str] = []
    blocks: List[str] = []

    if not selected_rows:
        return "", warnings, variants_included

    base_name = selected_rows[0].get("LigandBase", "Ligand")
    blocks.append(f"REMARK Name = {base_name}\n")

    for model_index, row in enumerate(selected_rows, start=1):
        out_path = Path(row.get("OutFile", "")).expanduser()
        if not out_path.is_absolute():
            out_path = out_path.resolve()
        if not out_path.exists():
            warnings.append(f"Missing OutFile for {row.get('LigandVariant', '')}: {out_path}")
            continue
        raw_text = out_path.read_text(encoding="utf-8", errors="ignore")
        pose_index = int(float(row.get("Pose", "1") or 1))
        pose_text, warning = extract_selected_pose_text(raw_text, pose_index)
        if warning:
            warnings.append(f"{row.get('LigandVariant', '')}: {warning}")
        blocks.append(compose_pose_block(row, pose_text, model_index))
        variants_included.append(row.get("LigandVariant", ""))

    return "".join(blocks), warnings, variants_included


def decorate_compacted_viewer_html(viewer_html: str) -> str:
    viewer_html = viewer_html.replace(
        ".pose-rmsd{font:12px var(--mono);color:var(--txt-muted)}",
        ".pose-rmsd{font:12px var(--mono);color:var(--txt-muted)}\n.pose-sub{margin-top:4px;font:11px var(--mono);color:var(--accent);word-break:break-word}",
        1,
    )
    viewer_html = viewer_html.replace(
        """      if (line.startsWith("REMARK VINA RESULT:")) {
        const p = line.replace("REMARK VINA RESULT:","").trim().split(/\\s+/);
        cur.score=parseFloat(p[0]); cur.rmsd_lb=parseFloat(p[1]); cur.rmsd_ub=parseFloat(p[2]);
      }""",
        """      if (line.startsWith("REMARK VINA RESULT:")) {
        const p = line.replace("REMARK VINA RESULT:","").trim().split(/\\s+/);
        cur.score=parseFloat(p[0]); cur.rmsd_lb=parseFloat(p[1]); cur.rmsd_ub=parseFloat(p[2]);
      }
      if (line.startsWith("REMARK COMPACTED_VARIANT:")) {
        cur.variant = line.replace("REMARK COMPACTED_VARIANT:","").trim();
      }
      if (line.startsWith("REMARK COMPACTED_STATE:")) {
        cur.state = line.replace("REMARK COMPACTED_STATE:","").trim();
      }
      if (line.startsWith("REMARK COMPACTED_SOURCE_POSE:")) {
        cur.sourcePose = line.replace("REMARK COMPACTED_SOURCE_POSE:","").trim();
      }""",
        1,
    )
    viewer_html = viewer_html.replace(
        """          <div class="pose-rmsd">lb ${pose.rmsd_lb.toFixed(2)} / ub ${pose.rmsd_ub.toFixed(2)} Å</div>""",
        """          <div class="pose-rmsd">lb ${pose.rmsd_lb.toFixed(2)} / ub ${pose.rmsd_ub.toFixed(2)} Å</div>
          <div class="pose-sub">${pose.variant || `Pose ${i+1}`} ${pose.state ? `· ${pose.state}` : ""}</div>""",
        1,
    )
    viewer_html = viewer_html.replace(
        '    document.getElementById("hud-pose").textContent  = `Pose ${i+1} of ${poses.length}`;',
        '    document.getElementById("hud-pose").textContent  = poses[i].variant ? `Pose ${i+1} of ${poses.length} · ${poses[i].variant}` : `Pose ${i+1} of ${poses.length}`;',
        1,
    )
    return viewer_html


def build_index_html(page_title: str, entries: List[Dict[str, Any]]) -> str:
    cards: List[str] = []
    for entry in entries:
        variants_preview = ", ".join(entry.get("variants_included", [])[:3])
        extra = ""
        if len(entry.get("variants_included", [])) > 3:
            extra = f" +{len(entry['variants_included']) - 3} more"
        cards.append(
            f"""
<article class="card">
  <div class="eyebrow">{entry["receptor"]}</div>
  <h2>{entry["ligand"]}</h2>
  <p>Best affinity: <strong>{entry["best_affinity"]}</strong> kcal/mol</p>
  <p>{entry["pose_count"]} compacted poses across {entry["variant_count"]} variant(s)</p>
  <p class="meta">{variants_preview}{extra}</p>
  <a href="{entry["viewer_file"]}">Open viewer</a>
</article>
"""
        )
    joined_cards = "\n".join(cards) or "<p>No viewer entries were generated.</p>"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{page_title}</title>
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
.card p {{ margin: 0 0 10px; color: #4d5b69; }}
.card .meta {{ font-family: "IBM Plex Mono", monospace; font-size: 12px; color: #006d77; }}
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
  <div class="eyebrow">Compacted Docking HTML Visualization</div>
  <h1>{page_title}</h1>
  <p class="lede">Each viewer below groups docking poses by original base ligand, then pools the top-ranked poses across conformers, protomers, and tautomers into a single compact review panel.</p>
  <section class="grid">
    {joined_cards}
  </section>
  <div class="credit">
    Adapted from <a href="{BASE.UPSTREAM_REPO_URL}" target="_blank" rel="noreferrer">py-VinaScope-Docking-Viewer</a> and
    <a href="{BASE.UPSTREAM_VIEWER_URL}" target="_blank" rel="noreferrer">VinaDock-Viz</a>.<br>
    Original source acknowledged to {BASE.UPSTREAM_AUTHOR}<br>
    {BASE.UPSTREAM_AUTHOR_AFFILIATION}
  </div>
</main>
</body>
</html>
"""


def build_project(
    csv_path: Path,
    outdir: Optional[Path] = None,
    receptor_roots: Optional[Iterable[Path]] = None,
    top_ligands: int = 25,
    top_poses: int = 25,
    project_name: str = "Docking_HTML_Viz_Project_Compacted",
    page_title: str = "Compacted Docking Visualization Project",
) -> Dict[str, Any]:
    csv_path = csv_path.resolve()
    base_dir = csv_path.parent
    output_root = (outdir or base_dir).resolve()
    project_dir = output_root / BASE.safe_slug(project_name)
    BASE.load_vinascope_generator_template()
    if project_dir.exists():
        shutil.rmtree(project_dir)
    viewers_dir = project_dir / "viewers"
    inputs_dir = project_dir / "inputs"
    project_dir.mkdir(parents=True, exist_ok=True)
    viewers_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)

    rows = load_rows(csv_path)
    groups = select_compacted_groups(rows, top_ligands=top_ligands, top_poses=top_poses)
    print_selection_summary(groups, top_ligands=top_ligands, top_poses=top_poses)

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
    receptor_copy_map: Dict[str, str] = {}
    used_input_names: set[str] = set()

    for group in groups:
        receptor_name = str(group["receptor"])
        ligand_base = str(group["ligand_base"])
        selected_rows = list(group["selected_rows"])
        receptor_file = BASE.find_receptor_file(receptor_name, search_roots)

        if receptor_file is None:
            missing_entries.append(
                {
                    "receptor": receptor_name,
                    "ligand": ligand_base,
                    "reason": "missing_receptor",
                }
            )
            continue

        receptor_text = receptor_file.read_text(encoding="utf-8", errors="ignore")
        ligand_text, build_warnings, variants_included = build_combined_pose_text(selected_rows)
        if not ligand_text.strip():
            missing_entries.append(
                {
                    "receptor": receptor_name,
                    "ligand": ligand_base,
                    "reason": "missing_pose_data",
                    "warnings": build_warnings,
                }
            )
            continue

        best_affinity = f"{min(safe_float(row.get('Binding_Affinity', '')) for row in selected_rows):.2f}"
        slug = BASE.safe_slug(f"{Path(receptor_name).stem}__{ligand_base}")

        receptor_key = str(receptor_file.resolve())
        receptor_copy_rel = receptor_copy_map.get(receptor_key)
        if receptor_copy_rel is None:
            receptor_filename = f"{BASE.safe_slug(Path(receptor_file).stem)}_receptor{receptor_file.suffix or '.pdbqt'}"
            receptor_copy = unique_output_path(inputs_dir, receptor_filename, used_input_names)
            receptor_copy.write_text(receptor_text, encoding="utf-8")
            receptor_copy_rel = str(receptor_copy.relative_to(project_dir))
            receptor_copy_map[receptor_key] = receptor_copy_rel

        ligand_copy = inputs_dir / f"{slug}_poses.pdbqt"
        viewer_file = viewers_dir / f"{slug}.html"

        ligand_copy.write_text(ligand_text, encoding="utf-8")
        viewer_html = BASE.build_viewer_html(
            page_title=page_title,
            receptor_label=Path(receptor_file).name,
            ligand_label=ligand_base,
            receptor_text=receptor_text,
            ligand_text=ligand_text,
            best_affinity=best_affinity,
            source_outfile="; ".join(str(Path(row.get("OutFile", "")).name) for row in selected_rows),
        )
        viewer_file.write_text(decorate_compacted_viewer_html(viewer_html), encoding="utf-8")

        manifest_entries.append(
            {
                "receptor": receptor_name,
                "ligand": ligand_base,
                "ligand_variant": ligand_base,
                "best_affinity": best_affinity,
                "pose_count": len(variants_included),
                "variant_count": len(sorted(set(variants_included))),
                "variants_included": sorted(dict.fromkeys(variants_included)),
                "source_outfiles": sorted(dict.fromkeys(str(row.get("OutFile", "")) for row in selected_rows)),
                "warnings": build_warnings,
                "receptor_file": str(receptor_file),
                "viewer_file": str(viewer_file.relative_to(project_dir)),
                "receptor_copy": receptor_copy_rel,
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
        "grouping_mode": "base_ligand_compacted",
        "attribution": {
            "repo": BASE.UPSTREAM_REPO_URL,
            "viewer": BASE.UPSTREAM_VIEWER_URL,
            "author": BASE.UPSTREAM_AUTHOR,
            "affiliation": BASE.UPSTREAM_AUTHOR_AFFILIATION,
        },
    }

    (project_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (project_dir / "UPSTREAM_LICENSE.txt").write_text(BASE.UPSTREAM_LICENSE_TEXT, encoding="utf-8")
    (project_dir / "index.html").write_text(build_index_html(page_title, manifest_entries), encoding="utf-8")
    zip_path = BASE.write_zip(project_dir)
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
