#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import importlib.util
import json
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


def _load_module(filename: str, module_name: str):
    module_path = Path(__file__).resolve().with_name(filename)
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader is not None
    spec.loader.exec_module(module)
    return module


COMPACT = _load_module("5_CompactedHTMLViz.py", "compacted_html_base")
PYMOL = _load_module("5C_BuildPymolSesh.py", "pymol_reconstruction_base")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compacted HTML project that renders corrected SDF bond orders."
    )
    parser.add_argument("--csv", help="Sorted docking CSV produced by 4_ParseScores.py or 4C_ConcatenateScores.py")
    parser.add_argument("--outdir", help="Output directory for the generated visualization project")
    parser.add_argument(
        "--receptor-roots",
        nargs="+",
        default=["Receptors", "Receptors_PDBQT", "."],
        help="Folders searched for receptor PDBQT files",
    )
    parser.add_argument("--top-ligands", type=int, default=25, help="Top base ligands per receptor to include")
    parser.add_argument("--top-poses", type=int, default=25, help="Top docking poses retained per base ligand")
    parser.add_argument(
        "--project-name",
        default="Docking_HTML_Viz_Project_Compacted_SDF",
        help="Output project folder name",
    )
    parser.add_argument(
        "--page-title",
        default="Compacted SDF Docking Visualization Project",
        help="Display title embedded in the HTML pages",
    )
    parser.add_argument("--obabel-bin", help="Path to obabel executable; default=OBABEL_BIN or PATH")
    parser.add_argument(
        "--hydrogen-mode",
        choices=["none", "nonpolar", "all"],
        default="none",
        help="Hydrogen handling for corrected SDF output. Default none matches the PyMOL workflow.",
    )
    parser.add_argument(
        "--allow-reference-fallback",
        action="store_true",
        help="If exact generated SDF is missing, allow base/original SDF fallback.",
    )
    parser.add_argument(
        "--include-previous-outputs-as-reference",
        action="store_true",
        help="Also search previous 11_AA_* output folders for reference SDFs.",
    )
    parser.add_argument(
        "--min-mcs-fraction",
        type=float,
        default=0.80,
        help="Minimum heavy-atom MCS fraction accepted for corrected-bond reconstruction.",
    )
    return parser.parse_args()


def resolve_args(args: argparse.Namespace) -> argparse.Namespace:
    if args.csv:
        return args

    print("\n5_COMPACTED_SDF_HTML interactive mode")
    print("No --csv flag was provided, so we'll gather the inputs step by step.")

    csv_path = COMPACT.BASE.prompt_csv_path()
    args.csv = str(csv_path)

    suggested_outdir = str(csv_path.resolve().parent)
    args.outdir = COMPACT.BASE.prompt_text("Output directory", suggested_outdir)

    roots_default = " ".join(COMPACT.BASE.default_receptor_roots_for(csv_path))
    receptor_roots_raw = COMPACT.BASE.prompt_text("Receptor search roots (space-separated)", roots_default)
    args.receptor_roots = receptor_roots_raw.split()

    args.top_ligands = COMPACT.BASE.prompt_int("Top base ligands per receptor", args.top_ligands)
    args.top_poses = COMPACT.BASE.prompt_int("Top poses per base ligand", args.top_poses)
    args.project_name = COMPACT.BASE.prompt_text("Project folder name", args.project_name)
    args.page_title = COMPACT.BASE.prompt_text("Viewer page title", args.page_title)
    return args


def b64encode_text(text: str) -> str:
    return base64.b64encode(text.encode("utf-8")).decode("ascii")


def build_corrected_pose_bundle(
    selected_rows: List[Dict[str, str]],
    resolver: Any,
    obabel_bin: str,
    hydrogen_mode: str,
    min_mcs_fraction: float,
    allow_reference_fallback: bool,
) -> Dict[str, Any]:
    if not getattr(PYMOL, "HAVE_RDKIT", False):
        raise SystemExit("❌ RDKit is required for corrected SDF bond-order reconstruction.")

    corrected_models: List[Dict[str, Any]] = []
    warnings: List[str] = []
    sdf_blocks: List[str] = []

    for display_index, row in enumerate(selected_rows, start=1):
        ligand_variant = row.get("LigandVariant", "")
        ligand_base = row.get("LigandBase", "")
        state_tag = row.get("StateTag", "")
        pose_index = int(COMPACT.safe_float(row.get("Pose", "1") or 1, 1))
        affinity = COMPACT.safe_float(row.get("Binding_Affinity", ""))
        outfile = Path(row.get("OutFile", "")).expanduser()
        if not outfile.is_absolute():
            outfile = outfile.resolve()
        if not outfile.exists():
            warnings.append(f"{ligand_variant}: missing OutFile {outfile}")
            continue

        with tempfile.TemporaryDirectory(prefix="compact_sdf_pose_") as tmpdir:
            selected_pose_pdbqt = Path(tmpdir) / f"{COMPACT.BASE.safe_slug(ligand_variant)}_pose{pose_index}.pdbqt"
            try:
                _, pose_warning = PYMOL.extract_selected_pose_pdbqt(outfile, pose_index, selected_pose_pdbqt)
            except Exception as exc:
                warnings.append(f"{ligand_variant}: pose extraction failed ({exc})")
                continue

            resolution = resolver.resolve(row, allow_fallback=allow_reference_fallback)
            reference_sdf = resolution.get("reference_sdf")
            ref_mol = PYMOL.load_reference_mol(reference_sdf)

            try:
                coords_mol = PYMOL.pdbqt_to_mol(selected_pose_pdbqt, obabel_bin)
            except Exception as exc:
                warnings.append(f"{ligand_variant}: Open Babel conversion failed ({exc})")
                continue

            fitted_mol = None
            fit_warning = pose_warning or ""
            if ref_mol is not None:
                fitted_mol, fit_status, _mcs_size, warning = PYMOL.rigid_fit_by_mcs(
                    coords_mol,
                    ref_mol,
                    min_mcs_fraction=min_mcs_fraction,
                )
                if warning:
                    fit_warning = (fit_warning + " " + warning).strip()
            else:
                fit_status = "missing_reference_sdf"

            if fitted_mol is None:
                fitted_mol = PYMOL.Chem.RemoveHs(PYMOL.Chem.Mol(coords_mol), sanitize=False)
                if fit_status == "missing_reference_sdf":
                    fit_warning = (fit_warning + " No reference SDF was available; used PDBQT coordinates only.").strip()
                else:
                    fit_warning = (fit_warning + " Falling back to PDBQT coordinates only.").strip()

            final_mol = PYMOL.apply_hydrogen_mode(fitted_mol, hydrogen_mode)
            molblock = PYMOL.Chem.MolToMolBlock(final_mol, kekulize=False)
            sdf_block = molblock
            if not sdf_block.rstrip().endswith("$$$$"):
                sdf_block = sdf_block.rstrip() + "\n$$$$\n"
            sdf_blocks.append(sdf_block)

            corrected_models.append(
                {
                    "pose_index": display_index,
                    "variant": ligand_variant,
                    "ligand_base": ligand_base,
                    "state": state_tag,
                    "source_pose": pose_index,
                    "score": affinity,
                    "fit_status": fit_status,
                    "fit_warning": fit_warning,
                    "reference_sdf": str(reference_sdf or ""),
                    "molblock": molblock,
                }
            )

    return {
        "models": corrected_models,
        "warnings": warnings,
        "combined_sdf_text": "".join(sdf_blocks),
    }


def decorate_sdf_viewer_html(viewer_html: str, corrected_models: List[Dict[str, Any]]) -> str:
    encoded_json = b64encode_text(json.dumps(corrected_models))
    viewer_html = viewer_html.replace(
        '  const poseData     = b64ToStr(POSE_B64);\n  const poses        = parsePoses(poseData);\n',
        '  const poseData     = b64ToStr(POSE_B64);\n  const poses        = parsePoses(poseData);\n'
        f'  const correctedPoseData = JSON.parse(b64ToStr("{encoded_json}"));\n',
        1,
    )
    viewer_html = viewer_html.replace(
        '  const poseMdls   = poses.map(p=>viewer.addModel(poseToPDB(p),"pdb"));',
        '  const poseMdls   = poses.map((p,i)=>{\n'
        '    const corrected = correctedPoseData[i];\n'
        '    if (corrected && corrected.molblock) return viewer.addModel(corrected.molblock,"sdf");\n'
        '    return viewer.addModel(poseToPDB(p),"pdb");\n'
        '  });',
        1,
    )
    viewer_html = viewer_html.replace(
        '          <div class="pose-sub">${pose.variant || `Pose ${i+1}`} ${pose.state ? `· ${pose.state}` : ""}</div>',
        '          <div class="pose-sub">${pose.variant || correctedPoseData[i]?.variant || `Pose ${i+1}`} ${(pose.state || correctedPoseData[i]?.state) ? `· ${pose.state || correctedPoseData[i]?.state}` : ""}</div>',
        1,
    )
    viewer_html = viewer_html.replace(
        '    document.getElementById("hud-pose").textContent  = poses[i].variant ? `Pose ${i+1} of ${poses.length} · ${poses[i].variant}` : `Pose ${i+1} of ${poses.length}`;',
        '    const activeVariant = poses[i].variant || correctedPoseData[i]?.variant;\n'
        '    document.getElementById("hud-pose").textContent  = activeVariant ? `Pose ${i+1} of ${poses.length} · ${activeVariant}` : `Pose ${i+1} of ${poses.length}`;',
        1,
    )
    viewer_html = viewer_html.replace(
        '      toast("Viewer ready");',
        '      toast("Viewer ready (corrected SDF bond orders enabled)");',
        1,
    )
    return viewer_html


def build_project(
    csv_path: Path,
    outdir: Optional[Path] = None,
    receptor_roots: Optional[Iterable[Path]] = None,
    top_ligands: int = 25,
    top_poses: int = 25,
    project_name: str = "Docking_HTML_Viz_Project_Compacted_SDF",
    page_title: str = "Compacted SDF Docking Visualization Project",
    obabel_bin_cli: Optional[str] = None,
    hydrogen_mode: str = "none",
    allow_reference_fallback: bool = False,
    include_previous_outputs_as_reference: bool = False,
    min_mcs_fraction: float = 0.80,
) -> Dict[str, Any]:
    csv_path = csv_path.resolve()
    base_dir = csv_path.parent
    output_root = (outdir or base_dir).resolve()
    project_dir = output_root / COMPACT.BASE.safe_slug(project_name)
    COMPACT.BASE.load_vinascope_generator_template()
    if project_dir.exists():
        shutil.rmtree(project_dir)
    viewers_dir = project_dir / "viewers"
    inputs_dir = project_dir / "inputs"
    project_dir.mkdir(parents=True, exist_ok=True)
    viewers_dir.mkdir(exist_ok=True)
    inputs_dir.mkdir(exist_ok=True)

    rows = COMPACT.load_rows(csv_path)
    groups = COMPACT.select_compacted_groups(rows, top_ligands=top_ligands, top_poses=top_poses)
    COMPACT.print_selection_summary(groups, top_ligands=top_ligands, top_poses=top_poses)

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

    obabel_bin = PYMOL.obabel_path(obabel_bin_cli)
    resolver = PYMOL.ReferenceResolver(
        base_dir.resolve(),
        include_previous_outputs_as_reference=include_previous_outputs_as_reference,
    )

    manifest_entries: List[Dict[str, Any]] = []
    missing_entries: List[Dict[str, Any]] = []
    receptor_copy_map: Dict[str, str] = {}
    used_input_names: set[str] = set()

    for group in groups:
        receptor_name = str(group["receptor"])
        ligand_base = str(group["ligand_base"])
        selected_rows = list(group["selected_rows"])
        receptor_file = COMPACT.BASE.find_receptor_file(receptor_name, search_roots)

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
        ligand_text, build_warnings, variants_included = COMPACT.build_combined_pose_text(selected_rows)
        corrected_bundle = build_corrected_pose_bundle(
            selected_rows=selected_rows,
            resolver=resolver,
            obabel_bin=obabel_bin,
            hydrogen_mode=hydrogen_mode,
            min_mcs_fraction=min_mcs_fraction,
            allow_reference_fallback=allow_reference_fallback,
        )
        build_warnings = build_warnings + corrected_bundle["warnings"]
        corrected_models = corrected_bundle["models"]

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
        if not corrected_models:
            missing_entries.append(
                {
                    "receptor": receptor_name,
                    "ligand": ligand_base,
                    "reason": "missing_corrected_sdf_models",
                    "warnings": build_warnings,
                }
            )
            continue

        best_affinity = f"{min(COMPACT.safe_float(row.get('Binding_Affinity', '')) for row in selected_rows):.2f}"
        slug = COMPACT.BASE.safe_slug(f"{Path(receptor_name).stem}__{ligand_base}")

        receptor_key = str(receptor_file.resolve())
        receptor_copy_rel = receptor_copy_map.get(receptor_key)
        if receptor_copy_rel is None:
            receptor_filename = f"{COMPACT.BASE.safe_slug(Path(receptor_file).stem)}_receptor{receptor_file.suffix or '.pdbqt'}"
            receptor_copy = COMPACT.unique_output_path(inputs_dir, receptor_filename, used_input_names)
            receptor_copy.write_text(receptor_text, encoding="utf-8")
            receptor_copy_rel = str(receptor_copy.relative_to(project_dir))
            receptor_copy_map[receptor_key] = receptor_copy_rel

        ligand_copy = inputs_dir / f"{slug}_poses.pdbqt"
        corrected_sdf_copy = inputs_dir / f"{slug}_corrected_poses.sdf"
        viewer_file = viewers_dir / f"{slug}.html"

        ligand_copy.write_text(ligand_text, encoding="utf-8")
        corrected_sdf_copy.write_text(corrected_bundle["combined_sdf_text"], encoding="utf-8")

        viewer_html = COMPACT.BASE.build_viewer_html(
            page_title=page_title,
            receptor_label=Path(receptor_file).name,
            ligand_label=ligand_base,
            receptor_text=receptor_text,
            ligand_text=ligand_text,
            best_affinity=best_affinity,
            source_outfile="; ".join(str(Path(row.get("OutFile", "")).name) for row in selected_rows),
        )
        viewer_html = COMPACT.decorate_compacted_viewer_html(viewer_html)
        viewer_html = decorate_sdf_viewer_html(viewer_html, corrected_models)
        viewer_file.write_text(viewer_html, encoding="utf-8")

        manifest_entries.append(
            {
                "receptor": receptor_name,
                "ligand": ligand_base,
                "ligand_variant": ligand_base,
                "best_affinity": best_affinity,
                "pose_count": len(corrected_models),
                "variant_count": len(sorted(set(variants_included))),
                "variants_included": sorted(dict.fromkeys(variants_included)),
                "source_outfiles": sorted(dict.fromkeys(str(row.get("OutFile", "")) for row in selected_rows)),
                "warnings": build_warnings,
                "bond_mode": "corrected_sdf_from_reference_fit",
                "receptor_file": str(receptor_file),
                "viewer_file": str(viewer_file.relative_to(project_dir)),
                "receptor_copy": receptor_copy_rel,
                "pose_copy": str(ligand_copy.relative_to(project_dir)),
                "corrected_sdf_copy": str(corrected_sdf_copy.relative_to(project_dir)),
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
        "bond_mode": "corrected_sdf_from_reference_fit",
        "hydrogen_mode": hydrogen_mode,
        "attribution": {
            "repo": COMPACT.BASE.UPSTREAM_REPO_URL,
            "viewer": COMPACT.BASE.UPSTREAM_VIEWER_URL,
            "author": COMPACT.BASE.UPSTREAM_AUTHOR,
            "affiliation": COMPACT.BASE.UPSTREAM_AUTHOR_AFFILIATION,
        },
    }

    (project_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (project_dir / "UPSTREAM_LICENSE.txt").write_text(COMPACT.BASE.UPSTREAM_LICENSE_TEXT, encoding="utf-8")
    (project_dir / "index.html").write_text(COMPACT.build_index_html(page_title, manifest_entries), encoding="utf-8")
    zip_path = COMPACT.BASE.write_zip(project_dir)
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
        obabel_bin_cli=args.obabel_bin,
        hydrogen_mode=args.hydrogen_mode,
        allow_reference_fallback=args.allow_reference_fallback,
        include_previous_outputs_as_reference=args.include_previous_outputs_as_reference,
        min_mcs_fraction=args.min_mcs_fraction,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
