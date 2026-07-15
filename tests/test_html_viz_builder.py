import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class HtmlVizBuilderTests(unittest.TestCase):
    def test_resolve_args_prompts_for_missing_csv(self):
        module = load_script_module("5_BuidlHTMLViz.py", "html_viz_builder_prompt_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "demo_vina_docking_scores_sorted.csv"
            csv_path.write_text("Receptor,Ligand,Pose,Binding_Affinity,OutFile\n", encoding="utf-8")

            args = module.argparse.Namespace(
                csv=None,
                outdir=None,
                receptor_roots=["Receptors", "Receptors_PDBQT", "."],
                top_ligands=25,
                top_poses=3,
                project_name="Docking_HTML_Viz_Project",
                page_title="Docking Visualization Project",
            )

            responses = iter(
                [
                    "1",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ]
            )

            with mock.patch.object(module.Path, "cwd", return_value=root):
                with mock.patch("builtins.input", side_effect=lambda _: next(responses)):
                    resolved = module.resolve_args(args)

            self.assertEqual(Path(resolved.csv).resolve(), csv_path.resolve())
            self.assertEqual(Path(resolved.outdir).resolve(), root.resolve())
            self.assertEqual(resolved.receptor_roots, ["."])

    def test_build_project_writes_manifest_viewers_and_zip(self):
        module = load_script_module("5_BuidlHTMLViz.py", "html_viz_builder_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            receptor_file = receptors / "recA.pdbqt"
            receptor_file.write_text("ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n")

            out_dir = root / "Docking_Results_demo" / "recA" / "LigA"
            out_dir.mkdir(parents=True)
            pose_file = out_dir / "out.pdbqt"
            pose_file.write_text(
                "MODEL 1\nREMARK VINA RESULT: -9.1 0.0 0.0\nATOM      1  C   LIG A   1      12.000  10.000   8.000  0.00  0.00    C\nENDMDL\n"
                "MODEL 2\nREMARK VINA RESULT: -8.2 0.0 0.0\nATOM      1  C   LIG A   1      12.500  10.500   8.500  0.00  0.00    C\nENDMDL\n"
            )
            csv_path.write_text(
                "Receptor,Ligand,LigandBase,LigandVariant,Pose,Binding_Affinity,OutFile\n"
                f"recA,LigA,LigA,LigA_p01_t01_c001,1,-9.1,{pose_file}\n"
                f"recA,LigA,LigA,LigA_p01_t01_c001,2,-8.2,{pose_file}\n",
                encoding="utf-8",
            )

            manifest = module.build_project(
                csv_path=csv_path,
                receptor_roots=[receptors],
                top_ligands=10,
                top_poses=1,
                project_name="DemoViz",
                page_title="Demo Viz",
            )

            manifest_path = Path(manifest["project_dir"]) / "manifest.json"
            viewer_path = Path(manifest["project_dir"]) / manifest["entries"][0]["viewer_file"]
            pose_copy = Path(manifest["project_dir"]) / manifest["entries"][0]["pose_copy"]
            zip_path = Path(manifest["zip_path"])

            self.assertTrue(manifest_path.exists())
            self.assertTrue(viewer_path.exists())
            self.assertTrue(pose_copy.exists())
            self.assertTrue(zip_path.exists())
            self.assertEqual(manifest["entry_count"], 1)
            self.assertIn("py-VinaScope-Docking-Viewer", viewer_path.read_text(encoding="utf-8"))
            pose_text = pose_copy.read_text(encoding="utf-8")
            self.assertIn("MODEL 1", pose_text)
            self.assertNotIn("MODEL 2", pose_text)

            parsed = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(parsed["entries"][0]["ligand"], "LigA")

    def test_build_project_keeps_multiple_variants_for_selected_ligand_base(self):
        module = load_script_module("5_BuidlHTMLViz.py", "html_viz_builder_variants_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            (receptors / "recA.pdbqt").write_text("ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n")

            variant_paths = []
            for variant, score in [("LigA_p01_t01_c001", "-9.1"), ("LigA_p01_t02_c001", "-8.8")]:
                out_dir = root / "Docking_Results_demo" / "recA" / variant
                out_dir.mkdir(parents=True)
                pose_file = out_dir / "out.pdbqt"
                pose_file.write_text(
                    f"MODEL 1\nREMARK VINA RESULT: {score} 0.0 0.0\nATOM      1  C   UNL A   1      12.000  10.000   8.000  0.00  0.00    C\nENDMDL\n",
                    encoding="utf-8",
                )
                variant_paths.append((variant, pose_file, score))

            csv_path.write_text(
                "Receptor,Ligand,LigandBase,LigandVariant,Pose,Binding_Affinity,OutFile\n"
                + "\n".join(
                    f"recA,LigA,LigA,{variant},1,{score},{pose_file}"
                    for variant, pose_file, score in variant_paths
                )
                + "\n",
                encoding="utf-8",
            )

            manifest = module.build_project(
                csv_path=csv_path,
                receptor_roots=[receptors],
                top_ligands=1,
                top_poses=1,
                project_name="DemoViz",
                page_title="Demo Viz",
            )

            self.assertEqual(manifest["entry_count"], 2)
            variants = [entry["ligand_variant"] for entry in manifest["entries"]]
            self.assertEqual(variants, ["LigA_p01_t01_c001", "LigA_p01_t02_c001"])

    def test_build_project_reuses_one_receptor_copy_per_unique_receptor(self):
        module = load_script_module("5_BuidlHTMLViz.py", "html_viz_builder_dedupe_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            (receptors / "recA.pdbqt").write_text(
                "ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n",
                encoding="utf-8",
            )

            entries = []
            for variant, score in [("LigA_p01_t01_c001", "-9.1"), ("LigB_p01_t01_c001", "-8.8")]:
                out_dir = root / "Docking_Results_demo" / "recA" / variant
                out_dir.mkdir(parents=True)
                pose_file = out_dir / "out.pdbqt"
                pose_file.write_text(
                    f"MODEL 1\nREMARK VINA RESULT: {score} 0.0 0.0\nATOM      1  C   UNL A   1      12.000  10.000   8.000  0.00  0.00    C\nENDMDL\n",
                    encoding="utf-8",
                )
                ligand = variant.split("_", 1)[0]
                entries.append(f"recA,{ligand},{ligand},{variant},1,{score},{pose_file}")

            csv_path.write_text(
                "Receptor,Ligand,LigandBase,LigandVariant,Pose,Binding_Affinity,OutFile\n"
                + "\n".join(entries)
                + "\n",
                encoding="utf-8",
            )

            manifest = module.build_project(
                csv_path=csv_path,
                receptor_roots=[receptors],
                top_ligands=10,
                top_poses=1,
                project_name="DemoViz",
                page_title="Demo Viz",
            )

            receptor_copies = {entry["receptor_copy"] for entry in manifest["entries"]}
            self.assertEqual(len(receptor_copies), 1)
            inputs_dir = Path(manifest["project_dir"]) / "inputs"
            self.assertEqual(len(list(inputs_dir.glob("*_receptor.pdbqt"))), 1)
            self.assertEqual(len(list(inputs_dir.glob("*_poses.pdbqt"))), 2)


if __name__ == "__main__":
    unittest.main()
