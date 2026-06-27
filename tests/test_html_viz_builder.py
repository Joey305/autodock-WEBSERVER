import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


class HtmlVizBuilderTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
