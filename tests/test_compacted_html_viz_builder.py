import importlib.util
import io
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


class CompactedHtmlVizBuilderTests(unittest.TestCase):
    def test_print_selection_summary_uses_rank_order_with_scores(self):
        module = load_script_module("5_CompactedHTMLViz.py", "compacted_html_viz_ranked_summary")
        groups = [
            {
                "receptor": "9ba9",
                "ligand_base": "CONIVAPTAN",
                "selected_rows": [
                    {"Binding_Affinity": "-9.344"},
                    {"Binding_Affinity": "-9.100"},
                ],
            },
            {
                "receptor": "9ba9",
                "ligand_base": "DIHYDROERGOTAMINE_MESYLATE",
                "selected_rows": [
                    {"Binding_Affinity": "-9.377"},
                ],
            },
        ]
        stream = io.StringIO()
        with mock.patch("sys.stdout", stream):
            module.print_selection_summary(groups, top_ligands=50, top_poses=15)

        output = stream.getvalue()
        first = output.index("[1] DIHYDROERGOTAMINE_MESYLATE")
        second = output.index("[2] CONIVAPTAN")
        self.assertLess(first, second)
        self.assertIn("best=-9.377  poses=1", output)
        self.assertIn("best=-9.344  poses=2", output)

    def test_build_project_compacts_multiple_variants_into_one_viewer(self):
        module = load_script_module("5_CompactedHTMLViz.py", "compacted_html_viz_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            (receptors / "recA.pdbqt").write_text(
                "ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n",
                encoding="utf-8",
            )

            rows = []
            for variant, score in [
                ("LigA_p01_t01_c001", "-9.1"),
                ("LigA_p01_t02_c001", "-8.8"),
                ("LigA_p01_t03_c001", "-8.2"),
            ]:
                out_dir = root / "Docking_Results_demo" / "recA" / variant
                out_dir.mkdir(parents=True)
                pose_file = out_dir / "out.pdbqt"
                pose_file.write_text(
                    f"MODEL 1\nREMARK VINA RESULT: {score} 0.0 0.0\nATOM      1  C   UNL A   1      12.000  10.000   8.000  0.00  0.00    C\nENDMDL\n",
                    encoding="utf-8",
                )
                rows.append(f"recA,LigA,LigA,{variant},1,{score},{pose_file}")

            csv_path.write_text(
                "Receptor,Ligand,LigandBase,LigandVariant,Pose,Binding_Affinity,OutFile\n"
                + "\n".join(rows)
                + "\n",
                encoding="utf-8",
            )

            manifest = module.build_project(
                csv_path=csv_path,
                receptor_roots=[receptors],
                top_ligands=1,
                top_poses=2,
                project_name="CompactViz",
                page_title="Compact Viz",
            )

            self.assertEqual(manifest["entry_count"], 1)
            entry = manifest["entries"][0]
            self.assertEqual(entry["ligand"], "LigA")
            self.assertEqual(entry["pose_count"], 2)
            self.assertEqual(entry["variant_count"], 2)
            self.assertEqual(
                entry["variants_included"],
                ["LigA_p01_t01_c001", "LigA_p01_t02_c001"],
            )

            pose_copy = Path(manifest["project_dir"]) / entry["pose_copy"]
            pose_text = pose_copy.read_text(encoding="utf-8")
            self.assertIn("REMARK Name = LigA", pose_text)
            self.assertIn("REMARK COMPACTED_VARIANT: LigA_p01_t01_c001", pose_text)
            self.assertIn("REMARK COMPACTED_VARIANT: LigA_p01_t02_c001", pose_text)
            self.assertNotIn("LigA_p01_t03_c001", pose_text)
            self.assertEqual(pose_text.count("MODEL "), 2)

            viewer_path = Path(manifest["project_dir"]) / entry["viewer_file"]
            viewer_html = viewer_path.read_text(encoding="utf-8")
            self.assertIn("pose.variant", viewer_html)
            self.assertIn("REMARK COMPACTED_VARIANT:", pose_text)

            parsed_manifest = json.loads((Path(manifest["project_dir"]) / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(parsed_manifest["grouping_mode"], "base_ligand_compacted")


if __name__ == "__main__":
    unittest.main()
