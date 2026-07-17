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


class CompactedSdfHtmlBuilderTests(unittest.TestCase):
    def test_build_project_uses_corrected_sdf_payload(self):
        module = load_script_module("5_COMPACTED_SDF_HTML.py", "compacted_sdf_html_module")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            csv_path = root / "scores.csv"
            receptors = root / "Receptors"
            receptors.mkdir()
            (receptors / "recA.pdbqt").write_text(
                "ATOM      1  N   GLY A   1      11.104  13.207  10.451  1.00 20.00           N\n",
                encoding="utf-8",
            )

            out_dir = root / "Docking_Results_demo" / "recA" / "LigA_p01_t01_c001"
            out_dir.mkdir(parents=True)
            pose_file = out_dir / "out.pdbqt"
            pose_file.write_text(
                "MODEL 1\nREMARK VINA RESULT: -9.1 0.0 0.0\nATOM      1  C   UNL A   1      12.000  10.000   8.000  0.00  0.00    C\nENDMDL\n",
                encoding="utf-8",
            )
            csv_path.write_text(
                "Receptor,Ligand,LigandBase,LigandVariant,Pose,Binding_Affinity,OutFile\n"
                f"recA,LigA,LigA,LigA_p01_t01_c001,1,-9.1,{pose_file}\n",
                encoding="utf-8",
            )

            fake_bundle = {
                "models": [
                    {
                        "pose_index": 1,
                        "variant": "LigA_p01_t01_c001",
                        "ligand_base": "LigA",
                        "state": "p01_t01",
                        "source_pose": 1,
                        "score": -9.1,
                        "fit_status": "og_mcs_rigid_fit",
                        "fit_warning": "",
                        "reference_sdf": "/tmp/LigA.sdf",
                        "molblock": "LigA\n  Mock\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n",
                    }
                ],
                "warnings": [],
                "combined_sdf_text": "LigA\n  Mock\n\n  1  0  0  0  0  0            999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
            }

            with mock.patch.object(module.PYMOL, "obabel_path", return_value="/usr/bin/obabel"), mock.patch.object(
                module.PYMOL,
                "ReferenceResolver",
                return_value=object(),
            ), mock.patch.object(
                module,
                "build_corrected_pose_bundle",
                return_value=fake_bundle,
            ):
                manifest = module.build_project(
                    csv_path=csv_path,
                    receptor_roots=[receptors],
                    top_ligands=1,
                    top_poses=1,
                    project_name="CompactSDFViz",
                    page_title="Compact SDF Viz",
                )

            self.assertEqual(manifest["entry_count"], 1)
            entry = manifest["entries"][0]
            self.assertEqual(entry["bond_mode"], "corrected_sdf_from_reference_fit")
            corrected_sdf = Path(manifest["project_dir"]) / entry["corrected_sdf_copy"]
            self.assertTrue(corrected_sdf.exists())
            self.assertIn("$$$$", corrected_sdf.read_text(encoding="utf-8"))

            viewer_html = (Path(manifest["project_dir"]) / entry["viewer_file"]).read_text(encoding="utf-8")
            self.assertIn("correctedPoseData", viewer_html)
            self.assertIn('viewer.addModel(corrected.molblock,"sdf")', viewer_html)

            parsed_manifest = json.loads((Path(manifest["project_dir"]) / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(parsed_manifest["bond_mode"], "corrected_sdf_from_reference_fit")


if __name__ == "__main__":
    unittest.main()
