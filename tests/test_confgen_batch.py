import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ConfgenBatchTests(unittest.TestCase):
    def test_build_run_cmd_mode4_includes_csv_folder_and_state_flags(self):
        module = load_script_module("1B_confgen_batch.py", "confgen_batch_mode4")

        class Args:
            enumerate_protomers = True
            ph_min = 6.2
            ph_max = 7.8
            ph_precision = 0.25
            max_protomers = 6
            max_tautomers = 5
            max_transforms = 250

        cmd = module.build_run_cmd(
            mode="4",
            target="Ligands_CPD_demo",
            poses=32,
            workers=8,
            filetype=None,
            csv_smiles_col="canonical_smiles",
            csv_id_col="compound_id",
            args=Args(),
        )

        self.assertIn("--mode 4", cmd)
        self.assertIn('--folder "Ligands_CPD_demo"', cmd)
        self.assertIn("--smiles-col canonical_smiles", cmd)
        self.assertIn("--id-col compound_id", cmd)
        self.assertIn("--enumerate-protomers", cmd)
        self.assertIn("--max-tautomers 5", cmd)
        self.assertIn("--max-transforms 250", cmd)

    def test_discover_csv_folders_finds_only_directories_with_csv(self):
        module = load_script_module("1B_confgen_batch.py", "confgen_batch_discovery")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "Ligands_CPD_1").mkdir()
            (root / "Ligands_CPD_1" / "ligands.csv").write_text("smiles,id\n", encoding="utf-8")
            (root / "Ligands_CPD_2").mkdir()
            (root / "notes.txt").write_text("x", encoding="utf-8")

            original_cwd = Path.cwd()
            try:
                import os
                os.chdir(root)
                found = module.discover_csv_folders()
            finally:
                os.chdir(original_cwd)

        self.assertEqual(found, ["Ligands_CPD_1"])

    def test_discover_folders_finds_split_batch_dirs_and_smi_inputs(self):
        module = load_script_module("1B_confgen_batch.py", "confgen_batch_split_discovery")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "Ligands_split_001").mkdir()
            (root / "Ligands_split_001" / "ligA.sdf").write_text("demo", encoding="utf-8")
            (root / "Ligands_split_002").mkdir()
            (root / "Ligands_split_002" / "ligB.sdf").write_text("demo", encoding="utf-8")
            (root / "Smiles_batch_001").mkdir()
            (root / "Smiles_batch_001" / "lig1.smi").write_text("CCO ethanol\n", encoding="utf-8")
            (root / "empty_dir").mkdir()

            original_cwd = Path.cwd()
            try:
                import os
                os.chdir(root)
                sdf_found = module.discover_folders("sdf")
                smiles_found = module.discover_folders("smiles")
            finally:
                os.chdir(original_cwd)

        self.assertEqual(sdf_found, ["Ligands_split_001", "Ligands_split_002"])
        self.assertEqual(smiles_found, ["Smiles_batch_001"])


if __name__ == "__main__":
    unittest.main()
