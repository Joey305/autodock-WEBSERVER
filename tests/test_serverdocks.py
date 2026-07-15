import importlib.util
import os
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


class ServerDocksTests(unittest.TestCase):
    def test_discovers_arbitrary_ligand_pdbqt_directories(self):
        module = load_script_module("3B_ServerDocks.py", "serverdocks_ligands")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            lig1 = root / "chembl_phase4_approved_smiles_part_001_Ligands_PDBQT_64Poses_20260714_0754"
            lig2 = root / "my_custom_batch_output"
            rec = root / "Receptors"
            other = root / "misc_folder"
            lig1.mkdir()
            lig2.mkdir()
            rec.mkdir()
            other.mkdir()

            original_cwd = Path.cwd()
            try:
                os.chdir(root)
                found = [p.name for p in module.ligand_dirs_only()]
            finally:
                os.chdir(original_cwd)

        self.assertEqual(found[:2], [lig1.name, lig2.name])
        self.assertIn(other.name, found)

    def test_discovers_receptor_directories_by_name_priority(self):
        module = load_script_module("3B_ServerDocks.py", "serverdocks_receptors")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            rec1 = root / "Receptors"
            rec2 = root / "protein_batch_a"
            rec3 = root / "receptor_set_alt"
            lig = root / "ligand_outputs"
            rec1.mkdir()
            rec2.mkdir()
            rec3.mkdir()
            lig.mkdir()

            original_cwd = Path.cwd()
            try:
                os.chdir(root)
                found = [p.name for p in module.receptor_dirs_only()]
            finally:
                os.chdir(original_cwd)

        self.assertEqual(found[:2], ["Receptors", "receptor_set_alt"])


if __name__ == "__main__":
    unittest.main()
