import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


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

    def test_single_receptor_lsf_targets_only_that_receptor(self):
        module = load_script_module("3B_ServerDocks.py", "serverdocks_single_receptor")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            receptors = root / "Receptors"
            receptors.mkdir()
            (receptors / "b_site.pdbqt").touch()
            (receptors / "A_site.PDBQT").touch()
            (receptors / "notes.txt").touch()

            self.assertEqual(
                [path.name for path in module.receptor_files_in_dir(receptors)],
                ["A_site.PDBQT", "b_site.pdbqt"],
            )

            original_here = module.HERE
            try:
                module.HERE = root
                lsf_path = module.write_lsf(
                    jobtag="Receptors_A_site_ligands",
                    receptors="Receptors",
                    receptor_file="Receptors/A_site.PDBQT",
                    ligands="Ligands",
                    centers_csv="vina_centers.csv",
                    poses=9,
                    exhaustiveness=8,
                    min_rmsd=1.0,
                    energy_range=100.0,
                    batch_seed=123,
                    retention_mode="full_ensemble",
                    queue="gpu_cheminfo",
                    project="brd",
                    walltime="24:00",
                    workers=4,
                    mem_per_core=2000,
                    email="",
                )
            finally:
                module.HERE = original_here

            contents = lsf_path.read_text(encoding="utf-8")
            self.assertIn('--receptor-file "Receptors/A_site.PDBQT"', contents)

    def test_main_generates_one_submission_per_receptor_and_ligand_directory(self):
        module = load_script_module("3B_ServerDocks.py", "serverdocks_main_per_receptor")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            receptors = root / "Receptors"
            ligands = root / "Ligands_PDBQT"
            receptors.mkdir()
            ligands.mkdir()
            (receptors / "site_a.pdbqt").touch()
            (receptors / "site_b.pdbqt").touch()
            (ligands / "ligand_a.pdbqt").touch()
            (root / "vina_centers.csv").write_text(
                "PDB_ID,X,Y,Z\nsite_a,1,2,3\nsite_b,4,5,6\n", encoding="utf-8"
            )

            original_cwd = Path.cwd()
            original_here = module.HERE
            try:
                os.chdir(root)
                module.HERE = root
                with patch("builtins.input", side_effect=["1", "1", "1", "9", "8", "1", "2", "123", "1", "1"]):
                    module.main()
            finally:
                module.HERE = original_here
                os.chdir(original_cwd)

            generated = sorted(root.glob("run_vina_*.lsf"))
            self.assertEqual(len(generated), 2)
            self.assertEqual(
                {path.read_text(encoding="utf-8").split('--receptor-file "')[1].split('"', 1)[0] for path in generated},
                {"Receptors/site_a.pdbqt", "Receptors/site_b.pdbqt"},
            )
            submitter = (root / "submit_all_vina.sh").read_text(encoding="utf-8")
            self.assertEqual(submitter.count("bsub <"), 2)

    def test_runner_build_jobs_filters_to_requested_receptor(self):
        module = load_script_module("3_Complete_batch_docking.py", "batch_docking_filter")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            receptors = root / "Receptors"
            ligands = root / "Ligands"
            results = root / "results"
            receptors.mkdir()
            ligands.mkdir()
            (receptors / "site_a.pdbqt").touch()
            target = receptors / "site_b.pdbqt"
            target.touch()
            (ligands / "ligand_a.pdbqt").touch()

            jobs = module.build_jobs(
                str(results),
                str(ligands),
                str(receptors),
                [
                    {"PDB_ID": "site_a", "X": "1", "Y": "2", "Z": "3"},
                    {"PDB_ID": "site_b", "X": "4", "Y": "5", "Z": "6"},
                ],
                9,
                "vina",
                8,
                1.0,
                3.0,
                123,
                str(target),
            )

            self.assertEqual(len(jobs), 1)
            self.assertEqual(Path(jobs[0]["receptor_file"]), target)
            self.assertEqual(jobs[0]["pdbid"], "site_b")


if __name__ == "__main__":
    unittest.main()
