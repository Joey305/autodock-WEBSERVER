import tempfile
import unittest
from pathlib import Path

from app import infer_ligand_workflow, normalize_package_mode
from packager import assemble_job_tree
from runner_templates import build_portable_runners


class PackagingModeTests(unittest.TestCase):
    RETIRED_SERVER = "4" + "_SERVERBATCH.py"
    RETIRED_PARSE = "4C" + "MULTIVINA.py"

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.ws = Path(self._tmp.name)
        (self.ws / "Receptors").mkdir()
        (self.ws / "Ligands").mkdir()
        (self.ws / "Receptors_PDBQT").mkdir()
        (self.ws / "Receptors" / "raw_receptor.pdb").write_text("ATOM raw\n")
        (self.ws / "Receptors_PDBQT" / "raw_receptor.pdbqt").write_text("RECEPTOR\n")
        (self.ws / "Ligands" / "ligand.sdf").write_text("ligand\n")
        (self.ws / "vina_centers.csv").write_text("PDB_ID,X,Y,Z,SIZE\nraw_receptor.pdbqt,1,2,3,20\n")

    def tearDown(self):
        self._tmp.cleanup()

    def test_normalize_package_mode(self):
        self.assertEqual(normalize_package_mode({"package_mode": "portable"}), "portable")
        self.assertEqual(normalize_package_mode({"package_mode": "lsf"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "joey_lsf"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "mainak_lsf"}), "mainak_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "custom_lsf"}), "custom_lsf")
        self.assertEqual(normalize_package_mode({"include_lsf": "1"}), "joey_lsf")
        self.assertEqual(normalize_package_mode({}, default_mode="lsf"), "joey_lsf")
        self.assertEqual(normalize_package_mode({"package_mode": "lsf"}, lsf_enabled=False), "portable")

    def test_infer_ligand_workflow(self):
        self.assertEqual(infer_ligand_workflow({"upload_mode": "single", "filename": "ligands.csv"}), ("1", "csv", None))
        self.assertEqual(infer_ligand_workflow({"upload_mode": "single", "filename": "ligands.sdf"}), ("3", "sdf", "Ligands/ligands.sdf"))
        self.assertEqual(infer_ligand_workflow({"upload_mode": "single", "filename": "ligands.smiles"}), ("2", "smiles", None))
        self.assertEqual(infer_ligand_workflow({"upload_mode": "zip", "filename": "Ligands.zip"}), ("2", "sdf", None))

    def test_portable_package_stages_public_runtime_only(self):
        jobroot, warnings = assemble_job_tree(self.ws, self.ws / "Receptors", self.ws / "Ligands", package_mode="portable")
        build_portable_runners(jobroot)

        self.assertTrue((jobroot / "1_ConformerGeneration.py").exists())
        self.assertTrue((jobroot / "3_Complete_batch_docking.py").exists())
        self.assertTrue((jobroot / "4_ParseScores.py").exists())
        self.assertTrue((jobroot / "4C_ConcatenateScores.py").exists())
        self.assertTrue((jobroot / "7_Graphs.py").exists())
        self.assertTrue((jobroot / "5_CompactedHTMLViz.py").exists())
        self.assertTrue((jobroot / "5_COMPACTED_SDF_HTML.py").exists())
        self.assertTrue((jobroot / "create_vina_env.sh").exists())
        self.assertTrue((jobroot / "docking.yaml").exists())
        self.assertTrue((jobroot / "ligand_naming.py").exists())
        self.assertTrue((jobroot / "ligand_manifest.py").exists())
        self.assertTrue((jobroot / "run_confgen_local.sh").exists())
        self.assertTrue((jobroot / "run_vina_local.sh").exists())
        self.assertTrue((jobroot / "run_all_local.sh").exists())
        self.assertTrue((jobroot / "README_RUN_LOCAL.md").exists())
        self.assertTrue((jobroot / "AutoDockTools_py3").is_dir())
        self.assertTrue((jobroot / "Receptors").is_dir())
        self.assertTrue((jobroot / "Receptors_PDB").is_dir())
        self.assertFalse((jobroot / "3B_ServerDocks.py").exists())
        self.assertFalse((jobroot / "runDOCKING-tmux.sh").exists())
        self.assertFalse((jobroot / self.RETIRED_SERVER).exists())
        self.assertFalse((jobroot / self.RETIRED_PARSE).exists())
        self.assertIsInstance(warnings, list)

    def test_lsf_package_includes_hpc_files(self):
        jobroot, warnings = assemble_job_tree(self.ws, self.ws / "Receptors", self.ws / "Ligands", package_mode="joey_lsf")

        self.assertTrue((jobroot / "3B_ServerDocks.py").exists())
        self.assertTrue((jobroot / "1B_confgen_batch.py").exists())
        self.assertTrue((jobroot / "3a_PDB2PDBQTbatch.py").exists())
        self.assertTrue((jobroot / "4B_LSFbatch.py").exists())
        self.assertTrue((jobroot / "hpc_profiles.py").exists())
        self.assertTrue((jobroot / "4_ParseScores.py").exists())
        self.assertTrue((jobroot / "4C_ConcatenateScores.py").exists())
        self.assertTrue((jobroot / "5_CompactedHTMLViz.py").exists())
        self.assertTrue((jobroot / "5_COMPACTED_SDF_HTML.py").exists())
        self.assertTrue((jobroot / "lsf_templates.py").exists())
        self.assertTrue((jobroot / "runDOCKING-tmux.sh").exists())
        self.assertIsInstance(warnings, list)

    def test_custom_lsf_package_includes_hpc_files(self):
        jobroot, warnings = assemble_job_tree(self.ws, self.ws / "Receptors", self.ws / "Ligands", package_mode="custom_lsf")

        self.assertTrue((jobroot / "3B_ServerDocks.py").exists())
        self.assertTrue((jobroot / "1B_confgen_batch.py").exists())
        self.assertTrue((jobroot / "4B_LSFbatch.py").exists())
        self.assertTrue((jobroot / "hpc_profiles.py").exists())
        self.assertTrue((jobroot / "lsf_templates.py").exists())
        self.assertIsInstance(warnings, list)

    def test_mainak_lsf_package_includes_hpc_files(self):
        jobroot, warnings = assemble_job_tree(self.ws, self.ws / "Receptors", self.ws / "Ligands", package_mode="mainak_lsf")

        self.assertTrue((jobroot / "3B_ServerDocks.py").exists())
        self.assertTrue((jobroot / "1B_confgen_batch.py").exists())
        self.assertTrue((jobroot / "4B_LSFbatch.py").exists())
        self.assertTrue((jobroot / "hpc_profiles.py").exists())
        self.assertTrue((jobroot / "lsf_templates.py").exists())
        self.assertIsInstance(warnings, list)

    def test_nested_ligands_workspace_is_flattened_in_package(self):
        nested = self.ws / "Ligands" / "Ligands"
        nested.mkdir(parents=True, exist_ok=True)
        (self.ws / "Ligands" / "ligand.sdf").unlink()
        (nested / "DR7.sdf").write_text("dr7\n")

        jobroot, warnings = assemble_job_tree(self.ws, self.ws / "Receptors", self.ws / "Ligands", package_mode="portable")

        self.assertTrue((jobroot / "Ligands" / "DR7.sdf").exists())
        self.assertFalse((jobroot / "Ligands" / "Ligands").exists())
        self.assertTrue(any("Flattened nested ligand path" in warning for warning in warnings))


if __name__ == "__main__":
    unittest.main()
