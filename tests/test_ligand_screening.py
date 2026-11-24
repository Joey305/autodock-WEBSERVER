import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from rdkit import Chem


ROOT = Path(__file__).resolve().parent.parent


def load_confgen():
    spec = importlib.util.spec_from_file_location("confgen_screening", ROOT / "1_ConformerGeneration.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_docking():
    spec = importlib.util.spec_from_file_location("docking_screening", ROOT / "3_Complete_batch_docking.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ParentLigandScreeningTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_confgen()
        cls.policy = {
            "enabled": True,
            "max_heavy_atoms": 150,
            "max_rotatable_bonds": 32,
            "max_molecular_weight": 1500.0,
            "allow_large_ligands": False,
        }

    def test_drug_like_and_protac_scale_ligands_pass(self):
        drug = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        # This flexible but ordinary 28-carbon scaffold stays below the policy's 32-rotor limit.
        protac_scale = Chem.MolFromSmiles("C" * 28)
        self.assertTrue(self.module.screen_parent_molecule(drug, self.policy)["eligible"])
        self.assertTrue(self.module.screen_parent_molecule(protac_scale, self.policy)["eligible"])

    def test_each_limit_and_multiple_limits_are_reported(self):
        heavy = Chem.MolFromSmiles("C" * 151)
        result = self.module.screen_parent_molecule(
            heavy, {**self.policy, "max_rotatable_bonds": 1000, "max_molecular_weight": 1_000_000}
        )
        self.assertIn("heavy_atom_count", result["rejection_reasons"][0])

        rotatable = Chem.MolFromSmiles("C" * 40)
        result = self.module.screen_parent_molecule(
            rotatable, {**self.policy, "max_heavy_atoms": 1000, "max_molecular_weight": 1_000_000}
        )
        self.assertTrue(any("rotatable_bond_count" in reason for reason in result["rejection_reasons"]))

        high_mw = Chem.MolFromSmiles(".".join(["[I-]"] * 12))
        result = self.module.screen_parent_molecule(
            high_mw, {**self.policy, "max_heavy_atoms": 1000, "max_rotatable_bonds": 1000}
        )
        self.assertTrue(any("molecular_weight" in reason for reason in result["rejection_reasons"]))

        peptide_like = Chem.MolFromSmiles("NCC(=O)" * 100)
        result = self.module.screen_parent_molecule(peptide_like, self.policy)
        self.assertFalse(result["eligible"])
        self.assertEqual(len(result["rejection_reasons"]), 3)
        self.assertGreater(result["heavy_atom_count"], 360)
        self.assertGreater(result["rotatable_bond_count"], 150)

    def test_csv_and_sdf_tasks_are_filtered_before_workers(self):
        peptide = Chem.AddHs(Chem.MolFromSmiles("NCC(=O)" * 100))
        task = {
            "ligand_id": "large", "molblock": Chem.MolToMolBlock(peptide),
            "source_input": "input.csv", "source_input_type": "csv",
        }
        eligible, excluded = self.module.screen_parent_tasks([task], self.policy)
        self.assertEqual(eligible, [])
        self.assertEqual(excluded[0]["ligand_id"], "large")

        with tempfile.TemporaryDirectory() as tmp:
            sdf = Path(tmp) / "input.sdf"
            writer = Chem.SDWriter(str(sdf))
            mol = Chem.MolFromSmiles("NCC(=O)" * 100)
            mol.SetProp("_Name", "peptide_record")
            writer.write(mol)
            writer.close()
            tasks = list(self.module.iter_mols_from_sdf_file(sdf))
        eligible, excluded = self.module.screen_parent_tasks(tasks, self.policy)
        self.assertEqual(eligible, [])
        self.assertEqual(excluded[0]["source_input_type"], "sdf")

    def test_exclusion_report_and_override_are_intentional(self):
        peptide = Chem.AddHs(Chem.MolFromSmiles("NCC(=O)" * 100))
        task = {
            "ligand_id": "large", "molblock": Chem.MolToMolBlock(peptide),
            "source_input": "input.csv", "source_input_type": "csv", "source_record_index": 1,
        }
        eligible, excluded = self.module.screen_parent_tasks([task], self.policy)
        self.assertFalse(eligible)
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            report = self.module.write_excluded_ligands_report(output, excluded, self.policy)
            self.assertEqual(list(output.glob("*.pdbqt")), [])
            with report.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["LigandID"], "large")
            self.assertIn("HeavyAtomCount", rows[0])
            self.assertIn("RejectionReasons", rows[0])

        override = {**self.policy, "enabled": False, "allow_large_ligands": True}
        eligible, excluded = self.module.screen_parent_tasks([task], override)
        self.assertEqual(len(eligible), 1)
        self.assertEqual(excluded, [])

    def test_excluded_csv_parent_never_enters_generation_or_writes_pdbqt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "ligands.csv"
            source.write_text("smiles,id\n" + "NCC(=O)" * 100 + ",peptide\n", encoding="utf-8")
            original_cwd = Path.cwd()
            try:
                import os
                os.chdir(root)
                self.module.run_csv_mode(
                    source,
                    SimpleNamespace(smiles_col="smiles", id_col="id", remove_tmp=False),
                    1, 1, "unused-obabel", {}, self.policy,
                )
            finally:
                os.chdir(original_cwd)
            outputs = list(root.glob("ligands_Ligands_PDBQT_*"))
            self.assertEqual(len(outputs), 1)
            self.assertEqual(list(outputs[0].glob("*.pdbqt")), [])
            with (outputs[0] / "excluded_ligands.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 1)


class DockingStageGuardTests(unittest.TestCase):
    def test_legacy_pdbqt_guard_rejects_extreme_input_and_can_be_overridden(self):
        module = load_docking()
        with tempfile.TemporaryDirectory() as tmp:
            ligand = Path(tmp) / "legacy_peptide.pdbqt"
            ligand.write_text("".join([f"ATOM  {i:5d}  C   LIG     1\n" for i in range(1, 361)]) + "TORSDOF 150\n", encoding="utf-8")
            eligible, excluded = module.filter_pdbqt_ligands([str(ligand)], max_atoms=150, max_torsions=32)
            self.assertEqual(eligible, [])
            self.assertEqual(excluded[0]["atom_count"], 360)
            self.assertEqual(excluded[0]["torsdof"], 150)
            self.assertEqual(len(excluded[0]["reasons"]), 2)
            eligible, excluded = module.filter_pdbqt_ligands(
                [str(ligand)], max_atoms=150, max_torsions=32, allow_large_ligands=True
            )
            self.assertEqual(eligible, [str(ligand)])
            self.assertEqual(excluded, [])


if __name__ == "__main__":
    unittest.main()
