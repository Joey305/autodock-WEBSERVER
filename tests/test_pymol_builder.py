import csv
import importlib
import importlib.util
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


class PymolBuilderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.have_rdkit = importlib.util.find_spec("rdkit") is not None

    def test_infer_tmp_sdf_relative_path(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_tmp_path")
        rel = module.infer_tmp_sdf_relative_path("DR7_p01_t02_c005")
        self.assertEqual(rel, Path("DR7") / "DR7_p01_t02" / "DR7_p01_t02_c005.sdf")

    def test_manifest_lookup_returns_exact_variant_sdf(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_manifest_lookup")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            manifest_dir = root / "Ligands_Ligands_PDBQT_8Poses_20260611_1718"
            manifest_dir.mkdir()
            tmp_root = root / "Ligands_TMP_SDF_8Poses_20260611_1718" / "DR7" / "DR7_p01_t02"
            tmp_root.mkdir(parents=True)
            sdf_path = tmp_root / "DR7_p01_t02_c005.sdf"
            sdf_path.write_text("fake sdf\n")
            with (manifest_dir / "ligand_state_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["LigandVariant", "TmpSDFFile", "SDFFile"])
                writer.writeheader()
                writer.writerow(
                    {
                        "LigandVariant": "DR7_p01_t02_c005",
                        "TmpSDFFile": "DR7/DR7_p01_t02/DR7_p01_t02_c005.sdf",
                        "SDFFile": "DR7_p01_t02_c005.sdf",
                    }
                )

            resolver = module.ReferenceResolver(root)
            resolved, method, warning, manifest_row = resolver.resolve({"LigandVariant": "DR7_p01_t02_c005", "LigandBase": "DR7"})
            self.assertEqual(resolved, sdf_path.resolve())
            self.assertEqual(method, "manifest:TmpSDFFile")
            self.assertEqual(warning, "")
            self.assertEqual(manifest_row["LigandVariant"], "DR7_p01_t02_c005")

    def test_extract_selected_pose_pdbqt_uses_requested_model(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_pose_extract")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            src = root / "out.pdbqt"
            src.write_text(
                "MODEL 1\nREMARK VINA RESULT: -7.0 0 0\nATOM A\nENDMDL\n"
                "MODEL 2\nREMARK VINA RESULT: -8.0 0 0\nATOM B\nENDMDL\n"
            )
            dest = root / "pose2.pdbqt"
            result, warning = module.extract_selected_pose_pdbqt(src, 2, dest)
            self.assertEqual(result, dest)
            text = dest.read_text()
            self.assertIn("MODEL 2", text)
            self.assertIn("ATOM B", text)
            self.assertNotIn("ATOM A", text)
            self.assertEqual(warning, "")

    def test_group_names_do_not_collide_with_receptor_object(self):
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_groups")
        names = module.build_pymol_names("3eky", "DR7", "DR7_p01_t02_c005", 1)
        self.assertEqual(names["receptor_obj"], "obj_3eky")
        self.assertEqual(names["receptor_group"], "grp_3eky")
        self.assertEqual(names["ligand_group"], "grp_3eky_DR7")
        self.assertNotEqual(names["receptor_obj"], names["receptor_group"])

    def _build_test_mol(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CC(O)N"))
        AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
        return mol

    def test_strip_all_hydrogens_preserves_heavy_atoms(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_strip_h")
        mol = self._build_test_mol()
        stripped = module.strip_all_hydrogens(mol)
        counts = module.atom_count_summary(stripped)
        self.assertEqual(counts["hydrogen"], 0)
        self.assertTrue(module.heavy_atoms_preserved(mol, stripped))

    def test_add_nonpolar_hydrogens_only_keeps_only_carbon_bound_h(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_nonpolar_h")
        mol = self._build_test_mol()
        processed = module.add_nonpolar_hydrogens_only(mol)
        self.assertTrue(module.heavy_atoms_preserved(mol, processed))
        for atom in processed.GetAtoms():
            if atom.GetAtomicNum() != 1:
                continue
            neighbor = atom.GetNeighbors()[0]
            self.assertEqual(neighbor.GetAtomicNum(), 6)

    def test_apply_hydrogen_policy_modes(self):
        if not self.have_rdkit:
            self.skipTest("rdkit not installed")
        module = load_script_module("5C_BuildPymolSesh.py", "pymol_builder_h_modes")
        mol = self._build_test_mol()
        all_counts = module.atom_count_summary(mol)

        none_mol, none_msg = module.apply_hydrogen_policy(mol, "none")
        self.assertEqual(module.atom_count_summary(none_mol)["hydrogen"], 0)
        self.assertIn("Removed all hydrogens", none_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, none_mol))

        nonpolar_mol, nonpolar_msg = module.apply_hydrogen_policy(mol, "nonpolar")
        self.assertLess(module.atom_count_summary(nonpolar_mol)["hydrogen"], all_counts["hydrogen"])
        self.assertIn("carbon-bound hydrogens", nonpolar_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, nonpolar_mol))

        all_mol, all_msg = module.apply_hydrogen_policy(mol, "all")
        self.assertEqual(module.atom_count_summary(all_mol)["hydrogen"], all_counts["hydrogen"])
        self.assertIn("Preserved all hydrogens", all_msg)
        self.assertTrue(module.heavy_atoms_preserved(mol, all_mol))


if __name__ == "__main__":
    unittest.main()
