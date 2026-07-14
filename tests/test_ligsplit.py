import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


def load_ligsplit_module():
    repo_root = Path(__file__).resolve().parents[1]
    module_path = repo_root / "0_LIGSPLIT.py"
    spec = importlib.util.spec_from_file_location("ligsplit_module", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ligsplit = load_ligsplit_module()


class LigSplitTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_split_single_csv_into_requested_group_size(self):
        csv_path = self.root / "ligands.csv"
        csv_path.write_text(
            "smiles,id\n"
            "CCO,lig1\n"
            "CCC,lig2\n"
            "CCN,lig3\n"
            "CCCl,lig4\n"
            "CCBr,lig5\n",
            encoding="utf-8",
        )

        result = ligsplit.split_ligands(
            input_path=csv_path,
            output_root=self.root / "csv_split",
            group_size=2,
            prefix="batch",
        )

        self.assertEqual(result["mode"], "csv")
        self.assertEqual(result["group_sizes"], [2, 2, 1])
        self.assertTrue((self.root / "csv_split" / "batch_001" / "ligands_part_001.csv").exists())
        self.assertTrue((self.root / "csv_split" / "batch_003" / "ligands_part_003.csv").exists())

    def test_split_multirecord_sdf_creates_individual_sdf_files(self):
        sdf_path = self.root / "library.sdf"
        sdf_path.write_text(
            "LigA\n  CDK     07142026\n\n  0  0  0  0  0  0            999 V2000\nM  END\n$$$$\n"
            "LigB\n  CDK     07142026\n\n  0  0  0  0  0  0            999 V2000\nM  END\n$$$$\n"
            "LigC\n  CDK     07142026\n\n  0  0  0  0  0  0            999 V2000\nM  END\n$$$$\n",
            encoding="utf-8",
        )

        result = ligsplit.split_ligands(
            input_path=sdf_path,
            output_root=self.root / "sdf_split",
            num_groups=2,
            prefix="group",
        )

        self.assertEqual(result["mode"], "sdf")
        self.assertEqual(result["group_sizes"], [2, 1])
        self.assertTrue((self.root / "sdf_split" / "group_001").exists())
        written = sorted((self.root / "sdf_split").rglob("*.sdf"))
        self.assertEqual(len(written), 3)

    def test_split_smiles_file_writes_one_record_per_file(self):
        smiles_path = self.root / "library.smi"
        smiles_path.write_text(
            "CCO ethanol\n"
            "CCC propane\n"
            "CCN ethylamine\n",
            encoding="utf-8",
        )

        result = ligsplit.split_ligands(
            input_path=smiles_path,
            output_root=self.root / "smiles_split",
            group_size=2,
            prefix="chunk",
        )

        self.assertEqual(result["mode"], "smiles")
        self.assertEqual(result["group_sizes"], [2, 1])
        written = sorted((self.root / "smiles_split").rglob("*.smi"))
        self.assertEqual(len(written), 3)

    def test_mixed_csv_and_sdf_is_rejected(self):
        lig_dir = self.root / "Ligands"
        lig_dir.mkdir()
        (lig_dir / "ligands.csv").write_text("smiles,id\nCCO,lig1\n", encoding="utf-8")
        (lig_dir / "ligand.sdf").write_text(
            "LigA\n  CDK     07142026\n\n  0  0  0  0  0  0            999 V2000\nM  END\n$$$$\n",
            encoding="utf-8",
        )

        with self.assertRaises(ValueError):
            ligsplit.split_ligands(
                input_path=lig_dir,
                output_root=self.root / "bad_split",
                group_size=1,
            )


if __name__ == "__main__":
    unittest.main()
