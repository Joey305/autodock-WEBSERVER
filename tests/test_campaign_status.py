import json
import tempfile
import unittest
from pathlib import Path

from campaign_status import campaign_progress


class CampaignStatusTests(unittest.TestCase):
    def test_uses_actual_vina_jobs_and_deduplicates_pose_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "campaign.json").write_text(json.dumps({
                "campaign": "example",
                "batches": [
                    {"name": "batch_001", "input_ligands": 2, "pdbqt_maximum": 32, "vina_maximum": 64},
                    {"name": "batch_002", "input_ligands": 2, "pdbqt_maximum": 32, "vina_maximum": 64},
                ],
            }))
            pdbqt = root / "batches" / "batch_001" / "input_Ligands_PDBQT_16Poses"
            pdbqt.mkdir(parents=True)
            (pdbqt / "one.pdbqt").write_text("REMARK\n")
            results = root / "results" / "batch_001"
            results.mkdir(parents=True)
            (results / "docking_jobs.csv").write_text("ResolvedReceptor,Ligand,Seed\nr,one,1\nr,two,2\nr,three,3\n")
            (results / "compact_results.csv").write_text("ResolvedReceptor,Ligand,Seed,Status\nr,one,1,success\nr,one,1,success\n")
            (results / "compact_failures.csv").write_text("ResolvedReceptor,Ligand,Seed,Status\nr,two,2,failed\n")

            report = campaign_progress(root)

        first, second = report["batches"]
        self.assertEqual(first["pdbqt"]["produced"], 1)
        self.assertEqual(first["pdbqt"]["maximum"], 32)
        self.assertTrue(first["vina"]["total_is_exact"])
        self.assertEqual(first["vina"]["total"], 3)
        self.assertEqual(first["vina"]["completed"], 2)
        self.assertEqual(first["vina"]["successful"], 1)
        self.assertEqual(first["vina"]["failed"], 1)
        self.assertFalse(second["vina"]["total_is_exact"])
        self.assertEqual(second["vina"]["total"], 64)
        self.assertFalse(report["vina"]["total_is_exact"])
        self.assertEqual(report["vina"]["total"], 67)


if __name__ == "__main__":
    unittest.main()
