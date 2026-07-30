import importlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


class BoundLigandExtractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import app as app_module
        except ModuleNotFoundError as exc:
            raise unittest.SkipTest(f"Flask app dependencies are not installed: {exc}") from exc

        cls._tmp = tempfile.TemporaryDirectory()
        tmp_root = Path(cls._tmp.name)
        db_path = tmp_root / "test.db"
        cls.workspace_root = tmp_root / "workspaces"

        app_module.Config.SQLALCHEMY_DATABASE_URI = f"sqlite:///{db_path}"
        app_module.Config.TMP_ROOT = str(cls.workspace_root)
        app_module.Config.PUBLIC_EMAIL = "public@autodock.local"
        cls.app_module = importlib.reload(app_module)
        cls.app_module.Config.SQLALCHEMY_DATABASE_URI = f"sqlite:///{db_path}"
        cls.app_module.Config.TMP_ROOT = str(cls.workspace_root)
        cls.app_module.Config.PUBLIC_EMAIL = "public@autodock.local"

        cls.app = cls.app_module.create_app()
        cls.client = cls.app.test_client()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "_tmp"):
            cls._tmp.cleanup()

    def _create_workspace_with_receptor(self, name="extract-test", extra_dr7: bool = False):
        workspace = self.client.post("/api/v1/workspaces", json={"workspace_name": name}).get_json()["data"]
        ws = self.workspace_root / workspace["jobname"]
        receptor_dir = ws / "Receptors"
        receptor_dir.mkdir(parents=True, exist_ok=True)
        pdb = receptor_dir / "3eky.pdb"
        lines = [
            "HETATM    1  CAA DR7 A 101      11.000  21.000  31.000  1.00 20.00           C  ",
            "HETATM    2  NBG DR7 A 101      12.000  21.000  31.000  1.00 20.00           N  ",
            "HETATM    3  C1  ATP B 201       1.000   2.000   3.000  1.00 20.00           C  ",
        ]
        if extra_dr7:
            lines.append("HETATM    4  C1  DR7 B 102      15.000  25.000  35.000  1.00 20.00           C  ")
        pdb.write_text("\n".join(lines + ["END", ""]), encoding="utf-8")
        state = {
            "receptors": [{"rel": "Receptors/3eky.pdb", "display": "3eky.pdb", "status": "new"}],
            "centers_csv": "vina_centers.csv",
            "ligands_uploaded": False,
            "ligand_info": {},
            "prep_job": None,
        }
        (ws / "_state.json").write_text(self.app_module.json.dumps(state), encoding="utf-8")
        return workspace, ws

    def test_extract_bound_hetatm_writes_sdf_and_marks_ligands_uploaded(self):
        workspace, ws = self._create_workspace_with_receptor("extract-dr7")

        def fake_run(cmd, stdout=None, stderr=None, text=None):
            pdb_lines = Path(cmd[2]).read_text(encoding="utf-8").splitlines()
            self.assertEqual(pdb_lines[0][12:16], " CAA")
            self.assertEqual(pdb_lines[0][17:20], "DR7")
            self.assertEqual(pdb_lines[0][21], "A")
            self.assertEqual(pdb_lines[0][22:26].strip(), "101")
            self.assertEqual(pdb_lines[0][76:78].strip(), "C")
            self.assertEqual(pdb_lines[1][12:16], " NBG")
            self.assertEqual(pdb_lines[1][76:78].strip(), "N")
            out_path = Path(cmd[-1])
            out_path.write_text("DR7\n  OpenAI\n\nM  END\n$$$$\n", encoding="utf-8")
            return self.app_module.subprocess.CompletedProcess(cmd, 0, "", "")

        with patch("app.subprocess.run", side_effect=fake_run):
            response = self.client.post(
                "/api/ligands/extract",
                json={
                    "jobname": workspace["jobname"],
                    "rel": "Receptors/3eky.pdb",
                    "resname": "DR7",
                    "chain": "A",
                    "resi": "101",
                },
            )

        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["extracted"]["target_name"], "DR7.sdf")
        self.assertEqual(payload["extracted"]["target_rel"], "Ligands/DR7.sdf")
        self.assertIsInstance(payload["extracted"]["target_path"], str)
        self.assertTrue((ws / "Ligands" / "DR7.sdf").exists())

        state = self.app_module.json.loads((ws / "_state.json").read_text(encoding="utf-8"))
        self.assertTrue(state["ligands_uploaded"])
        self.assertEqual(state["ligand_info"]["upload_mode"], "extracted")
        self.assertEqual(state["ligand_info"]["accepted_files"], ["DR7.sdf"])
        self.assertEqual(state["ligand_info"]["source"]["resname"], "DR7")

        summary = self.client.get(f"/api/summary?jobname={workspace['jobname']}").get_json()
        self.assertTrue(summary["ligands_uploaded"])
        self.assertEqual(summary["ligand_info"]["filename"], "DR7.sdf")

    def test_extract_bound_hetatm_reports_ambiguous_instances(self):
        workspace, _ws = self._create_workspace_with_receptor("extract-ambiguous", extra_dr7=True)
        response = self.client.post(
            "/api/ligands/extract",
            json={
                "jobname": workspace["jobname"],
                "rel": "Receptors/3eky.pdb",
                "resname": "DR7",
            },
        )

        self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
        payload = response.get_json()
        self.assertEqual(payload["error"], "ambiguous_selection")
        self.assertEqual(len(payload["details"]["candidates"]), 2)


if __name__ == "__main__":
    unittest.main()
