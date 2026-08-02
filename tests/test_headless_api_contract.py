import importlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


class HeadlessApiContractTests(unittest.TestCase):
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

    def test_health(self):
        response = self.client.get("/api/v1/health")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["data"]["api_version"], "v1")

    def test_workspace_create(self):
        response = self.client.post("/api/v1/workspaces", json={"workspace_name": "api-test"})
        self.assertEqual(response.status_code, 201)
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["data"]["jobname"], "api-test")

    def test_center_resolve_explicit_xyz(self):
        workspace = self.client.post("/api/v1/workspaces", json={"workspace_name": "xyz-test"}).get_json()["data"]
        response = self.client.post(
            f"/api/v1/workspaces/{workspace['jobname']}/centers/resolve",
            json={"method": "xyz", "center": [1, 2, 3], "size": 20},
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["data"]["center"], [1.0, 2.0, 3.0])

    def test_center_resolve_error_response(self):
        workspace = self.client.post("/api/v1/workspaces", json={"workspace_name": "error-test"}).get_json()["data"]
        response = self.client.post(
            f"/api/v1/workspaces/{workspace['jobname']}/centers/resolve",
            json={"method": "xyz", "center": [1, 2]},
        )
        self.assertEqual(response.status_code, 400)
        payload = response.get_json()
        self.assertFalse(payload["ok"])
        self.assertEqual(payload["error"], "invalid_xyz")

    def test_hetatm_listing_endpoint(self):
        workspace = self.client.post("/api/v1/workspaces", json={"workspace_name": "hetatm-list"}).get_json()["data"]
        receptor = "\n".join(
            [
                "HETATM    1  C1  DR7 A 100      10.000  20.000  30.000  1.00 20.00           C",
                "HETATM    2  C2  DR7 A 100      12.000  22.000  32.000  1.00 20.00           C",
                "HETATM    3  O   HOH A 200       1.000   1.000   1.000  1.00 20.00           O",
                "END",
                "",
            ]
        ).encode("utf-8")
        upload = self.client.post(
            f"/api/v1/workspaces/{workspace['jobname']}/receptors/upload",
            data={"mode": "single", "file": (io.BytesIO(receptor), "fixture.pdb")},
            content_type="multipart/form-data",
        )
        self.assertEqual(upload.status_code, 200, upload.get_data(as_text=True))

        response = self.client.get(f"/api/v1/workspaces/{workspace['jobname']}/hetatms")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["data"]["count"], 1)
        self.assertEqual(payload["data"]["instances"][0]["resname"], "DR7")
        self.assertEqual(payload["data"]["instances"][0]["center"], [11.0, 21.0, 31.0])

    def test_headless_package_from_bound_ligand(self):
        workspace = self.client.post("/api/v1/workspaces", json={"workspace_name": "one-shot"}).get_json()["data"]
        receptor = "\n".join(
            [
                "ATOM      1  CA  GLU A  82       1.000   2.000   3.000  1.00 20.00           C",
                "HETATM    2  C1  DR7 A 100      10.000  20.000  30.000  1.00 20.00           C",
                "HETATM    3  C2  DR7 A 100      12.000  22.000  32.000  1.00 20.00           C",
                "END",
                "",
            ]
        ).encode("utf-8")
        upload = self.client.post(
            f"/api/v1/workspaces/{workspace['jobname']}/receptors/upload",
            data={"mode": "single", "file": (io.BytesIO(receptor), "fixture.pdb")},
            content_type="multipart/form-data",
        )
        self.assertEqual(upload.status_code, 200, upload.get_data(as_text=True))

        def fake_run(cmd, stdout=None, stderr=None, text=None):
            if "-osdf" in cmd:
                Path(cmd[-1]).write_text("DR7\n  OpenAI\n\nM  END\n$$$$\n", encoding="utf-8")
            else:
                Path(cmd[3]).write_text("RECEPTOR PDBQT\n", encoding="utf-8")
            return self.app_module.subprocess.CompletedProcess(cmd, 0, "", "")

        with patch("app.subprocess.run", side_effect=fake_run):
            response = self.client.post(
                "/api/v1/headless/package",
                json={
                    "workspace_name": workspace["jobname"],
                    "reuse": True,
                    "receptor": {"rel": "fixture.pdb"},
                    "bound_ligand": {"resname": "DR7", "chain": "A", "resi": "100"},
                    "center": {"method": "same_as_bound_ligand", "size": 20},
                    "prep": {"remove_hets": "all"},
                    "package": {"package_mode": "portable", "poses_conf": 8, "poses_vina": 3},
                },
            )

        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        payload = response.get_json()
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["data"]["center"]["center"], [11.0, 21.0, 31.0])
        self.assertTrue(Path(payload["data"]["artifact"]["zip"]).exists())


if __name__ == "__main__":
    unittest.main()
