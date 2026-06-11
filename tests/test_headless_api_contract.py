import importlib
import tempfile
import unittest
from pathlib import Path


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


if __name__ == "__main__":
    unittest.main()
