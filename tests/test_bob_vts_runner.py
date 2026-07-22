import importlib
import io
import tempfile
import time
import unittest
import zipfile
from pathlib import Path


class BobVtsRunnerTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        tmp_root = Path(self._tmp.name)
        db_path = tmp_root / "test.db"

        import app as app_module

        self.app_module = importlib.reload(app_module)
        self.app_module.Config.SQLALCHEMY_DATABASE_URI = f"sqlite:///{db_path}"
        self.app_module.Config.TMP_ROOT = str(tmp_root / "workspaces")
        self.app_module.Config.VTS_RUNNER_ROOT = str(tmp_root / "vts_runs")
        self.app_module.Config.VTS_RUNNER_CONDA_ENV = ""
        self.app_module.Config.VTS_RUNNER_ENV_LINE = ""
        self.app_module.Config.VTS_RUNNER_TIMEOUT_SECONDS = 30
        self.app_module.Config.VTS_RUNNER_MAX_ZIP_MB = 10
        self.app_module.Config.VTS_RUNNER_TOKEN = ""

        self.app = self.app_module.create_app()
        self.client = self.app.test_client()

    def tearDown(self):
        self._tmp.cleanup()

    def _make_package(self):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(
                "Tiny_VTS_Package/run_tiny_macos_linux.sh",
                "#!/usr/bin/env bash\nset -euo pipefail\necho run-start\nprintf 'done\\n' > result.txt\necho run-end\n",
            )
            zf.writestr("Tiny_VTS_Package/README.md", "tiny package")
        buf.seek(0)
        return buf

    def test_hidden_page_renders(self):
        response = self.client.get("/BOB-VTS.html")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"BOB VTS Runner", response.data)
        self.assertEqual(response.headers["X-Robots-Tag"], "noindex, nofollow")

    def test_hidden_page_lowercase_alias_renders(self):
        response = self.client.get("/bob-vts.html")
        self.assertEqual(response.status_code, 200)
        self.assertIn(b"BOB VTS Runner", response.data)
        self.assertEqual(response.headers["X-Robots-Tag"], "noindex, nofollow")

    def test_zip_upload_run_and_download(self):
        response = self.client.post(
            "/api/bob-vts/runs",
            data={"package": (self._make_package(), "tiny.zip")},
            content_type="multipart/form-data",
        )
        self.assertEqual(response.status_code, 202)
        payload = response.get_json()
        run_id = payload["run_id"]

        status = payload
        for _ in range(30):
            if status["status"] in {"succeeded", "failed"}:
                break
            time.sleep(0.2)
            status = self.client.get(f"/api/bob-vts/runs/{run_id}").get_json()

        self.assertEqual(status["status"], "succeeded")
        self.assertIn("run-end", status["log_tail"])
        self.assertIn("download_url", status)

        download = self.client.get(status["download_url"])
        self.assertEqual(download.status_code, 200)
        downloaded = io.BytesIO(download.data)
        with zipfile.ZipFile(downloaded) as zf:
            self.assertIn("Tiny_VTS_Package/result.txt", zf.namelist())
            self.assertEqual(zf.read("Tiny_VTS_Package/result.txt"), b"done\n")
            self.assertIn("Tiny_VTS_Package/bob_vts_terminal_output.txt", zf.namelist())


if __name__ == "__main__":
    unittest.main()
