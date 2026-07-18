import tempfile
import unittest
import urllib.error
from pathlib import Path
from unittest.mock import patch

from packager import fetch_pdb_and_prep


class _FakeResponse:
    def __init__(self, payload: bytes):
        self.payload = payload

    def read(self):
        return self.payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


class ReceptorFetchTests(unittest.TestCase):
    def test_fetch_falls_back_to_cif_when_pdb_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)

            def fake_urlopen(url):
                if url.endswith(".pdb"):
                    raise urllib.error.HTTPError(url, 404, "Not Found", hdrs=None, fp=None)
                if url.endswith(".cif"):
                    return _FakeResponse(b"data_9BA9\n#\n")
                raise AssertionError(url)

            with patch("packager.urllib.request.urlopen", side_effect=fake_urlopen):
                out = fetch_pdb_and_prep("9BA9", dest)

            self.assertTrue(out["pdb_path"].endswith("9ba9.cif"))
            self.assertEqual((dest / "9ba9.cif").read_text(), "data_9BA9\n#\n")

    def test_fetch_rejects_chain_filter_for_cif_only_entries(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp)

            def fake_urlopen(url):
                if url.endswith(".pdb"):
                    raise urllib.error.HTTPError(url, 404, "Not Found", hdrs=None, fp=None)
                if url.endswith(".cif"):
                    return _FakeResponse(b"data_9BA9\n#\n")
                raise AssertionError(url)

            with patch("packager.urllib.request.urlopen", side_effect=fake_urlopen):
                with self.assertRaises(ValueError) as ctx:
                    fetch_pdb_and_prep("9BA9", dest, chains="A")

            self.assertIn("only available as mmCIF", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
