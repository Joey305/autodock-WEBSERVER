import tempfile
import unittest
from pathlib import Path


class ParserWorkflowTests(unittest.TestCase):
    def test_lsf_batch_references_new_parse_script(self):
        text = Path("4B_LSFbatch.py").read_text()
        self.assertIn("python 4_ParseScores.py", text)
        retired_call = "python " + "4C" + "MULTIVINA.py"
        self.assertNotIn(retired_call, text)

    def test_runner_readme_mentions_new_parse_flow(self):
        text = Path("runner_templates.py").read_text()
        self.assertIn("4_ParseScores.py", text)
        self.assertIn("4C_ConcatenateScores.py", text)
        self.assertIn("ligand_state_manifest.csv", text)

    def test_retired_parser_files_are_gone(self):
        for retired in [
            "4A" + "_Parse_VinaResults.py",
            "4A" + "_PARSEvinascores.py",
            "4B" + "_Parse_VinaResults.py",
            "4C" + "MULTIVINA.py",
            "4F" + "_Combine_PerDirCSVs_with_Provenance.py",
            "4" + "_SERVERBATCH.py",
        ]:
            self.assertFalse(Path(retired).exists(), retired)


if __name__ == "__main__":
    unittest.main()
