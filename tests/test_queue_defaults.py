import importlib.util
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, REPO_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class QueueDefaultTests(unittest.TestCase):
    def test_confgen_batch_default_queue_is_general(self):
        module = load_script_module("1B_confgen_batch.py", "queue_default_confgen_batch")
        old_argv = sys.argv[:]
        try:
            sys.argv = ["1B_confgen_batch.py"]
            args = module.parse_args()
        finally:
            sys.argv = old_argv
        self.assertEqual(args.queue, "general")

    def test_parse_batch_default_prompt_queue_is_general(self):
        text = (REPO_ROOT / "4B_LSFbatch.py").read_text(encoding="utf-8")
        self.assertIn('input_default("Queue", "general")', text)

    def test_shared_lsf_template_default_queue_is_general(self):
        module = load_script_module("lsf_templates.py", "queue_default_lsf_templates")
        self.assertEqual(module.DEFAULT_QUEUE, "general")


if __name__ == "__main__":
    unittest.main()
