import importlib
import subprocess
import sys
from unittest import TestCase


class TestCommandLineSmoke(TestCase):
    """Check that the command line starts on a fresh install."""

    def test_modules_imported_by_main_load(self):
        for module in ['moonstone.main', 'moonstone.analysis.classify', 'moonstone.analysis.randomForest',
                       'moonstone.normalization.reads.read_downsize', 'moonstone.normalization.reads.downsize_dir']:
            with self.subTest(module=module):
                importlib.import_module(module)

    def run_moonstone(self, *args):
        return subprocess.run(
            [sys.executable, '-c', 'from moonstone.main import run; run()', *args],
            capture_output=True, text=True, timeout=120,
        )

    def test_help(self):
        result = self.run_moonstone('--help')
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('usage: moonstone', result.stdout)
        self.assertIn('--random_forest', result.stdout)
        self.assertEqual(result.returncode, 0)

    def test_missing_arguments_exit_1(self):
        result = self.run_moonstone()
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('the following arguments are required', result.stderr)
        self.assertEqual(result.returncode, 1)
