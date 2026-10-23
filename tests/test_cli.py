import importlib
import subprocess
import sys
from unittest import TestCase


class TestCommandLineSmoke(TestCase):
    """Check that the command line starts on a fresh install."""

    def test_modules_imported_by_main_load(self):
        for module in ['moonstone.main', 'moonstone.analysis.classify', 'moonstone.analysis.randomForest']:
            with self.subTest(module=module):
                importlib.import_module(module)

    def test_help(self):
        result = subprocess.run(
            [sys.executable, '-c', 'from moonstone.main import run; run()', '--help'],
            capture_output=True, text=True, timeout=120,
        )
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('usage: moonstone', result.stdout)
        self.assertIn('--random_forest', result.stdout)
        # parse_arguments() turns every argparse exit, including --help, into exit code 1.
        self.assertEqual(result.returncode, 1)
