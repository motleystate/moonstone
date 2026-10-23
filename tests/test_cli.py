import importlib
import os
import subprocess
import sys
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd


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

    def test_random_forest(self):
        tmpdir = tempfile.mkdtemp()
        rng = np.random.default_rng(0)
        samples = [f's{i}' for i in range(40)]
        counts = pd.DataFrame(rng.poisson(20, (12, 40)), index=[f'taxon{i}' for i in range(12)], columns=samples)
        counts.index.name = 'OTU'
        counts.to_csv(os.path.join(tmpdir, 'counts.csv'))
        metadata = pd.DataFrame({'group': ['A'] * 20 + ['B'] * 20}, index=pd.Index(samples, name='sample'))
        metadata.to_csv(os.path.join(tmpdir, 'metadata.csv'))
        outdir = os.path.join(tmpdir, 'output')

        result = self.run_moonstone(os.path.join(tmpdir, 'counts.csv'), os.path.join(tmpdir, 'metadata.csv'),
                                    outdir, '--random_forest', 'group')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Gradient Boost:', result.stdout)
        for output_file in ['rf_Allfeatures.csv', 'group_rfFeatures.pdf']:
            with self.subTest(output_file=output_file):
                self.assertTrue(os.path.isfile(os.path.join(outdir, output_file)))

    def test_svm_classifier(self):
        tmpdir = tempfile.mkdtemp()
        rng = np.random.default_rng(0)
        samples = [f's{i}' for i in range(40)]
        counts = pd.DataFrame(rng.poisson(20, (12, 40)), index=[f'taxon{i}' for i in range(12)], columns=samples)
        counts.index.name = 'OTU'
        counts.to_csv(os.path.join(tmpdir, 'counts.csv'))
        metadata = pd.DataFrame({'group': ['A'] * 20 + ['B'] * 20}, index=pd.Index(samples, name='sample'))
        metadata.to_csv(os.path.join(tmpdir, 'metadata.csv'))
        outdir = os.path.join(tmpdir, 'output')

        result = self.run_moonstone(os.path.join(tmpdir, 'counts.csv'), os.path.join(tmpdir, 'metadata.csv'),
                                    outdir, '--svm_classifier', 'group')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        for output_file in ['group_variable_coefs.pdf', 'group_SVM_components.csv', 'group_rfecv.pdf']:
            with self.subTest(output_file=output_file):
                self.assertTrue(os.path.isfile(os.path.join(outdir, output_file)))
