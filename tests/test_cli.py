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
        for module in ['moonstone.main', 'moonstone.analysis.classification',
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
        self.assertIn('--classify', result.stdout)
        self.assertNotIn('--random_forest', result.stdout)
        self.assertEqual(result.returncode, 0)

    def test_missing_arguments_exit_1(self):
        result = self.run_moonstone()
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('the following arguments are required', result.stderr)
        self.assertEqual(result.returncode, 1)

    def write_inputs(self):
        tmpdir = tempfile.mkdtemp()
        rng = np.random.default_rng(0)
        samples = [f's{i}' for i in range(40)]
        counts = pd.DataFrame(rng.poisson(20, (12, 40)), index=[f'taxon{i}' for i in range(12)], columns=samples)
        counts.index.name = 'OTU'
        counts.to_csv(os.path.join(tmpdir, 'counts.csv'))
        metadata = pd.DataFrame({'group': ['A'] * 20 + ['B'] * 20}, index=pd.Index(samples, name='sample'))
        metadata.to_csv(os.path.join(tmpdir, 'metadata.csv'))
        return os.path.join(tmpdir, 'counts.csv'), os.path.join(tmpdir, 'metadata.csv'), os.path.join(tmpdir, 'output')

    def test_classify_both_models(self):
        counts, metadata, outdir = self.write_inputs()
        # Feature selection with the random forest is slow, so it is covered in the module tests.
        result = self.run_moonstone(counts, metadata, outdir, '--classify', 'group', '--model', 'svm',
                                    'random_forest', '--roc')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Balanced accuracy on test folds', result.stdout)
        for model in ['svm', 'random_forest']:
            for suffix in ['cv_scores.csv', 'permutation_importance.csv', 'permutation_importance.pdf', 'roc.pdf']:
                with self.subTest(model=model, suffix=suffix):
                    self.assertTrue(os.path.isfile(os.path.join(outdir, f'group_{model}_{suffix}')))

    def test_classify_svm_feature_selection(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--classify', 'group', '--feature_selection')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Features kept per fold', result.stdout)
        for suffix in ['selected_features.csv', 'rfecv.pdf']:
            with self.subTest(suffix=suffix):
                self.assertTrue(os.path.isfile(os.path.join(outdir, f'group_svm_{suffix}')))

    def test_explore_kernels_then_auto_kernel(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--classify', 'group', '--explore_kernels',
                                    '--kernel', 'auto')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Exploratory SVM kernel comparison', result.stdout)
        self.assertIn('Kernel chosen per fold', result.stdout)
        for output_file in ['group_svm_kernel_comparison.csv', 'group_svm_cv_scores.csv']:
            with self.subTest(output_file=output_file):
                self.assertTrue(os.path.isfile(os.path.join(outdir, output_file)))

    def test_classify_error_exits_1(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--classify', 'group', '--kernel', 'rbf',
                                    '--feature_selection')
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('linear kernel', result.stderr)
        self.assertEqual(result.returncode, 1)

    def test_removed_flags_name_replacement(self):
        counts, metadata, outdir = self.write_inputs()
        for flag, replacement in [('--svm', '--classify group'), ('-sr', '--classify group --roc'),
                                  ('--svm_classifier', '--classify group --feature_selection'),
                                  ('-rf', '--classify group --model random_forest')]:
            with self.subTest(flag=flag):
                result = self.run_moonstone(counts, metadata, outdir, flag, 'group')
                self.assertNotIn('Traceback', result.stderr)
                self.assertIn(f'{flag} has been removed. Use: {replacement}', result.stderr)
                self.assertEqual(result.returncode, 1)

    def test_removed_pca_and_kmeans_flags(self):
        counts, metadata, outdir = self.write_inputs()
        for args in [('-p',), ('--pca_plot',), ('-k', '3'), ('--k_means', '3')]:
            with self.subTest(flag=args[0]):
                result = self.run_moonstone(counts, metadata, outdir, *args)
                self.assertNotIn('Traceback', result.stderr)
                self.assertIn(f'{args[0]} has been removed.', result.stderr)
                self.assertEqual(result.returncode, 1)
                self.assertFalse(os.path.exists(os.path.join(outdir, 'metaData_withKClusters.csv')))

    def test_help_does_not_list_removed_flags(self):
        result = self.run_moonstone('--help')
        self.assertNotIn('--pca_plot', result.stdout)
        self.assertNotIn('--k_means', result.stdout)


class TestCommandLinePermanova(TestCase):

    def run_moonstone(self, *args):
        return subprocess.run(
            [sys.executable, '-c', 'from moonstone.main import run; run()', *args],
            capture_output=True, text=True, timeout=300,
        )

    def write_inputs(self):
        from tests.analysis.test_permanova import make_taxonomy_counts
        tmpdir = tempfile.mkdtemp()
        counts_df, metadata = make_taxonomy_counts()
        counts_df.reset_index().to_csv(os.path.join(tmpdir, 'counts.tsv'), sep='\t', index=False)
        metadata.index.name = 'sample'
        metadata.to_csv(os.path.join(tmpdir, 'metadata.csv'))
        return os.path.join(tmpdir, 'counts.tsv'), os.path.join(tmpdir, 'metadata.csv'), os.path.join(tmpdir, 'out')

    def test_permanova_default_species(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', '--permutations', '199')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('species / group', result.stdout)
        table = pd.read_csv(os.path.join(outdir, 'permanova.csv'))
        self.assertEqual(list(table['variable']), ['group', 'random'])
        self.assertEqual(set(table['level']), {'species'})
        self.assertIn('p_value_se', table.columns)

    def test_permanova_threads_same_results(self):
        counts, metadata, outdir = self.write_inputs()
        tables = []
        for threads in ['1', '2']:
            out = outdir + threads
            result = self.run_moonstone(counts, metadata, out, '--permanova', '--permutations', '199',
                                        '--threads', threads)
            self.assertEqual(result.returncode, 0)
            tables.append(pd.read_csv(os.path.join(out, 'permanova.csv')))
        pd.testing.assert_frame_equal(tables[0], tables[1])

    def test_permanova_rejects_negative_threads(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', '--permutations', '9',
                                    '--threads', '-3')
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('Traceback', result.stderr)

    def test_permanova_all_levels_one_variable(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', 'group', '--taxa_level', 'all',
                                    '--permutations', '99')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('Warning (genus)', result.stdout)
        table = pd.read_csv(os.path.join(outdir, 'permanova.csv'))
        self.assertEqual(table.shape[0], 7)
        self.assertTrue(table.loc[table['level'] != 'species', 'warnings'].str.contains('species level only').all())
        self.assertEqual(table.loc[table['level'] == 'kingdom', 'skipped'].iloc[0],
                         'fewer than 2 taxa, composition cannot be compared')

    def test_permanova_no_read_threshold(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', 'group', '--no_read_threshold',
                                    '--permutations', '9')
        self.assertEqual(result.returncode, 0)
        table = pd.read_csv(os.path.join(outdir, 'permanova.csv'))
        self.assertEqual(table.loc[0, 'taxa_kept'], table.loc[0, 'taxa_total'])

    def test_permanova_blank_metadata_is_missing(self):
        counts, metadata, outdir = self.write_inputs()
        lines = open(metadata).read().splitlines()
        lines[1] = lines[1].replace(',A,', ',,', 1)     # s0: blank group
        lines[2] = lines[2].replace(',A,', ',NA,', 1)   # s1: NA group
        open(metadata, 'w').write('\n'.join(lines) + '\n')
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', 'group', '--permutations', '9')
        self.assertEqual(result.returncode, 0)
        self.assertIn('group: 2 sample(s) with a missing value left out', result.stdout)
        table = pd.read_csv(os.path.join(outdir, 'permanova.csv'))
        self.assertEqual(table.loc[0, 'n_missing'], 2)
        self.assertEqual(table.loc[0, 'n_samples'], 28)
        self.assertEqual(table.loc[0, 'n_groups'], 2)

    def test_permanova_mostly_missing_variable_flagged(self):
        counts, metadata, outdir = self.write_inputs()
        lines = open(metadata).read().splitlines()
        lines = [line if i == 0 or i > 16 else line.rsplit(',', 1)[0] + ',' for i, line in enumerate(lines)]
        open(metadata, 'w').write('\n'.join(lines) + '\n')     # 'random' blank for 16 of 30 samples
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', 'random', '--permutations', '9')
        self.assertEqual(result.returncode, 0)
        self.assertIn('random: 16 sample(s) with a missing value left out - mostly missing', result.stdout)

    def test_permanova_rejects_other_analyses(self):
        counts, metadata, outdir = self.write_inputs()
        result = self.run_moonstone(counts, metadata, outdir, '--permanova', '-f', '5')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 1)
        self.assertIn('cannot be combined with -f/--filtering', result.stderr)
