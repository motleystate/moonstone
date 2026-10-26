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

    def test_filtering_writes_filtered_table(self):
        counts, metadata, outdir = self.write_inputs()
        table = pd.read_csv(counts, index_col=0)
        table.loc['taxon0'] = 1   # mean of 1 read per sample, below the threshold
        table.to_csv(counts)
        result = self.run_moonstone(counts, metadata, outdir, '-f', '5')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        filtered = pd.read_csv(os.path.join(outdir, 'filteredCountFile.csv'), index_col=0)
        self.assertEqual(filtered.shape, (40, 11))
        self.assertNotIn('taxon0', filtered.columns)
        self.assertTrue(os.path.isfile(os.path.join(outdir, 'filtered_variables.csv')))

    def test_sample_without_metadata_warns_and_continues(self):
        counts, metadata, outdir = self.write_inputs()
        pd.read_csv(metadata, index_col=0).iloc[:-1].to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group', '--folds', '3')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('1 sample(s) of the count table have no metadata and are left out: s39', result.stdout)
        self.assertIn('39 samples in both the count and metadata files', result.stdout)
        self.assertIn('39 samples: A 20', result.stdout)

    def test_extra_metadata_samples_are_ignored(self):
        # Used to exit with "sample numbers don't match" although every counted sample had metadata.
        counts, metadata, outdir = self.write_inputs()
        extra = pd.DataFrame({'group': ['A'] * 3}, index=pd.Index(['x1', 'x2', 'x3'], name='sample'))
        pd.concat([pd.read_csv(metadata, index_col=0), extra]).to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group', '--folds', '3')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('3 sample(s) of the metadata are not in the count table and are ignored: x1, x2, x3',
                      result.stdout)
        self.assertIn('40 samples in both the count and metadata files', result.stdout)

    def test_same_number_of_samples_with_different_ids_warns(self):
        # Used to print "Nice! 40 samples" and silently train on the 30 matching samples.
        counts, metadata, outdir = self.write_inputs()
        df = pd.read_csv(metadata, index_col=0)
        df.index = [s.upper() if i >= 30 else s for i, s in enumerate(df.index)]
        df.rename_axis('sample').to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group', '--folds', '3')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('10 sample(s) of the count table have no metadata and are left out', result.stdout)
        self.assertIn("differ only by case or spaces (not matched): 's30' / 'S30'", result.stdout)
        self.assertIn('30 samples in both the count and metadata files', result.stdout)

    def test_no_common_sample_exits_1(self):
        counts, metadata, outdir = self.write_inputs()
        df = pd.read_csv(metadata, index_col=0)
        df.index = [f'p{i}' for i in range(len(df))]
        df.rename_axis('sample').to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group')
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('No sample ID is found in both the count table and the metadata', result.stderr)
        self.assertEqual(result.returncode, 1)

    def test_duplicated_metadata_sample_exits_1(self):
        counts, metadata, outdir = self.write_inputs()
        df = pd.read_csv(metadata, index_col=0)
        pd.concat([df, df.iloc[[0]]]).to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group')
        self.assertNotIn('Traceback', result.stderr)
        self.assertIn('Sample IDs used more than once in the metadata: s0', result.stderr)
        self.assertEqual(result.returncode, 1)

    def test_numeric_sample_ids_keep_leading_zeros(self):
        counts, metadata, outdir = self.write_inputs()
        ids = [f'{i:03d}' for i in range(40)]
        counts_df = pd.read_csv(counts, index_col=0)
        counts_df.columns = ids
        counts_df.to_csv(counts)
        df = pd.read_csv(metadata, index_col=0)
        df.index = pd.Index(ids, name='sample')
        df.sample(frac=1, random_state=0).to_csv(metadata)
        result = self.run_moonstone(counts, metadata, outdir, '-c', 'group', '--folds', '3')
        self.assertNotIn('Traceback', result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn('not in the same order in the two files; they are matched by ID', result.stdout)
        self.assertIn('40 samples in both the count and metadata files', result.stdout)

    def test_existing_output_directory_asks_before_overwriting(self):
        counts, metadata, outdir = self.write_inputs()
        os.mkdir(outdir)
        for answer, returncode in [('n\n', 1), ('y\n', 0)]:
            with self.subTest(answer=answer):
                result = subprocess.run(
                    [sys.executable, '-c', 'from moonstone.main import run; run()', counts, metadata, outdir],
                    capture_output=True, text=True, timeout=120, input=answer,
                )
                self.assertIn('already exists! Continue anyway', result.stdout)
                self.assertNotIn('Traceback', result.stderr)
                self.assertEqual(result.returncode, returncode)
        self.assertIn('Not overwriting...aborting.', subprocess.run(
            [sys.executable, '-c', 'from moonstone.main import run; run()', counts, metadata, outdir],
            capture_output=True, text=True, timeout=120, input='n\n').stderr)


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
