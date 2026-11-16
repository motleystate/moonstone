import contextlib
import io
import os
import subprocess
import sys
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone import __version__, examples
from moonstone.analysis.permanova import read_taxonomy_counts
from moonstone.filtering.read_threshold_analysis import load_species_table
from moonstone.normalization import table
from moonstone.normalization.counts.geometric_mean import GeometricMeanNormalization
from moonstone.plot import PlotTaxonomyCounts
from tests.reads.helpers import run_and_measure

LEVELS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']


class TestNormalizeTable(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.source = examples.path('taxa_counts.csv')
        self.output = os.path.join(self.tmp, 'normalized.csv')
        self.counts = pd.read_csv(self.source).set_index(LEVELS + ['NCBI_taxonomy_ID'])

    def call(self, *arguments):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            try:
                code = table.main(list(arguments))
            except SystemExit as stop:
                code = stop.code
        return code, output.getvalue()

    def test_geometric_mean_equals_the_existing_class(self):
        expected = GeometricMeanNormalization(self.counts.copy()).normalized_df
        pd.testing.assert_frame_equal(table.normalize(self.counts, 'geometric_mean'), expected)

    def test_relative_abundance_sums_to_100(self):
        result = table.normalize(self.counts, 'relative_abundance')
        np.testing.assert_allclose(result.sum(), 100)
        pd.testing.assert_frame_equal(result, PlotTaxonomyCounts(self.counts).compute_relative_abundance_dataframe())

    def test_unknown_method(self):
        with self.assertRaisesRegex(ValueError, 'geometric_mean, relative_abundance'):
            table.normalize(self.counts, 'total_counts')

    def test_identifying_columns_are_kept(self):
        for method in table.METHODS:
            with self.subTest(method=method):
                output = os.path.join(self.tmp, method + '.csv')
                summary = table.normalize_table(self.source, output, method)
                self.assertEqual((summary['samples'], summary['features']), (60, 8))
                written, source = pd.read_csv(output), pd.read_csv(self.source)
                pd.testing.assert_frame_equal(written[LEVELS + ['NCBI_taxonomy_ID']], source[LEVELS + ['NCBI_taxonomy_ID']])
                self.assertEqual(list(written.columns), list(source.columns))

    def test_plain_feature_table(self):
        table.normalize_table(examples.path('resistance_genes.csv'), self.output, 'relative_abundance')
        written = pd.read_csv(self.output)
        self.assertEqual(list(written['gene']), ['tetM', 'ermB', 'lnuA', 'blaTEM'])
        np.testing.assert_allclose(written.drop(columns='gene').sum(), 100)

    def test_tab_separated_input(self):
        table.normalize_table(examples.path('species_counts.tsv'), self.output, 'geometric_mean')
        self.assertEqual(read_taxonomy_counts(self.output).shape, (8, 60))

    def test_output_is_read_by_the_analyses(self):
        table.normalize_table(self.source, self.output, 'geometric_mean')
        self.assertEqual(read_taxonomy_counts(self.output).shape, (8, 60))
        self.assertEqual(load_species_table(self.output).shape, (8, 60))

    def test_method_is_required(self):
        code, output = self.call(self.source, self.output)
        self.assertEqual(code, 2)
        self.assertIn('--method', output)
        self.assertFalse(os.path.exists(self.output))

    def test_sample_without_reads_is_named(self):
        source = os.path.join(self.tmp, 'empty_sample.csv')
        with open(source, 'w') as handle:
            handle.write('gene,S1,S2,S3\ntetM,1,0,3\nermB,4,0,6\n')
        for method in table.METHODS:
            with self.subTest(method=method), self.assertRaisesRegex(ValueError, 'S2'):
                table.normalize_table(source, os.path.join(self.tmp, method + '.csv'), method)

    def test_text_in_a_sample_column_is_named(self):
        source = os.path.join(self.tmp, 'text.csv')
        with open(source, 'w') as handle:
            handle.write('gene,S1,S2\ntetM,1,many\n')
        with self.assertRaisesRegex(ValueError, 'S2'):
            table.normalize_table(source, self.output, 'relative_abundance')

    def test_methods_file_states_what_was_done_and_what_it_suits(self):
        table.normalize_table(self.source, self.output, 'geometric_mean')
        with open(self.output + '.methods.txt') as handle:
            methods = handle.read()
        for expected in ['taxa_counts.csv', 'geometric mean', 'non-zero in at least 80%', 'of 8', 'every Moonstone analysis',
                         __version__]:
            self.assertIn(expected, methods)
        other = os.path.join(self.tmp, 'relative.csv')
        table.normalize_table(self.source, other, 'relative_abundance')
        with open(other + '.methods.txt') as handle:
            self.assertIn('not for PERMANOVA or read threshold', handle.read())

    def test_relaxed_zero_filter_is_recorded(self):
        """With sparse data the existing class lowers its zero filter by itself: the value used is written down."""
        rng = np.random.default_rng(0)
        counts = pd.DataFrame(rng.poisson(0.5, (60, 12)) * rng.integers(1, 50, (60, 12)),
                              index=[f'taxon{i}' for i in range(60)], columns=[f'S{i}' for i in range(12)])
        counts.iloc[1, :9] = 7       # the most widespread taxon is non-zero in 75% of the samples: below the default 80%
        source = os.path.join(self.tmp, 'sparse.csv')
        counts.rename_axis('taxon').to_csv(source)
        summary = table.normalize_table(source, self.output, 'geometric_mean')
        self.assertLess(summary['zero_threshold'], 80)
        with open(self.output + '.methods.txt') as handle:
            self.assertIn(f"non-zero in at least {summary['zero_threshold']}%", handle.read())

    def test_never_overwrites(self):
        table.normalize_table(self.source, self.output, 'relative_abundance')
        code, output = self.call(self.source, self.output, '--method', 'relative_abundance')
        self.assertEqual(code, 1)
        self.assertIn('already exists', output)

    def test_output_may_be_a_directory(self):
        code, output = self.call(self.source, self.tmp, '--method', 'relative_abundance')
        self.assertEqual(code, 0, output)
        self.assertTrue(os.path.exists(os.path.join(self.tmp, 'counts_relative_abundance.csv')))

    def test_does_not_need_the_reads_package(self):
        code = "import sys; sys.modules['moonstone.reads'] = None; import moonstone.normalization.table"
        self.assertEqual(subprocess.run([sys.executable, '-c', code]).returncode, 0)

    def test_large_table(self):
        rng = np.random.default_rng(0)
        counts = pd.DataFrame(rng.poisson(3, (4000, 2000)), columns=[f'S{i}' for i in range(2000)])
        counts.insert(0, 'species', [f'species {i}' for i in range(4000)])
        source = os.path.join(self.tmp, 'large.csv')
        counts.to_csv(source, index=False)
        for method in table.METHODS:
            with self.subTest(method=method):
                code, output, peak_kb = run_and_measure('moonstone.normalization.table', source,
                                                        os.path.join(self.tmp, method + '.csv'), '--method', method)
                self.assertEqual(code, 0, output)
                self.assertTrue(0 < peak_kb < 2 * 1024 * 1024, peak_kb)       # under 2 GB
