import contextlib
import io
import os
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone import __version__, examples
from moonstone.analysis.compositional_pipeline import read_feature_table
from moonstone.analysis.permanova import read_taxonomy_counts
from moonstone.filtering.read_threshold_analysis import load_species_table
from moonstone.reads import cli
from moonstone.reads.table import read_moonbeam_table, to_moonstone_table
from moonstone.validation import match_samples
from tests.reads.helpers import LINEAGE, run_and_measure


class TestTable(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.source = examples.path('moonbeam_bracken.tsv')
        self.output = os.path.join(self.tmp, 'species_counts.csv')

    def write(self, text, name='table.tsv'):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_counts_are_unchanged(self):
        summary = to_moonstone_table(self.source, self.output)
        self.assertEqual(summary, {'samples': 60, 'taxa': 8})
        written = pd.read_csv(self.output).set_index('NCBI_taxonomy_ID')
        expected = pd.read_csv(examples.path('taxa_counts.csv')).set_index('NCBI_taxonomy_ID')
        samples = [f'S{i:03d}' for i in range(1, 61)]
        levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        self.assertEqual(list(written.columns), levels + samples)
        pd.testing.assert_frame_equal(written.loc[expected.index, samples], expected[samples])
        self.assertEqual(written.loc[47770, 'species'], 'Lactobacillus_crispatus')

    def test_both_lineage_headers(self):
        with open(self.source) as handle:
            sunbeam = self.write(handle.read().replace('\ttaxonomy\n', '\tConsensus Lineage\n'))
        pd.testing.assert_frame_equal(read_moonbeam_table(sunbeam), read_moonbeam_table(self.source))

    def test_sample_suffixes_are_removed(self):
        table = self.write('# Constructed from biom file\n#OTU ID\tA.brackenCorrected_taxa\tB.taxa\tC\ttaxonomy\n'
                           f'562\t1.0\t2.0\t3.0\t{LINEAGE}\n')
        self.assertEqual(list(read_moonbeam_table(table).columns), ['NCBI_taxonomy_ID', 'A', 'B', 'C'])

    def test_duplicate_sample_after_cleaning_is_an_error(self):
        table = self.write('# Constructed from biom file\n#OTU ID\tA.brackenCorrected_taxa\tA.taxa\ttaxonomy\n'
                           f'562\t1.0\t2.0\t{LINEAGE}\n')
        with self.assertRaisesRegex(ValueError, 'A'):
            read_moonbeam_table(table)

    def test_a_table_without_lineage_is_refused(self):
        table = self.write('gene,S1,S2\ntetM,1,2\n', 'genes.csv')
        with self.assertRaisesRegex(ValueError, 'taxonomy'):
            read_moonbeam_table(table)

    def test_output_is_read_by_the_analyses(self):
        to_moonstone_table(self.source, self.output)
        metadata = pd.read_csv(examples.path('metadata.csv'), index_col=0, dtype={0: str})
        for name, table in [('permanova', read_taxonomy_counts(self.output)),
                            ('correlation', read_feature_table(self.output)[0]),
                            ('read threshold', load_species_table(self.output))]:
            with self.subTest(reader=name):
                self.assertEqual(table.shape, (8, 60))
                self.assertEqual(len(match_samples(table.columns, metadata.index).matched), 60)

    def test_methods_file(self):
        to_moonstone_table(self.source, self.output)
        with open(self.output + '.methods.txt') as handle:
            methods = handle.read()
        for expected in ['moonbeam_bracken.tsv', 'not normalized', '60 samples', '8 taxa', __version__]:
            self.assertIn(expected, methods)

    def test_never_overwrites(self):
        to_moonstone_table(self.source, self.output)
        with self.assertRaises(FileExistsError):
            to_moonstone_table(self.source, self.output)
        with open(self.source) as handle:
            self.assertIn('S001.brackenCorrected_taxa', handle.read())       # the input is left as it was

    def test_command(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            code = cli.main(['table', self.source, self.output])
            again = cli.main(['table', self.source, self.output])
            missing = cli.main(['table', '/nonexistent.tsv', self.output + '2'])
        self.assertEqual((code, again, missing), (0, 1, 1))
        self.assertIn('60 samples', output.getvalue())
        self.assertIn('already exists', output.getvalue())

    def test_large_table(self):
        """2,000 samples x 4,000 taxa, the size of a real cohort: memory stays ordinary."""
        rng = np.random.default_rng(0)
        table = pd.DataFrame(rng.poisson(3, (4000, 2000)).astype(float),
                             columns=[f'S{i}.brackenCorrected_taxa' for i in range(2000)])
        table.insert(0, '#OTU ID', range(1, 4001))
        table['taxonomy'] = [f'k__Bacteria; p__P; c__C; o__O; f__F; g__G{i}; s__s{i}' for i in range(4000)]
        source = os.path.join(self.tmp, 'large.tsv')
        with open(source, 'w') as handle:
            handle.write('# Constructed from biom file\n')
            table.to_csv(handle, sep='\t', index=False)
        code, output, peak_kb = run_and_measure('moonstone.reads.table', source, self.output)
        self.assertEqual(code, 0, output)
        self.assertTrue(0 < peak_kb < 2 * 1024 * 1024, peak_kb)       # under 2 GB
        self.assertEqual(pd.read_csv(self.output, nrows=5).shape, (5, 2008))


class TestTableOutputFolder(TestCase):

    def test_output_may_be_a_folder(self):
        folder = tempfile.mkdtemp()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            code = cli.main(['table', examples.path('moonbeam_bracken.tsv'), folder])
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(os.path.join(folder, 'counts_raw.csv')))
        self.assertTrue(os.path.exists(os.path.join(folder, 'counts_raw.csv.methods.txt')))
