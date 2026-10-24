import os
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from moonstone.analysis.permanova import (
    RESULT_COLUMNS, read_taxonomy_counts, read_threshold_filter, run_permanova, sum_to_level
)


def make_taxonomy_counts(n_per_group=15, n_species=40, shift=8.0, seed=0):
    """Counts with a 7-level taxonomy index; group B has more of the first 8 species."""
    rng = np.random.default_rng(seed)
    taxonomy = pd.MultiIndex.from_tuples(
        [('Bacteria', f'p{i % 2}', f'c{i % 3}', f'o{i % 4}', f'f{i % 6}', f'g{i % 10}', f'sp{i}')
         for i in range(n_species)],
        names=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
    )
    base = rng.gamma(1, 100, n_species) * (rng.random(n_species) < 0.7)   # some rare/absent species
    shifted = base.copy()
    shifted[:8] *= shift
    samples = [f's{i}' for i in range(2 * n_per_group)]
    counts = np.array([rng.poisson(base * rng.lognormal(0, 0.2, n_species)) for _ in range(n_per_group)]
                      + [rng.poisson(shifted * rng.lognormal(0, 0.2, n_species)) for _ in range(n_per_group)])
    counts_df = pd.DataFrame(counts.T.astype(float), index=taxonomy, columns=samples)
    metadata = pd.DataFrame({'group': ['A'] * n_per_group + ['B'] * n_per_group,
                             'random': rng.choice(['x', 'y'], 2 * n_per_group)}, index=samples)
    return counts_df, metadata


class TestReadTaxonomyCounts(TestCase):

    def test_tab_and_comma(self):
        counts_df, _ = make_taxonomy_counts()
        tmpdir = tempfile.mkdtemp()
        for sep, name in [('\t', 'counts.tsv'), (',', 'counts.csv')]:
            path = os.path.join(tmpdir, name)
            counts_df.reset_index().assign(NCBI_taxonomy_ID=range(counts_df.shape[0])).to_csv(path, sep=sep,
                                                                                              index=False)
            with self.subTest(sep=sep):
                read = read_taxonomy_counts(path)
                self.assertEqual(list(read.index.names), list(counts_df.index.names))
                self.assertEqual(list(read.columns), list(counts_df.columns))
                np.testing.assert_allclose(read.to_numpy(), counts_df.to_numpy())

    def test_no_taxonomy_column_raises(self):
        path = os.path.join(tempfile.mkdtemp(), 'flat.csv')
        pd.DataFrame({'OTU': ['a', 'b'], 's1': [1, 2]}).to_csv(path, index=False)
        with self.assertRaises(ValueError):
            read_taxonomy_counts(path)


class TestSumAndFilter(TestCase):

    def test_sum_to_level_keeps_all_reads(self):
        counts_df, _ = make_taxonomy_counts()
        genus, warnings = sum_to_level(counts_df, 'genus')
        self.assertEqual(genus.shape[0], 10)
        np.testing.assert_allclose(genus.sum().to_numpy(), counts_df.sum().to_numpy())
        self.assertEqual(warnings, [])

    def test_sum_to_level_pools_unnamed(self):
        counts_df, _ = make_taxonomy_counts()
        index = counts_df.index.to_frame()
        index.loc[index.index[:3], 'genus'] = np.nan
        counts_df.index = pd.MultiIndex.from_frame(index)
        genus, warnings = sum_to_level(counts_df, 'genus')
        self.assertIn('unclassified genus', genus.index)
        np.testing.assert_allclose(genus.sum().to_numpy(), counts_df.sum().to_numpy())
        self.assertEqual(len(warnings), 1)

    def test_sum_then_filter_order(self):
        """A genus made of rare species must be judged on its summed reads."""
        samples = [f's{i}' for i in range(10)]
        rows = {('gA', f'spA{i}'): [3] * 10 for i in range(30)}       # 30 rare species of one genus
        rows.update({('gB', 'spB'): [50] * 10, ('gC', 'spC'): [1] * 10, ('gD', 'spD'): [0] * 9 + [1]})
        counts_df = pd.DataFrame.from_dict(rows, orient='index', columns=samples)
        counts_df.index = pd.MultiIndex.from_tuples(counts_df.index, names=['genus', 'species'])
        genus, _ = sum_to_level(counts_df, 'genus')
        filtered, summary, _ = read_threshold_filter(genus, 'genus')
        self.assertIn('gA', filtered.index)    # 90 reads per sample once summed

    def test_species_has_no_knee_warning(self):
        counts_df, _ = make_taxonomy_counts()
        species, _ = sum_to_level(counts_df, 'species')
        filtered, summary, warnings = read_threshold_filter(species, 'species')
        self.assertEqual(warnings, [])
        self.assertLess(summary['taxa_kept'], summary['taxa_total'])
        self.assertEqual(summary['taxa_kept'], filtered.shape[0])

    def test_other_levels_warn(self):
        counts_df, _ = make_taxonomy_counts()
        for level in ['genus', 'family', 'phylum']:
            with self.subTest(level=level):
                level_df, _ = sum_to_level(counts_df, level)
                _, _, warnings = read_threshold_filter(level_df, level)
                self.assertTrue(any('validated at species level only' in w for w in warnings))

    def test_too_few_taxa_not_filtered(self):
        counts_df, _ = make_taxonomy_counts()
        kingdom, _ = sum_to_level(counts_df, 'kingdom')
        filtered, summary, warnings = read_threshold_filter(kingdom, 'kingdom')
        self.assertEqual(filtered.shape[0], 1)
        self.assertTrue(any('not applied' in w for w in warnings))


class TestRunPermanova(TestCase):

    def test_species_default(self):
        counts_df, metadata = make_taxonomy_counts()
        results = run_permanova(counts_df, metadata, permutations=199)
        self.assertEqual(list(results.columns), RESULT_COLUMNS)
        self.assertEqual(set(results['level']), {'species'})
        self.assertEqual(list(results['variable']), ['group', 'random'])
        group = results.set_index('variable').loc['group']
        self.assertLess(group['p_value_corrected'], 0.05)
        self.assertEqual(group['warnings'], '')

    def test_all_levels_corrected_together(self):
        counts_df, metadata = make_taxonomy_counts()
        results = run_permanova(counts_df, metadata, taxa_level='all', permutations=299)
        self.assertEqual(list(results['level'].unique()),
                         ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])
        tested = results.dropna(subset=['p_value'])
        np.testing.assert_allclose(tested['p_value_corrected'],
                                   multipletests(tested['p_value'], method='fdr_bh')[1])
        np.testing.assert_allclose(tested['permdisp_p_value_corrected'],
                                   multipletests(tested['permdisp_p_value'], method='fdr_bh')[1])
        genus_warnings = results.loc[results['level'] == 'genus', 'warnings'].iloc[0]
        self.assertIn('validated at species level only', genus_warnings)
        self.assertEqual(results.loc[results['level'] == 'species', 'warnings'].iloc[0], '')
        kingdom = results[results['level'] == 'kingdom']
        self.assertTrue((kingdom['skipped'] == 'fewer than 2 taxa, composition cannot be compared').all())
        self.assertTrue(kingdom['p_value'].isna().all())

    def test_warns_when_filter_removes_many_reads(self):
        level_df = pd.DataFrame([[100] * 10] + [[30, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20,
                                columns=[f's{i}' for i in range(10)])
        filtered, summary, warnings = read_threshold_filter(level_df, 'species')
        self.assertLess(summary['reads_pct_kept'], 90)
        self.assertTrue(any('% of reads' in w for w in warnings))

    def test_no_read_threshold(self):
        counts_df, metadata = make_taxonomy_counts()
        results = run_permanova(counts_df, metadata, read_threshold=False, permutations=49)
        row = results.iloc[0]
        self.assertEqual(row['taxa_kept'], row['taxa_total'])
        self.assertEqual(row['reads_pct_kept'], 100.0)
        self.assertEqual(row['warnings'], '')

    def test_p_value_standard_error(self):
        counts_df, metadata = make_taxonomy_counts()
        results = run_permanova(counts_df, metadata, permutations=199)
        np.testing.assert_allclose(results['p_value_se'],
                                   np.sqrt(results['p_value'] * (1 - results['p_value']) / 199))
        np.testing.assert_allclose(results['permdisp_p_value_se'],
                                   np.sqrt(results['permdisp_p_value'] * (1 - results['permdisp_p_value']) / 199))

    def test_too_few_permutations_for_correction_warned(self):
        counts_df, metadata = make_taxonomy_counts()
        with self.assertLogs('moonstone.analysis.permanova', level='WARNING') as log:
            results = run_permanova(counts_df, metadata, taxa_level='all', permutations=99)
        message = '99 permutations are too few for 10 tests'
        self.assertTrue(any(message in line for line in log.output))
        self.assertTrue(results['warnings'].str.contains(message).all())
        self.assertIn('Use at least 200 permutations', results['warnings'].iloc[0])
        # no warning without correction
        results = run_permanova(counts_df, metadata, taxa_level='all', permutations=99, correction_method=None)
        self.assertFalse(results['warnings'].str.contains('permutations are too few').any())

    def test_default_permutations(self):
        counts_df, metadata = make_taxonomy_counts()
        results = run_permanova(counts_df, metadata, columns=['group'])
        self.assertEqual(results.loc[0, 'permutations'], 9999)

    def test_parallel_matches_serial(self):
        counts_df, metadata = make_taxonomy_counts()
        serial = run_permanova(counts_df, metadata, taxa_level='all', permutations=299, n_jobs=1)
        parallel = run_permanova(counts_df, metadata, taxa_level='all', permutations=299, n_jobs=2)
        pd.testing.assert_frame_equal(serial, parallel)

    def test_selected_columns_and_seed(self):
        counts_df, metadata = make_taxonomy_counts()
        first = run_permanova(counts_df, metadata, columns=['random'], permutations=99, seed=3)
        second = run_permanova(counts_df, metadata, columns=['random'], permutations=99, seed=3)
        self.assertEqual(list(first['variable']), ['random'])
        pd.testing.assert_frame_equal(first, second)

    def test_unmatched_and_empty_samples(self):
        counts_df, metadata = make_taxonomy_counts()
        counts_df['extra'] = 1.0                 # no metadata
        counts_df['s0'] = 0.0                    # no reads
        results = run_permanova(counts_df, metadata, columns=['group'], read_threshold=False, permutations=9)
        self.assertEqual(results.loc[0, 'n_samples'], 29)
        self.assertIn('no reads left out', results.loc[0, 'warnings'])

    def test_missing_values_counted_at_every_level(self):
        counts_df, metadata = make_taxonomy_counts()
        metadata.loc[['s0', 's20'], 'group'] = [np.nan, '']
        results = run_permanova(counts_df, metadata, columns=['group', 'random'], taxa_level='all',
                                read_threshold=False, permutations=9)
        group = results[results['variable'] == 'group'].set_index('level')
        tested = group['skipped'] == ''
        self.assertTrue((group.loc[tested, 'n_missing'] == 2).all())
        self.assertTrue((group.loc[tested, 'n_samples'] == 28).all())
        random = results[(results['variable'] == 'random') & (results['skipped'] == '')]
        self.assertTrue((random['n_missing'] == 0).all())

    def test_unknown_level_raises(self):
        counts_df, metadata = make_taxonomy_counts()
        with self.assertRaises(ValueError):
            run_permanova(counts_df, metadata, taxa_level='strain')
