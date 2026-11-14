import os
from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.analysis.diversity.beta import BrayCurtis
from moonstone.analysis.permanova import NOT_APPLICABLE, permanova_model, run_permanova

REFERENCE = os.path.join(os.path.dirname(__file__), '..', 'reference')
PERMUTATIONS = 9999


class TestPermanovaModel(TestCase):
    """permanova_model against vegan (R script and outputs in tests/reference/permanova_adonis2.R)."""

    @classmethod
    def setUpClass(cls):
        cls.counts = pd.read_csv(os.path.join(REFERENCE, 'permanova_counts.csv'), index_col=0)
        cls.meta = pd.read_csv(os.path.join(REFERENCE, 'permanova_metadata.csv'), index_col=0)
        cls.distmat = BrayCurtis(cls.counts).beta_diversity
        cls.expected = pd.read_csv(os.path.join(REFERENCE, 'permanova_adonis2_expected.csv'))
        cls.terms = ['group', 'sex', 'age', 'bmi']

    def fit(self, meta=None, **kwargs):
        kwargs.setdefault('permutations', 999)
        return permanova_model(self.distmat, self.meta if meta is None else meta, self.terms, **kwargs)

    def test_matches_adonis2(self):
        for strata, by in [('none', 'terms'), ('none', 'margin'), ('batch', 'terms'), ('batch', 'margin')]:
            with self.subTest(strata=strata, by=by):
                expected = self.expected[(self.expected.strata == strata) & (self.expected.by == by)].set_index('term')
                result = self.fit(by=by, strata=None if strata == 'none' else 'batch', permutations=PERMUTATIONS,
                                  dispersion=False)
                self.assertEqual(result['df'].tolist(), expected['df'].tolist())
                np.testing.assert_allclose(result['R2'], expected['R2'], rtol=0, atol=1e-8)
                np.testing.assert_allclose(result['pseudo_F'], expected['F'], rtol=0, atol=1e-8)
                # both p-values are Monte Carlo estimates: allow 5 standard errors of their difference
                tolerance = 5 * np.sqrt(2) * result['p_value_se'].clip(lower=1e-4)
                self.assertTrue((abs(result['p_value'] - expected['p']) < tolerance).all(),
                                (result['p_value'].tolist(), expected['p'].tolist()))

    def test_permdisp_against_betadisper(self):
        # Measured gap on this data: F differs by 0.006 (group, betadisper F = 0.29) and 0.003 (sex, F near 0), from
        # scikit-bio leaving out the negative-eigenvalue axes. A p-value moves most when F is near 0 (0.95 vs 0.996).
        expected = pd.read_csv(os.path.join(REFERENCE, 'permanova_betadisper_expected.csv')).set_index('term')
        result = self.fit(permutations=PERMUTATIONS)
        gap = (result['permdisp_F'] - expected['F']).dropna()
        self.assertLess(gap.abs().max(), 0.01)
        np.testing.assert_allclose(result.loc[['group', 'sex'], 'permdisp_p_value'], expected['p'], atol=0.05)
        self.assertTrue((result.loc[['age', 'bmi'], 'permdisp_F'].isna()).all())

    def test_single_categorical_term_matches_skbio(self):
        single = permanova_model(self.distmat, self.meta, ['group'], permutations=999, dispersion=False)
        screen = BrayCurtis(self.counts).permanova(self.meta, 'group', permutations=999, dispersion=False)
        self.assertAlmostEqual(single.loc['group', 'pseudo_F'], screen['pseudo_F'], places=8)
        self.assertAlmostEqual(single.loc['group', 'R2'], screen['R2'], places=8)
        self.assertAlmostEqual(single.loc['group', 'p_value'], screen['p_value'], delta=5 * screen['p_value_se'] * 2)

    def test_strata_never_move_a_sample_out_of_its_block(self):
        from moonstone.analysis.permanova import _permutation_chunks
        blocks = [np.arange(0, 7), np.arange(7, 10), np.array([10])]
        for perms in _permutation_chunks(np.random.default_rng(0), blocks, 11, 250):
            for block in blocks:
                self.assertTrue(np.isin(perms[:, block], block).all())
            self.assertTrue((perms[:, 10] == 10).all())

    def test_order_seed_and_threads(self):
        first = self.fit(strata='batch', dispersion=False)
        shuffled = self.fit(self.meta.sample(frac=1, random_state=3), strata='batch', dispersion=False)
        pd.testing.assert_frame_equal(first, shuffled)
        pd.testing.assert_frame_equal(first, self.fit(strata='batch', dispersion=False))
        pd.testing.assert_frame_equal(first, self.fit(strata='batch', dispersion=False, n_jobs=2))
        other_seed = self.fit(strata='batch', dispersion=False, seed=1)
        self.assertNotEqual(first['p_value'].tolist(), other_seed['p_value'].tolist())

    def test_missing_values_are_dropped_and_counted(self):
        meta = self.meta.copy()
        meta.loc[['S01', 'S02'], 'age'] = np.nan
        meta.loc['S03', 'sex'] = ' '
        result = self.fit(meta, dispersion=False)
        self.assertEqual(result['n_samples'].unique().tolist(), [57])
        self.assertEqual(result['n_missing'].unique().tolist(), [3])
        kept = self.distmat.filter([s for s in self.distmat.ids if s not in ('S01', 'S02', 'S03')])
        complete = permanova_model(kept, self.meta, self.terms, permutations=999, dispersion=False)
        pd.testing.assert_frame_equal(result.drop(columns='n_missing'), complete.drop(columns='n_missing'))

    def test_errors(self):
        for kwargs in [dict(by='sequence'), dict(strata='group'), dict(categorical=['age', 'nope'])]:
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                self.fit(**kwargs)
        meta = self.meta.assign(twin=self.meta['age'] * 2)
        with self.assertRaisesRegex(ValueError, 'collinear'):
            permanova_model(self.distmat, meta, ['age', 'twin'], permutations=9)
        with self.assertRaisesRegex(ValueError, 'single sample'):
            self.fit(self.meta.assign(id=self.meta.index), strata='id')

    def test_categorical_codes(self):
        meta = self.meta.assign(code=self.meta['group'].map({'a': 1, 'b': 2, 'c': 3}))
        as_slope = permanova_model(self.distmat, meta, ['code'], permutations=99)
        as_group = permanova_model(self.distmat, meta, ['code'], categorical=['code'], permutations=99)
        self.assertEqual((as_slope.loc['code', 'df'], as_group.loc['code', 'df']), (1, 2))
        self.assertEqual(as_group.loc['code', 'term_type'], 'categorical')


class TestRunPermanovaModel(TestCase):
    def test_levels_and_corrected_p_values(self):
        counts = pd.read_csv(os.path.join(REFERENCE, 'permanova_counts.csv'), index_col=0)
        counts.index = pd.MultiIndex.from_arrays([['Bacteria'] * 12, counts.index], names=['kingdom', 'species'])
        meta = pd.read_csv(os.path.join(REFERENCE, 'permanova_metadata.csv'), index_col=0)
        result = run_permanova(counts, meta, ['group'], covariates=['age'], strata='batch', taxa_level='all',
                               permutations=999, read_threshold=False)
        self.assertEqual(result['level'].tolist(), ['kingdom', 'kingdom', 'species', 'species'])
        self.assertTrue(result.loc[result['level'] == 'kingdom', 'skipped'].notna().all())
        species = result[result['level'] == 'species'].set_index('variable')
        self.assertEqual(species.loc['age', 'permdisp_F'], NOT_APPLICABLE)
        self.assertEqual(species.loc['age', 'permdisp_p_value_corrected'], NOT_APPLICABLE)
        self.assertTrue(species['p_value_corrected'].between(0, 1).all())
        with self.assertRaisesRegex(ValueError, 'Name the variable'):
            run_permanova(counts, meta, covariates=['age'], taxa_level='species')
