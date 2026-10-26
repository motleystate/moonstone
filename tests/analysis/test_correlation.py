import logging

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from moonstone.analysis.correlation import (
    check_one_sample_per_subject, edges_to_matrix, feature_correlation, filter_edges, to_gephi
)
from moonstone.transformers.compositional import clr


def simulated_counts(n_features=20, n_samples=200, correlated=((0, 1, 0.8), (2, 3, -0.7)), depth=20000, seed=3):
    """Counts drawn from log-normal absolute abundances with planted correlations between pairs of features."""
    rng = np.random.default_rng(seed)
    cov = np.eye(n_features)
    for i, j, r in correlated:
        cov[i, j] = cov[j, i] = r
    log_abundance = rng.multivariate_normal(rng.normal(0, 0.3, n_features), cov, size=n_samples).T
    proportions = np.exp(log_abundance) / np.exp(log_abundance).sum(axis=0)
    counts = np.vstack([rng.multinomial(depth, proportions[:, s]) for s in range(n_samples)]).T
    return pd.DataFrame(counts, index=[f'f{i}' for i in range(n_features)],
                        columns=[f's{i}' for i in range(n_samples)])


def pair(edges, a, b):
    row = edges[(edges['feature_a'] == a) & (edges['feature_b'] == b)]
    assert len(row) == 1
    return row.iloc[0]


@pytest.fixture
def small_counts():
    rng = np.random.default_rng(0)
    return pd.DataFrame(rng.integers(0, 200, size=(5, 12)), index=[f'f{i}' for i in range(5)],
                        columns=[f's{i}' for i in range(12)])


class TestClassicCoefficients:

    @pytest.mark.parametrize('method, scipy_test', [
        ('spearman', stats.spearmanr), ('pearson', stats.pearsonr), ('kendall', stats.kendalltau)])
    def test_matches_scipy_on_clr_values(self, small_counts, method, scipy_test):
        edges = feature_correlation(small_counts, method=method, permutations=100)
        clr_df = clr(small_counts)
        for _, row in edges.iterrows():
            expected = scipy_test(clr_df.loc[row['feature_a']], clr_df.loc[row['feature_b']])
            assert row['r'] == pytest.approx(expected.statistic)

    def test_one_row_per_pair(self, small_counts):
        edges = feature_correlation(small_counts)
        assert len(edges) == 5 * 4 / 2
        assert list(edges.columns) == ['feature_a', 'feature_b', 'method', 'r', 'z', 'p_value', 'q_value', 'n_samples',
                                       'prevalence_a', 'prevalence_b']
        assert (edges['n_samples'] == 12).all()

    def test_benjamini_hochberg_by_default(self, small_counts):
        edges = feature_correlation(small_counts)
        expected = stats.false_discovery_control(edges['p_value'], method='bh')
        np.testing.assert_allclose(edges['q_value'], expected)

    def test_no_correction(self, small_counts):
        assert feature_correlation(small_counts, correction_method=None)['q_value'].isna().all()

    def test_recovers_planted_correlations(self):
        edges = feature_correlation(simulated_counts())
        assert pair(edges, 'f0', 'f1')['r'] > 0.6
        assert pair(edges, 'f2', 'f3')['r'] < -0.5
        assert pair(edges, 'f0', 'f1')['q_value'] < 0.001

    def test_permutation_null_corrects_the_clr_bias(self):
        # with 20 independent features the CLR values are correlated about -1/19 by construction
        edges = feature_correlation(simulated_counts(correlated=(), seed=11), permutations=200)
        assert edges['r'].mean() == pytest.approx(-1 / 19, abs=0.02)
        assert edges['z'].mean() == pytest.approx(0, abs=0.2)

    def test_false_discoveries_controlled_without_correlation(self):
        # textbook p-values on CLR values find a "significant" pair in about 30% of such datasets; the permutation
        # null brings it back to about 5%
        with_discovery = sum(
            (feature_correlation(simulated_counts(correlated=(), seed=seed), permutations=200)['q_value'] < 0.05).any()
            for seed in range(100, 110))
        assert with_discovery <= 2

    def test_prevalence(self):
        counts = pd.DataFrame({'s1': [0, 5, 1], 's2': [3, 5, 1], 's3': [0, 5, 1], 's4': [2, 1, 1]},
                              index=['a', 'b', 'c'])
        edges = feature_correlation(counts)
        assert pair(edges, 'a', 'b')['prevalence_a'] == 0.5
        assert pair(edges, 'a', 'b')['prevalence_b'] == 1.0

    def test_constant_feature_is_not_tested(self, small_counts, caplog):
        # two proportional genes: their CLR values are the same in every sample
        genes = pd.DataFrame([np.arange(1, 13), 2 * np.arange(1, 13)], index=['g0', 'g1'],
                             columns=small_counts.columns)
        with caplog.at_level(logging.WARNING):
            edges = feature_correlation(small_counts, genes, method='pearson')
        assert edges['r'].isna().all() and edges['q_value'].isna().all()
        assert 'constant feature' in caplog.text


class TestRho:

    def test_formula(self, small_counts):
        edges = feature_correlation(small_counts, method='rho', permutations=100)
        clr_df = clr(small_counts)
        a, b = clr_df.loc['f0'], clr_df.loc['f3']
        expected = 1 - np.var(a - b, ddof=1) / (np.var(a, ddof=1) + np.var(b, ddof=1))
        assert pair(edges, 'f0', 'f3')['r'] == pytest.approx(expected)

    def test_proportional_features_give_one(self, small_counts):
        counts = small_counts + 1
        counts.loc['f4'] = counts.loc['f0'] * 7
        edges = feature_correlation(counts, method='rho', permutations=100)
        assert pair(edges, 'f0', 'f4')['r'] == pytest.approx(1)

    def test_permutation_p_values_are_seeded(self, small_counts):
        first = feature_correlation(small_counts, method='rho', permutations=100, seed=1)
        again = feature_correlation(small_counts, method='rho', permutations=100, seed=1)
        other = feature_correlation(small_counts, method='rho', permutations=100, seed=2)
        pd.testing.assert_frame_equal(first, again)
        assert not first['p_value'].equals(other['p_value'])

    def test_strong_proportionality_is_significant(self):
        edges = feature_correlation(simulated_counts(n_features=8), method='rho', permutations=200)
        assert pair(edges, 'f0', 'f1')['p_value'] < 1e-6


class TestSparcc:

    def test_recovers_planted_correlations(self):
        edges = feature_correlation(simulated_counts(), method='sparcc', permutations=100)
        assert pair(edges, 'f0', 'f1')['r'] == pytest.approx(0.8, abs=0.1)
        assert pair(edges, 'f2', 'f3')['r'] == pytest.approx(-0.7, abs=0.1)
        others = edges[~edges['feature_a'].isin(['f0', 'f2']) & ~edges['feature_b'].isin(['f1', 'f3'])]
        assert others['r'].abs().mean() < 0.1

    def test_values_are_correlations(self, small_counts):
        edges = feature_correlation(small_counts, method='sparcc', permutations=100)
        assert edges['r'].between(-1, 1).all()

    def test_seeded(self, small_counts):
        first = feature_correlation(small_counts, method='sparcc', permutations=100, seed=4)
        again = feature_correlation(small_counts, method='sparcc', permutations=100, seed=4)
        pd.testing.assert_frame_equal(first, again)

    def test_warns_on_proportions(self, small_counts, caplog):
        with caplog.at_level(logging.WARNING):
            feature_correlation(small_counts / small_counts.sum(axis=0), method='sparcc', permutations=100)
        assert 'look like proportions' in caplog.text

    def test_zero_arguments_refused(self, small_counts):
        with pytest.raises(ValueError, match="SparCC handles zeros itself"):
            feature_correlation(small_counts, method='sparcc', pseudocount=1)


class TestTwoTables:

    @pytest.fixture
    def genes(self, small_counts):
        rng = np.random.default_rng(5)
        return pd.DataFrame(rng.integers(0, 50, size=(3, 12)), index=['g0', 'g1', 'g2'],
                            columns=small_counts.columns)

    def test_only_cross_pairs(self, small_counts, genes):
        edges = feature_correlation(small_counts, genes)
        assert len(edges) == 5 * 3
        assert set(edges['feature_a']) == set(small_counts.index)
        assert set(edges['feature_b']) == set(genes.index)

    def test_each_table_transformed_on_its_own(self, small_counts, genes):
        edges = feature_correlation(small_counts, genes, method='pearson')
        expected = stats.pearsonr(clr(small_counts).loc['f1'], clr(genes).loc['g2']).statistic
        assert pair(edges, 'f1', 'g2')['r'] == pytest.approx(expected)

    def test_samples_matched_by_name(self, small_counts, genes):
        shuffled = genes[genes.columns[::-1]]
        pd.testing.assert_frame_equal(feature_correlation(small_counts, shuffled),
                                      feature_correlation(small_counts, genes))

    def test_unshared_samples_left_out(self, small_counts, genes, caplog):
        with caplog.at_level(logging.WARNING):
            edges = feature_correlation(small_counts, genes.drop(columns=['s0', 's1']))
        assert (edges['n_samples'] == 10).all()
        assert 'left out' in caplog.text

    def test_no_shared_sample(self, small_counts, genes):
        with pytest.raises(ValueError, match="No sample name"):
            feature_correlation(small_counts, genes.set_axis([f'x{i}' for i in range(12)], axis=1))

    @pytest.mark.parametrize('method', ['rho', 'sparcc'])
    def test_single_composition_methods_refused(self, small_counts, genes, method):
        with pytest.raises(ValueError, match="can correlate two tables"):
            feature_correlation(small_counts, genes, method=method)

    def test_matrix_is_rectangular(self, small_counts, genes):
        matrix = edges_to_matrix(feature_correlation(small_counts, genes))
        assert matrix.shape == (5, 3)


class TestRepeatedSubjects:

    @pytest.fixture
    def metadata(self, small_counts):
        return pd.DataFrame({'patient': [f'p{i}' for i in range(12)], 'time': ['inclusion'] * 12},
                            index=small_counts.columns)

    def test_one_sample_per_subject_passes(self, small_counts, metadata):
        feature_correlation(small_counts, metadata=metadata, group_column='patient')

    def test_repeated_subject_stops(self, small_counts, metadata):
        metadata.loc['s5', 'patient'] = 'p0'
        with pytest.raises(ValueError, match="more than one sample"):
            feature_correlation(small_counts, metadata=metadata, group_column='patient')

    def test_missing_metadata_stops(self, small_counts, metadata):
        with pytest.raises(ValueError, match="no metadata row"):
            check_one_sample_per_subject(small_counts.columns, metadata.drop(index='s3'), 'patient')

    def test_blank_subject_stops(self, small_counts, metadata):
        metadata.loc['s3', 'patient'] = np.nan
        with pytest.raises(ValueError, match="no value"):
            check_one_sample_per_subject(small_counts.columns, metadata, 'patient')

    def test_unknown_column(self, small_counts, metadata):
        with pytest.raises(ValueError, match="not a metadata column"):
            check_one_sample_per_subject(small_counts.columns, metadata, 'subject')

    def test_needs_both_arguments(self, small_counts, metadata):
        with pytest.raises(ValueError, match="both metadata and group_column"):
            feature_correlation(small_counts, metadata=metadata)


class TestArguments:

    def test_unknown_method(self, small_counts):
        with pytest.raises(ValueError, match="method must be one of"):
            feature_correlation(small_counts, method='mic')

    def test_unknown_correction(self, small_counts):
        with pytest.raises(ValueError, match="correction_method"):
            feature_correlation(small_counts, correction_method='holm')

    def test_too_few_permutations(self, small_counts):
        with pytest.raises(ValueError, match="At least 100 permutations"):
            feature_correlation(small_counts, permutations=99)

    def test_too_few_samples(self, small_counts):
        with pytest.raises(ValueError, match="At least 3 samples"):
            feature_correlation(small_counts[['s0', 's1']])

    def test_multiindex_labels(self, small_counts):
        counts = small_counts.set_axis(pd.MultiIndex.from_tuples(
            [('Lactobacillus', 'L_crispatus'), ('Lactobacillus', 'L_iners'), ('Gardnerella', 'G_vaginalis'),
             ('Prevotella', 'P_bivia'), ('Prevotella', 'P_timonensis')], names=['genus', 'species']))
        assert set(feature_correlation(counts)['feature_a']) <= {'L_crispatus', 'L_iners', 'G_vaginalis', 'P_bivia'}
        duplicated = counts.set_axis(pd.MultiIndex.from_tuples(
            [('A', 'sp'), ('B', 'sp'), ('C', 'x'), ('D', 'y'), ('E', 'z')]))
        assert 'A;sp' in set(feature_correlation(duplicated)['feature_a'])


class TestHelpers:

    @pytest.fixture
    def edges(self):
        return pd.DataFrame({
            'feature_a': ['a', 'a', 'b'], 'feature_b': ['b', 'c', 'c'], 'method': 'spearman',
            'r': [0.8, -0.6, 0.1], 'p_value': [0.001, 0.002, 0.5], 'q_value': [0.003, 0.003, 0.5],
            'n_samples': 30, 'prevalence_a': 1.0, 'prevalence_b': 1.0})

    def test_filter_keeps_both_signs(self, edges):
        kept = filter_edges(edges)
        assert kept['r'].tolist() == [0.8, -0.6]

    def test_filter_min_abs_r(self, edges):
        assert filter_edges(edges, min_abs_r=0.7)['r'].tolist() == [0.8]

    def test_filter_without_correction_warns(self, edges, caplog):
        edges['q_value'] = np.nan
        with caplog.at_level(logging.WARNING):
            assert len(filter_edges(edges)) == 2
        assert 'uncorrected' in caplog.text

    def test_square_matrix(self, edges):
        matrix = edges_to_matrix(edges)
        assert matrix.loc['c', 'a'] == -0.6
        assert matrix.loc['a', 'c'] == -0.6
        assert matrix.loc['b', 'b'] == 1.0

    def test_gephi(self, edges):
        gephi = to_gephi(filter_edges(edges))
        assert list(gephi.columns) == ['Source', 'Target', 'Weight', 'Sign', 'r', 'q_value']
        assert gephi['Weight'].tolist() == [0.8, 0.6]
        assert gephi['Sign'].tolist() == ['positive', 'negative']


def test_too_sparse_to_permute():
    from moonstone.analysis.correlation import _shuffle_features
    with pytest.raises(ValueError, match="too sparse to permute"):
        _shuffle_features(np.array([[1.0, 0.0, 0.0]]), np.random.default_rng(0))
