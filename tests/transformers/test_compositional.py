import logging

import numpy as np
import pandas as pd
import pytest
from scipy import stats
from skbio.stats import composition

from moonstone.transformers.compositional import clr, multiplicative_replacement, rclr


@pytest.fixture
def counts():
    # features as rows, samples as columns, with a MultiIndex like Moonstone's taxonomy tables
    index = pd.MultiIndex.from_tuples(
        [('Lactobacillus', 'crispatus'), ('Lactobacillus', 'iners'), ('Gardnerella', 'vaginalis'),
         ('Prevotella', 'bivia')],
        names=['genus', 'species'],
    )
    return pd.DataFrame(
        {'s1': [120, 30, 0, 5], 's2': [0, 80, 40, 10], 's3': [60, 60, 15, 0], 's4': [3, 0, 200, 25]},
        index=index,
    )


@pytest.fixture
def positive_counts():
    rng = np.random.default_rng(1)
    return pd.DataFrame(rng.integers(1, 500, size=(6, 8)), columns=[f's{i}' for i in range(8)])


class TestClr:

    def test_matches_scikit_bio_on_positive_values(self, positive_counts):
        expected = composition.clr(positive_counts.to_numpy().T).T
        np.testing.assert_allclose(clr(positive_counts).to_numpy(), expected)
        np.testing.assert_allclose(clr(positive_counts, zeros='pseudocount', pseudocount=0.5).to_numpy(),
                                   composition.clr((positive_counts + 0.5).to_numpy().T).T)

    def test_keeps_shape_labels_and_input(self, counts):
        before = counts.copy()
        result = clr(counts)
        assert result.shape == counts.shape
        assert result.index.equals(counts.index)
        assert result.columns.equals(counts.columns)
        pd.testing.assert_frame_equal(counts, before)

    def test_each_sample_sums_to_zero(self, counts):
        np.testing.assert_allclose(clr(counts).sum(axis=0), 0, atol=1e-12)

    def test_unchanged_by_scaling_a_sample(self, counts):
        scaled = counts.astype(float).copy()
        scaled['s2'] *= 0.001  # e.g. relative abundance instead of counts
        pd.testing.assert_frame_equal(clr(scaled), clr(counts))

    def test_scaling_a_feature_does_not_change_correlations(self, positive_counts):
        scaled = positive_counts.astype(float).copy()
        scaled.iloc[2] *= 1 / 1500  # e.g. RPKM gene-length scaling
        before = clr(positive_counts).T.corr(method='spearman')
        after = clr(scaled).T.corr(method='spearman')
        pd.testing.assert_frame_equal(before, after)
        before = clr(positive_counts).T.corr(method='pearson')
        after = clr(scaled).T.corr(method='pearson')
        np.testing.assert_allclose(before, after, atol=1e-12)

    def test_does_not_depend_on_sample_order(self, counts):
        shuffled = counts[['s3', 's1', 's4', 's2']]
        pd.testing.assert_frame_equal(clr(shuffled)[counts.columns], clr(counts))

    def test_removes_closure_artefact(self):
        # independent absolute abundances, then closed to proportions: proportions of the minor features are
        # negatively correlated with the dominant feature by construction; CLR removes most of it
        rng = np.random.default_rng(0)
        absolute = pd.DataFrame(rng.lognormal(mean=[[5], [1], [1], [1], [1]], sigma=1, size=(5, 300)))
        proportions = absolute / absolute.sum(axis=0)
        raw_r = stats.spearmanr(proportions.iloc[0], proportions.iloc[1]).statistic
        abs_r = stats.spearmanr(absolute.iloc[0], absolute.iloc[1]).statistic
        assert raw_r < -0.2
        assert abs(abs_r) < 0.15
        clr_df = clr(proportions)
        r_minor = stats.spearmanr(clr_df.iloc[1], clr_df.iloc[2]).statistic
        assert abs(r_minor) < 0.3

    def test_pseudocount_is_required_and_positive(self, counts):
        with pytest.raises(ValueError, match="explicit pseudocount"):
            clr(counts, zeros='pseudocount')
        with pytest.raises(ValueError, match="explicit pseudocount"):
            clr(counts, zeros='pseudocount', pseudocount=0)

    def test_arguments_must_match_the_zero_method(self, counts):
        with pytest.raises(ValueError, match="only used with zeros='pseudocount'"):
            clr(counts, pseudocount=1)
        with pytest.raises(ValueError, match="only used with zeros='multiplicative'"):
            clr(counts, zeros='pseudocount', pseudocount=1, delta=0.01)
        with pytest.raises(ValueError, match="zeros must be one of"):
            clr(counts, zeros='drop')


class TestMultiplicativeReplacement:

    def test_no_zeros_and_samples_sum_to_one(self, counts):
        result = multiplicative_replacement(counts)
        assert (result > 0).all().all()
        np.testing.assert_allclose(result.sum(axis=0), 1)

    def test_keeps_ratios_between_non_zero_values(self, counts):
        result = multiplicative_replacement(counts)
        assert result.loc[('Lactobacillus', 'crispatus'), 's1'] / result.loc[('Prevotella', 'bivia'), 's1'] == \
            pytest.approx(120 / 5)

    def test_default_delta_is_below_every_observed_value(self, counts):
        proportions = counts / counts.sum(axis=0)
        result = multiplicative_replacement(counts)
        for sample in counts.columns:
            smallest = proportions.loc[proportions[sample] > 0, sample].min()
            replaced = result.loc[counts[sample] == 0, sample]
            np.testing.assert_allclose(replaced, 0.65 * smallest)

    def test_fixed_delta_matches_scikit_bio(self, counts):
        expected = composition.multi_replace(counts.to_numpy().T, delta=0.01).T
        np.testing.assert_allclose(multiplicative_replacement(counts, delta=0.01).to_numpy(), expected)

    def test_delta_out_of_range(self, counts):
        with pytest.raises(ValueError, match="between 0 and 1"):
            multiplicative_replacement(counts, delta=0)
        two_zeros = pd.DataFrame({'s1': [10, 0, 0], 's2': [1, 2, 3]})
        with pytest.raises(ValueError, match="too large"):
            multiplicative_replacement(two_zeros, delta=0.5)


class TestRclr:

    def test_zeros_stay_missing(self, counts):
        result = rclr(counts)
        pd.testing.assert_frame_equal(result.isna(), counts == 0)

    def test_matches_scikit_bio(self, counts):
        expected = composition.rclr(counts.to_numpy().T.astype(float)).T
        np.testing.assert_allclose(rclr(counts).to_numpy(), expected)

    def test_non_zero_values_centred(self, counts):
        np.testing.assert_allclose(rclr(counts).sum(axis=0), 0, atol=1e-12)


class TestInputChecks:

    @pytest.mark.parametrize('transform', [clr, rclr, multiplicative_replacement])
    def test_missing_values_are_refused(self, counts, transform):
        broken = counts.astype(float)
        broken.iloc[0, 1] = np.nan
        with pytest.raises(ValueError, match="missing values"):
            transform(broken)

    @pytest.mark.parametrize('transform', [clr, rclr])
    def test_negative_values_are_refused(self, counts, transform):
        broken = counts.copy()
        broken.iloc[0, 0] = -1
        with pytest.raises(ValueError, match="negative"):
            transform(broken)

    def test_empty_sample_is_refused(self, counts):
        broken = counts.copy()
        broken['s5'] = 0
        with pytest.raises(ValueError, match="only zeros"):
            clr(broken)

    def test_text_column_is_refused(self, counts):
        broken = counts.copy()
        broken['note'] = 'x'
        with pytest.raises(ValueError, match="not numeric"):
            clr(broken)

    def test_not_a_dataframe(self, counts):
        with pytest.raises(TypeError):
            clr(counts.to_numpy())

    def test_all_zero_feature_warns_but_is_kept(self, counts, caplog):
        with_empty = pd.concat([counts, pd.DataFrame([[0, 0, 0, 0]], columns=counts.columns,
                                                     index=pd.MultiIndex.from_tuples([('X', 'y')]))])
        with caplog.at_level(logging.WARNING):
            result = clr(with_empty)
        assert "zero in every sample" in caplog.text
        assert result.shape == with_empty.shape
