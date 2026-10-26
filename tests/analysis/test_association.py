import logging

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from moonstone.analysis.association import ancombc_check, feature_association
from moonstone.transformers.compositional import clr


def simulated(n=80, d=15, seed=0, bv_effect=1.5, age_effect=0.08):
    """Counts where f0 is higher with bv (vs healthy) and f1 rises with age; no other feature is associated."""
    rng = np.random.default_rng(seed)
    samples = [f's{i}' for i in range(n)]
    metadata = pd.DataFrame({
        'disease': rng.choice(['healthy', 'bv', 'candida'], n),
        'age': rng.normal(35, 8, n),
        'smoker': rng.choice([True, False], n),
    }, index=samples)
    log_abundance = rng.normal(3, 1, (d, n))
    log_abundance[0] += bv_effect * (metadata['disease'] == 'bv').to_numpy()
    log_abundance[1] += age_effect * (metadata['age'] - 35).to_numpy()
    proportions = np.exp(log_abundance) / np.exp(log_abundance).sum(axis=0)
    counts = np.stack([rng.multinomial(20000, proportions[:, j]) for j in range(n)]).T
    return pd.DataFrame(counts, index=[f'f{i}' for i in range(d)], columns=samples), metadata


@pytest.fixture
def data():
    return simulated()


class TestFeatureAssociation:

    def test_matches_statsmodels_ols(self, data):
        counts, metadata = data
        results = feature_association(counts, metadata, ['disease', 'age'], covariates='smoker',
                                      reference={'disease': 'healthy'}, correction_method=None)
        transformed = clr(counts)
        design = pd.DataFrame({
            'const': 1.0,
            'bv': (metadata['disease'] == 'bv').astype(float),
            'candida': (metadata['disease'] == 'candida').astype(float),
            'age': metadata['age'],
            'smoker': metadata['smoker'].astype(float),
        })
        for feature in ['f0', 'f1', 'f7']:
            fit = sm.OLS(transformed.loc[feature], design).fit()
            ci = fit.conf_int()
            for variable, level, term in [('disease', 'bv', 'bv'), ('disease', 'candida', 'candida'),
                                          ('age', None, 'age')]:
                row = results[(results['feature'] == feature) & (results['variable'] == variable)
                              & ((results['level'] == level) if level else results['level'].isna())].iloc[0]
                assert row['coef'] == pytest.approx(fit.params[term])
                assert row['std_error'] == pytest.approx(fit.bse[term])
                assert row['t_value'] == pytest.approx(fit.tvalues[term])
                assert row['p_value'] == pytest.approx(fit.pvalues[term])
                assert row['ci_low'] == pytest.approx(ci.loc[term, 0])
                assert row['ci_high'] == pytest.approx(ci.loc[term, 1])

    def test_finds_planted_associations_only(self, data):
        counts, metadata = data
        results = feature_association(counts, metadata, ['disease', 'age'], covariates='smoker',
                                      reference={'disease': 'healthy'})
        significant = results[results['q_value'] < 0.05]
        found = set(zip(significant['feature'], significant['variable'], significant['level'].fillna('')))
        assert found == {('f0', 'disease', 'bv'), ('f1', 'age', '')}
        bv = significant[significant['feature'] == 'f0'].iloc[0]
        assert bv['coef'] == pytest.approx(1.5, abs=0.6)

    def test_false_discoveries_are_controlled(self):
        datasets_with_discovery = 0
        for seed in range(20):
            counts, metadata = simulated(seed=seed, bv_effect=0, age_effect=0)
            results = feature_association(counts, metadata, ['disease', 'age'])
            datasets_with_discovery += (results['q_value'] < 0.05).any()
        assert datasets_with_discovery <= 3

    def test_only_variables_are_reported_and_corrected_together(self, data):
        counts, metadata = data
        results = feature_association(counts, metadata, 'disease', covariates=['age', 'smoker'])
        assert set(results['variable']) == {'disease'}
        assert len(results) == 15 * 2
        assert (results['q_value'] >= results['p_value']).all()
        bonferroni = feature_association(counts, metadata, 'disease', correction_method='bonferroni')
        np.testing.assert_allclose(bonferroni['q_value'], np.minimum(bonferroni['p_value'] * 30, 1))

    def test_default_reference_and_categorical_order(self, data):
        counts, metadata = data
        results = feature_association(counts, metadata, 'disease')
        assert set(results['reference']) == {'bv'}  # first in sorted order
        assert set(results['level']) == {'candida', 'healthy'}
        ordered = metadata.assign(disease=pd.Categorical(metadata['disease'], ['healthy', 'bv', 'candida']))
        results = feature_association(counts, ordered, 'disease')
        assert set(results['reference']) == {'healthy'}

    def test_boolean_is_categorical(self, data):
        counts, metadata = data
        results = feature_association(counts, metadata, 'smoker')
        assert set(results['reference']) == {'False'}
        assert set(results['level']) == {'True'}

    def test_missing_metadata_is_dropped_and_reported(self, data, caplog):
        counts, metadata = data
        metadata = metadata.copy()
        metadata.loc['s3', 'age'] = np.nan
        metadata = metadata.drop(index='s5')
        with caplog.at_level(logging.WARNING):
            results = feature_association(counts, metadata, 'age')
        assert "no metadata row" in caplog.text
        assert "missing value" in caplog.text
        assert set(results['n_samples']) == {78}
        assert set(results['n_dropped']) == {2}
        kept = [s for s in counts.columns if s not in ('s3', 's5')]
        expected = feature_association(counts[kept], metadata.loc[kept], 'age')
        pd.testing.assert_series_equal(results['coef'], expected['coef'])

    def test_samples_matched_by_name(self, data):
        counts, metadata = data
        shuffled = metadata.sample(frac=1, random_state=1)
        pd.testing.assert_frame_equal(feature_association(counts, shuffled, 'age'),
                                      feature_association(counts, metadata, 'age'))

    def test_standardize(self, data):
        counts, metadata = data
        raw = feature_association(counts, metadata, 'age')
        scaled = feature_association(counts, metadata, 'age', standardize=True)
        np.testing.assert_allclose(scaled['coef'], raw['coef'] * metadata['age'].std())
        np.testing.assert_allclose(scaled['p_value'], raw['p_value'])

    def test_categorical_codes(self, data, caplog):
        counts, metadata = data
        coded = metadata.assign(group=metadata['disease'].map({'healthy': 1, 'bv': 2, 'candida': 3}))
        with caplog.at_level(logging.WARNING):
            as_number = feature_association(counts, coded, 'group')
        assert "group codes" in caplog.text
        assert as_number['level'].isna().all()
        as_groups = feature_association(counts, coded, 'group', categorical='group')
        assert set(as_groups['level']) == {'2', '3'}
        reference = feature_association(counts, coded, 'group', categorical='group', reference={'group': '2'})
        assert set(reference['reference']) == {'2'}

    def test_constant_feature_is_not_tested(self, data, caplog):
        counts, metadata = data
        values = np.ones((3, 80))
        values[0] *= 5
        constant = pd.DataFrame(values, index=['a', 'b', 'c'], columns=counts.columns)
        with caplog.at_level(logging.WARNING):
            results = feature_association(constant, metadata, 'age')
        assert results['p_value'].isna().all()
        assert results['q_value'].isna().all()
        assert "fits the feature perfectly" in caplog.text

    def test_repeated_subjects_are_refused(self, data):
        counts, metadata = data
        metadata = metadata.assign(patient=[f'p{i // 2}' for i in range(80)])
        with pytest.raises(ValueError, match="more than one sample"):
            feature_association(counts, metadata, 'age', group_column='patient')
        metadata['patient'] = [f'p{i}' for i in range(80)]
        feature_association(counts, metadata, 'age', group_column='patient')

    def test_multiindex_feature_labels(self, data):
        counts, metadata = data
        counts.index = pd.MultiIndex.from_tuples([('Genus', f'species_{i}') for i in range(15)])
        results = feature_association(counts, metadata, 'age')
        assert results['feature'].iloc[0] == 'species_0'


class TestModelErrors:

    @pytest.mark.parametrize('kwargs, message', [
        (dict(variables=[]), "at least one variable"),
        (dict(variables='nope'), "Not metadata columns"),
        (dict(variables='age', covariates='age'), "listed twice"),
        (dict(variables='age', categorical='disease'), "not in the model"),
        (dict(variables='age', reference={'age': 30}), "is numeric"),
        (dict(variables='disease', reference={'disease': 'flu'}), "is not a value"),
        (dict(variables='age', correction_method='holm'), "correction_method"),
        (dict(variables='age', zeros='drop'), "zeros must be one of"),
    ])
    def test_argument_errors(self, data, kwargs, message):
        counts, metadata = data
        with pytest.raises(ValueError, match=message):
            feature_association(counts, metadata, **kwargs)

    def test_single_value(self, data):
        counts, metadata = data
        with pytest.raises(ValueError, match="single value"):
            feature_association(counts, metadata.assign(disease='bv'), 'disease')
        with pytest.raises(ValueError, match="single value"):
            feature_association(counts, metadata.assign(age=30.0), 'age')

    def test_collinear(self, data):
        counts, metadata = data
        metadata = metadata.assign(age_months=metadata['age'] * 12)
        with pytest.raises(ValueError, match="collinear"):
            feature_association(counts, metadata, 'age', covariates='age_months')

    def test_too_few_samples(self, data):
        counts, metadata = data
        with pytest.raises(ValueError, match="too few"):
            feature_association(counts.iloc[:, :3], metadata, ['disease', 'age'])

    def test_duplicate_metadata_rows(self, data):
        counts, metadata = data
        with pytest.raises(ValueError, match="more than once in the metadata"):
            feature_association(counts, pd.concat([metadata, metadata.iloc[:1]]), 'age')


class TestAncombcCheck:

    def test_agrees_on_planted_associations(self, data):
        counts, metadata = data
        results = ancombc_check(counts, metadata, ['disease', 'age'], covariates='smoker',
                                reference={'disease': 'healthy'})
        assert list(results.columns) == ['feature', 'variable', 'level', 'reference', 'log_fold_change',
                                         'std_error', 'p_value', 'q_value', 'n_samples', 'n_dropped']
        assert len(results) == 15 * 3
        top = results.sort_values('p_value').head(2)
        assert set(zip(top['feature'], top['variable'])) == {('f0', 'disease'), ('f1', 'age')}
        bv = results[(results['feature'] == 'f0') & (results['level'] == 'bv')].iloc[0]
        assert bv['q_value'] < 0.05
        assert bv['log_fold_change'] == pytest.approx(1.5, abs=0.6)
        assert results['q_value'].isna().equals(results['p_value'].isna())

    def test_uses_same_samples(self, data):
        counts, metadata = data
        metadata = metadata.copy()
        metadata.loc['s3', 'age'] = np.nan
        results = ancombc_check(counts, metadata, 'age')
        assert set(results['n_samples']) == {79}
        assert set(results['n_dropped']) == {1}
