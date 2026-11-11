
import pandas as pd
import pytest

from moonstone.analysis import compositional_pipeline as pipeline


@pytest.fixture
def table():
    return pd.DataFrame({'s1': [10, 0, 5, 0], 's2': [20, 1, 0, 0], 's3': [30, 2, 0, 1], 's4': [40, 0, 0, 0]},
                        index=['a', 'b', 'c', 'd'])


class TestReadFeatureTable:

    def test_taxonomy_table_is_summed(self, tmp_path):
        df = pd.DataFrame({'genus': ['G1', 'G1', 'G2'], 'species': ['G1 x', 'G1 y', 'G2 z'],
                           's1': [1, 2, 3], 's2': [4, 5, 6]})
        path = tmp_path / 'species.tsv'
        df.to_csv(path, sep='\t', index=False)
        table, notes = pipeline.read_feature_table(str(path), 'genus')
        assert table.loc['G1'].tolist() == [3, 9]
        assert notes == ['the counts were summed to genus level']

    def test_plain_table(self, tmp_path):
        path = tmp_path / 'genes.csv'
        pd.DataFrame({'gene': ['g1', 'g2'], '01': [1, 2], '02': [3, 4]}).to_csv(path, index=False)
        table, notes = pipeline.read_feature_table(str(path))
        assert list(table.columns) == ['01', '02']
        assert list(table.index) == ['g1', 'g2']
        assert notes == []

    @pytest.mark.parametrize('content, message', [
        ('gene,s1,s2\ng1,1,x\n', 'must be numeric'),
        ('gene,s1,s2\ng1,1,\n', 'blank cell'),
        ('gene,s1,s2\ng1,1,2\ng1,3,4\n', 'more than once'),
    ])
    def test_plain_table_errors(self, tmp_path, content, message):
        path = tmp_path / 'genes.csv'
        path.write_text(content)
        with pytest.raises(ValueError, match=message):
            pipeline.read_feature_table(str(path))


class TestSelectSamples:

    def test_selection(self):
        metadata = pd.DataFrame({'Type': ['inclusion', 'delivery', 'inclusion', None], 'site': [1, 1, 2, 1]},
                                index=['a', 'b', 'c', 'd'])
        selections = [pipeline.parse_selection('Type=inclusion'), pipeline.parse_selection(' site = 1 ')]
        selected, notes = pipeline.select_samples(metadata, selections)
        assert list(selected.index) == ['a']
        assert notes == ['Type = inclusion (2 of 4 samples)', 'site = 1 (1 of 2 samples)']

    @pytest.mark.parametrize('text', ['Type', '=x', 'Type='])
    def test_bad_selection_text(self, text):
        with pytest.raises(ValueError, match='COLUMN=VALUE'):
            pipeline.parse_selection(text)

    def test_unknown_column_or_value(self):
        metadata = pd.DataFrame({'Type': ['inclusion']}, index=['a'])
        with pytest.raises(ValueError, match='not a metadata column'):
            pipeline.select_samples(metadata, [('Visit', '1')])
        with pytest.raises(ValueError, match="Values found: \\['inclusion'\\]"):
            pipeline.select_samples(metadata, [('Type', 'delivery')])


class TestFilterFeatures:

    def test_prevalence(self, table):
        filtered, note = pipeline.filter_features(table, min_prevalence=0.5)
        assert list(filtered.index) == ['a', 'b']
        assert note == 'features present in at least 50% of samples were kept (2 of 4)'

    def test_top(self, table):
        filtered, note = pipeline.filter_features(table, min_prevalence=0.25, top=2)
        assert list(filtered.index) == ['a', 'c']
        assert 'then the 2 with the highest mean relative abundance (2 of 4)' in note

    def test_errors(self, table):
        with pytest.raises(ValueError, match='between 0 and 1'):
            pipeline.filter_features(table, min_prevalence=10)
        with pytest.raises(ValueError, match='at least 2'):
            pipeline.filter_features(table, top=1)
        with pytest.raises(ValueError, match='feature\\(s\\) left'):
            pipeline.filter_features(table, min_prevalence=1)


class TestMethodsText:

    def test_correlation_sparcc_and_no_correction(self):
        text = pipeline.correlation_methods('sparcc', 30, [('count table', [], 'all kept')], [], 200, 3, None, 0.05)
        assert 'SparCC (Friedman and Alm, 2012)' in text
        assert 'multiplicative replacement' not in text
        assert '200 permutations (seed 3)' in text
        assert 'not corrected' in text

    def test_association_options(self):
        text = pipeline.association_methods(['disease'], ['age', 'bmi'], 50, 2, ('count table', [], 'all kept'),
                                            [], 'bonferroni', 0.1, standardize=True, ancombc=True)
        assert '2 further sample(s) with missing metadata' in text
        assert 'CLR ~ disease + age + bmi' in text
        assert 'age, bmi were included as covariates' in text
        assert 'unit standard deviation' in text
        assert 'Bonferroni method, and results with q < 0.1' in text
        assert 'ANCOM-BC' in text
