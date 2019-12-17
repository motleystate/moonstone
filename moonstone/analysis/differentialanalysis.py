import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

from moonstone.filtering.merging_meta_and_readcounts import MergingMetaAndReadCounts


class DifferentialAnalysis:

    def __init__(self, metadata_dataframe, reads_dataframe):
        self.read_count_df = reads_dataframe
        self.metadata_df = metadata_dataframe
        print('removing rows with only zeros')
        self.read_count_df = self.read_count_df[self.read_count_df.sum(axis=1) != 0.0]
        instance = MergingMetaAndReadCounts(self.metadata_df, self.read_count_df)
        self.full_table = instance.full_dataframe_with_features_in_columns
        self.available_tests = {'dichotomic': ['t_test', 'wilcoxon_rank_test'],
                                'multiple': ['one_way_anova', 'kruskal_test']}

    @property
    def number_columns_to_skip(self):
        if getattr(self, "_number_columns_to_skip", None) is None:
            setattr(self, "_number_columns_to_skip", len(self.metadata_df))
        return self._number_columns_to_skip

    def test_t_test(self, feature):
        features = []
        taxons = []
        static_value = []
        pvalue = []
        variance_group1 = []
        variance_group2 = []
        cat1 = self.full_table[self.full_table[feature] == self.full_table[feature][0]]
        cat2 = self.full_table[self.full_table[feature] != self.full_table[feature][0]]
        for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
            t_test = st.ttest_ind(cat1[self.full_table.columns[family]], cat2[self.full_table.columns[family]])
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_value.append(round(t_test[0], 6))
            pvalue.append(round(t_test[1], 6))
            variance_group1.append(cat1[self.full_table.columns[family]].var())
            variance_group2.append(cat2[self.full_table.columns[family]].var())
        t_test_results = pd.DataFrame(list(zip(features, taxons, static_value, pvalue, variance_group1,
                                               variance_group2)), columns=['features', 'taxons', 'static_value',
                                                                           'p-value', 'variance_group1',
                                                                           'variance_group2'])
        return t_test_results

    def test_wilcoxon_rank_test(self, feature):
        features = []
        taxons = []
        static_value = []
        pvalue = []
        variance_group1 = []
        variance_group2 = []
        cat1 = self.full_table[self.full_table[feature] == self.full_table[feature][0]]
        cat2 = self.full_table[self.full_table[feature] != self.full_table[feature][0]]
        for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
            ranksums_test = st.ranksums(cat1[self.full_table.columns[family]],
                                        cat2[self.full_table.columns[family]])
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_value.append(round(ranksums_test[0], 6))
            pvalue.append(round(ranksums_test[1], 6))
            variance_group1.append(cat1[self.full_table.columns[family]].var())
            variance_group2.append(cat2[self.full_table.columns[family]].var())
            wilcoxon_ranksums_results = pd.DataFrame(list(zip(features, taxons, static_value, pvalue,
                                                              variance_group1, variance_group2)),
                                                     columns=['features', 'taxons', 'static_value', 'p-value',
                                                              'variance_group1', 'variance_group2'])
        return wilcoxon_ranksums_results

    def test_one_way_anova(self, feature):
        features = []
        taxons = []
        static_values = []
        pvalues = []
        variable_dic = {}
        for variable in self.full_table[feature].unique():
            variable_dic[variable] = self.full_table[self.full_table[feature] == variable]
        for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
            list_ofgroups = []
            for variable in variable_dic:
                list_ofgroups.append(variable_dic[variable][self.full_table.columns[family]])
            oneway_anova = st.f_oneway(*np.asarray(list_ofgroups))
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_values.append(round(oneway_anova[0], 6))
            pvalues.append(round(oneway_anova[1], 6))

        oneway_anova_results = pd.DataFrame(list(zip(features, taxons, static_values, pvalues)),
                                            columns=['features', 'taxons', 'static_value', 'p-value'])
        return oneway_anova_results

    def test_kruskal_test(self, feature):
        features = []
        taxons = []
        static_values = []
        pvalues = []
        variable_dic = {}
        for variable in self.full_table[feature].unique():
            variable_dic[variable] = self.full_table[self.full_table[feature] == variable]
        for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
            list_ofgroups = []
            for variable in variable_dic:
                list_ofgroups.append(variable_dic[variable][self.full_table.columns[family]])
            kruskal_test = st.kruskal(*np.asarray(list_ofgroups))
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_values.append(round(kruskal_test[0], 6))
            pvalues.append(kruskal_test[1])

        kruskal_test_results = pd.DataFrame(list(zip(features, taxons, static_values, pvalues)),
                                            columns=['features', 'taxons', 'static_value', 'p-value'])
        return kruskal_test_results

    def test_default(self, *args, **kwargs):
        raise Exception("For dichotomic features use {} and for multiple option features use {}."
                        .format(', '.join(self.available_tests['dichotomic']),
                                ', '.join(self.available_tests['multiple'])))

    def corrected_p_values(self, p_value_serie, correction_method_used):
        """
        Some available methods are:
        - 'bonferroni' : one-step correction
        - 'fdr_bh' : Benjamini/Hochberg (non-negative)
        """
        corrected_p_values = multipletests(p_value_serie, method=correction_method_used)
        return corrected_p_values[1]

    def differential_analysis_by_feature(self, test, features, correction_method_used):
        '''
        Features should be provided in a list splited by dichotomic or multiple option features.
        example:
        dicotomic_features = ['SEX', 'GRIPPE']
        multiple_option_features = ['Season', Age_Group']
        '''
        final_table = pd.DataFrame()
        for feature in features:
            test_result = getattr(self, f"test_{test}", self.test_default)(feature)
            test_result['corrected_p-value'] = self.corrected_p_values(test_result['p-value'], correction_method_used)
            final_table = final_table.append(test_result)
        return final_table
