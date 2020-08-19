import logging

import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

from moonstone.transformers.mergers import MergeCountsAndMetadata
from moonstone.filtering.basics_filtering import NoCountsFiltering

logger = logging.getLogger(__name__)


class DifferentialAnalysis:

    available_tests = {
        'dichotomic': ['t_test', 'wilcoxon_rank_test'],
        'multiple': ['one_way_anova', 'kruskal_test']
    }
    tests_functions_used = {
        't_test': st.ttest_ind,
        'wilcoxon_rank_test': st.ranksums,
        'one_way_anova': st.f_oneway,
        'kruskal_test': st.kruskal

    }

    def __init__(self, metadata_dataframe, read_count_dataframe):
        self.read_count_df = read_count_dataframe
        self.metadata_df = metadata_dataframe
        logger.info('removing rows with only zeros')
        filtering_instance = NoCountsFiltering(self.read_count_df)
        self.read_count_df = filtering_instance.filtered_df
        instance = MergeCountsAndMetadata(self.metadata_df, self.read_count_df)
        self.full_table = instance.full_dataframe_with_features_in_columns

    @property
    def number_columns_to_skip(self):
        if getattr(self, "_number_columns_to_skip", None) is None:
            setattr(self, "_number_columns_to_skip", len(self.metadata_df))
        return self._number_columns_to_skip

    def test_dichotomic_features(self, feature, test_to_use):
        features = []
        taxons = []
        static_value = []
        pvalue = []
        variance_group1 = []
        variance_group2 = []
        cat1 = self.full_table[self.full_table[feature] == self.full_table[feature][0]]
        cat2 = self.full_table[self.full_table[feature] != self.full_table[feature][0]]
        for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
            test = self.tests_functions_used[test_to_use](cat1[self.full_table.columns[family]],
                                                          cat2[self.full_table.columns[family]])
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_value.append(round(test[0], 6))
            pvalue.append(round(test[1], 6))
            variance_group1.append(cat1[self.full_table.columns[family]].var())
            variance_group2.append(cat2[self.full_table.columns[family]].var())
        test_results = pd.DataFrame(list(zip(features, taxons, static_value, pvalue, variance_group1,
                                             variance_group2)), columns=['features', 'taxons', 'static_value',
                                                                         'p-value', 'variance_group1',
                                                                         'variance_group2'])
        return test_results

    def test_multiple_features(self, feature, test_to_use):
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
            test = self.tests_functions_used[test_to_use](*np.asarray(list_ofgroups))
            features.append(feature)
            taxons.append(self.full_table.columns[family])
            static_values.append(round(test[0], 6))
            pvalues.append(round(test[1], 6))

        test_results = pd.DataFrame(list(zip(features, taxons, static_values, pvalues)),
                                    columns=['features', 'taxons', 'static_value', 'p-value'])
        return test_results

    def test_default(self, *args, **kwargs):
        raise Exception("For 'dichotomic_features' use {} and for 'multiple_features' use {}."
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

    def differential_analysis_by_feature(self, features, type_of_features, test_to_use, correction_method_used):
        '''
        Features should be provided in a list splited by dichotomic or multiple option features.
        example:
        dicotomic_features = ['SEX', 'GRIPPE']
        multiple_option_features = ['Season', Age_Group']
        '''
        final_table = pd.DataFrame()
        for feature in features:
            test_result = getattr(self, f"test_{type_of_features}", self.test_default)(feature, test_to_use)
            test_result['corrected_p-value'] = self.corrected_p_values(test_result['p-value'], correction_method_used)
            final_table = final_table.append(test_result)
        return final_table
