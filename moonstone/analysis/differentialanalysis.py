import pandas as pd
import numpy as np
import scipy.stats as st

from moonstone.parsers.Filtering.concat_meta_and_readcounts import ConcatMetaAndReadCounts


class DifferentialAnalysis:

    def __init__(self, metadata_dataframe, reads_dataframe):
        self.read_count_df = reads_dataframe
        self.metadata_df = metadata_dataframe

    @property
    def full_table(self):
        if getattr(self, "_full_table", None) is None:
            instance = ConcatMetaAndReadCounts(self.metadata_df, self.read_count_df)
            setattr(self, "_full_table", instance.full_dataframe_with_features_in_columns)
        return self._full_table

    @property
    def number_columns_to_skip(self):
        if getattr(self, "_number_columns_to_skip", None) is None:
            setattr(self, "_number_columns_to_skip", len(self.metadata_df))
        return self._number_columns_to_skip

    def t_test(self, dichotomic_features, significance_level):
        features = []
        taxons = []
        static_value = []
        pvalue = []
        variance_group1 = []
        variance_group2 = []

        for feature in dichotomic_features:
            cat1 = self.full_table[self.full_table[feature] == self.full_table[feature][0]]
            cat2 = self.full_table[self.full_table[feature] != self.full_table[feature][0]]
            for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
                t_test = st.ttest_ind(cat1[self.full_table.columns[family]], cat2[self.full_table.columns[family]])
                if t_test[1] < significance_level:
                    features.append(feature)
                    taxons.append(self.full_table.columns[family])
                    static_value.append(t_test[0])
                    pvalue.append(round(t_test[1], 6))
                    variance_group1.append(cat1[self.full_table.columns[family]].var())
                    variance_group2.append(cat2[self.full_table.columns[family]].var())
        significant_differences_t_test = pd.DataFrame(list(zip(features, taxons, static_value, pvalue, variance_group1,
                                                      variance_group2)), columns=['features', 'taxons', 'static_value',
                                                      'p-value', 'variance_group1', 'variance_group2'])
        return significant_differences_t_test

    def wilcoxon_rank_test(self, dichotomic_features, significance_level):
        features = []
        taxons = []
        static_value = []
        pvalue = []
        variance_group1 = []
        variance_group2 = []

        for feature in dichotomic_features:
            cat1 = self.full_table[self.full_table[feature] == self.full_table[feature][0]]
            cat2 = self.full_table[self.full_table[feature] != self.full_table[feature][0]]
            for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
                ranksums_test = st.ranksums(cat1[self.full_table.columns[family]],
                                            cat2[self.full_table.columns[family]])
                if ranksums_test[1] < significance_level:
                    features.append(feature)
                    taxons.append(self.full_table.columns[family])
                    static_value.append(ranksums_test[0])
                    pvalue.append(round(ranksums_test[1], 6))
                    variance_group1.append(cat1[self.full_table.columns[family]].var())
                    variance_group2.append(cat2[self.full_table.columns[family]].var())
        significant_differences_ranksums = pd.DataFrame(list(zip(features, taxons, static_value, pvalue,
                                                        variance_group1, variance_group2)),
                                                        columns=['features', 'taxons', 'static_value', 'p-value',
                                                                 'variance_group1', 'variance_group2'])
        return significant_differences_ranksums

    def one_way_anova(self, multiple_option_features, significance_level):
        features = []
        taxons = []
        static_values = []
        pvalues = []

        for characteristic in multiple_option_features:
            variable_dic = {}
            for variable in self.full_table[characteristic].unique():
                variable_dic[variable] = self.full_table[self.full_table[characteristic] == variable]
            for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
                list_ofgroups = []
                for variable in variable_dic:
                    list_ofgroups.append(variable_dic[variable][self.full_table.columns[family]])
                oneway_anova = st.f_oneway(*np.asarray(list_ofgroups))
                if oneway_anova[1] < significance_level:
                    features.append(characteristic)
                    taxons.append(self.full_table.columns[family])
                    static_values.append(oneway_anova[0])
                    pvalues.append(round(oneway_anova[1], 6))

        signigicant_differences_oneway_anova = pd.DataFrame(list(zip(features, taxons, static_values, pvalues)),
                                                            columns=['features', 'taxons', 'static_value', 'p-value'])
        return signigicant_differences_oneway_anova

    def kruskal_test(self, multiple_option_features, significance_level):
        features = []
        taxons = []
        static_values = []
        pvalues = []

        for characteristic in multiple_option_features:
            variable_dic = {}
            for variable in self.full_table[characteristic].unique():
                variable_dic[variable] = self.full_table[self.full_table[characteristic] == variable]
            for family in range(self.number_columns_to_skip, self.full_table.shape[1]):
                list_ofgroups = []
                for variable in variable_dic:
                    list_ofgroups.append(variable_dic[variable][self.full_table.columns[family]])
                oneway_anova = st.kruskal(*np.asarray(list_ofgroups))
                if oneway_anova[1] < significance_level:
                    features.append(characteristic)
                    taxons.append(self.full_table.columns[family])
                    static_values.append(oneway_anova[0])
                    pvalues.append(oneway_anova[1])

        signigicant_differences_kruskal = pd.DataFrame(list(zip(features, taxons, static_values, pvalues)),
                                                       columns=['features', 'taxons', 'static_value', 'p-value'])
        return signigicant_differences_kruskal
