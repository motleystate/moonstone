import logging
import re
from abc import ABC, abstractmethod
import skbio
from string import capwords
from typing import Union

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from moonstone.analysis.statistical_test import statistical_test_groups_comparison
from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.filtering.basics_filtering import NamesFiltering
from moonstone.plot.graphs.box import GroupBoxGraph, BoxGraph
from moonstone.plot.graphs.heatmap import HeatmapGraph
from moonstone.plot.graphs.histogram import Histogram
from moonstone.plot.graphs.violin import GroupViolinGraph, ViolinGraph

logger = logging.getLogger(__name__)


class DiversityBase(BaseModule, BaseDF, ABC):
    DIVERSITY_INDEXES_NAME = "index"
    DEF_TITLE = "(index diversity) distribution across the samples"
    AVAILABLE_GROUP_VIZ = ['violin', 'boxplot']
    DEF_GROUP_VIZ = "boxplot"
    DEF_PVAL_COLORSCALE = [
        [0, '#FFFF00'], [0.001, '#f5962a'], [0.05, '#FF0000'], [0.050001, '#000055'], [1, '#000000']
    ]

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame]):
        """
        :param dataframe: taxonomy count table normalized to have each samples with the same number of reads
        """
        super().__init__(dataframe)
        self.index_name = capwords(" ".join(re.findall('[A-Z][^A-Z]*', self.__class__.__name__)))

    @abstractmethod
    def compute_diversity(self) -> pd.Series:
        """
        method that compute the diversity
        """
        pass

    @property
    def diversity_indexes(self):
        # call compute_alpha_diversity and store into self._alpha_indexes
        if getattr(self, '_diversity_indexes', None) is None:
            self._diversity_indexes = self.compute_diversity()
            self._diversity_indexes.name = self.DIVERSITY_INDEXES_NAME
        return self._diversity_indexes

    def _get_default_title(self) -> str:
        return f"{self.index_name} {self.DEF_TITLE}"

    def _get_default_samples_label(self) -> str:
        return "Samples"

    def _visualize_histogram(self, bins_size, plotting_options, show, output_file, log_scale: bool, **kwargs):
        title = self._get_default_title()
        xlabel = self.index_name
        ylabel = self._get_default_samples_label()
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        graph = Histogram(self.diversity_indexes)
        graph.plot_one_graph(
            bins_size,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            **kwargs
        )

    def _visualize_violin(self, plotting_options: dict, show: bool, output_file: str, log_scale: bool, **kwargs):
        title = self._get_default_title()
        xlabel = self._get_default_samples_label()
        ylabel = self.index_name
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        graph = ViolinGraph(self.diversity_indexes)
        graph.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            **kwargs
        )

    def _visualize_boxplot(self, plotting_options: dict, show: bool, output_file: str, log_scale: bool, **kwargs):
        title = self._get_default_title()
        xlabel = self._get_default_samples_label()
        ylabel = self.index_name
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        graph = BoxGraph(self.diversity_indexes)
        graph.plot_one_graph(
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            **kwargs
        )

    def visualize(
        self, mode: str = 'histogram', bins_size: Union[int, float] = 0.1, log_scale: bool = False,
        show: bool = True, output_file: str = False, plotting_options: dict = None, **kwargs
    ):
        """
        :param mode: how to display (histogram, boxplot, or violin)
        :param bins_size: [mode histo only] size of the histo bins
        :param show: display your graph
        :param output_file: file path to output your html graph
        :param plotting_options: plotly plotting_options
        """
        if mode not in ['histogram', 'violin', 'boxplot']:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = "histogram"

        if mode == "histogram":
            self._visualize_histogram(bins_size, plotting_options, show, output_file, log_scale, **kwargs)
        elif mode == "violin":
            self._visualize_violin(plotting_options, show, output_file, log_scale, **kwargs)
        elif mode == "boxplot":
            self._visualize_boxplot(plotting_options, show, output_file, log_scale, **kwargs)

    def _get_grouped_df(self, metadata_df):
        return pd.concat([metadata_df, self.diversity_indexes], axis=1).dropna()

    def _get_filtered_df_from_metadata(self, metadata_df):
        return NamesFiltering(metadata_df, list(self.df.columns)).filtered_df

    def _make_graph(
        self, df, mode: str, group_col: str, group_col2: str, plotting_options: dict, log_scale: bool,
        show: bool, output_file: str, colors: dict, groups: list, groups2: list, **kwargs
    ):
        if mode not in self.AVAILABLE_GROUP_VIZ:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = self.DEF_GROUP_VIZ

        title = self._get_default_title()
        xlabel = f"{group_col}"
        ylabel = f"{self.index_name.capitalize()}"
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, log_scale
        )
        if mode == "violin":
            graph = GroupViolinGraph(df)
        elif mode == "boxplot":
            graph = GroupBoxGraph(df)
        graph.plot_one_graph(
            self.DIVERSITY_INDEXES_NAME, group_col,
            group_col2=group_col2,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            colors=colors,
            groups=groups,
            groups2=groups2,
            **kwargs
        )

    def _structure_remodelling(self, datastruct: Union[pd.Series, pd.DataFrame], structure: str, sym: bool):
        if sym:
            datastruct = pd.concat([datastruct, datastruct.reorder_levels([1, 0])])
        if structure == 'dataframe':
            datastruct = datastruct.unstack(level=1)
            datastruct.index.name = None
            datastruct.columns.name = None
        return datastruct

    def _run_statistical_test_groups(
        self, df: pd.DataFrame, group_col: str, stats_test: str, correction_method: str,
        structure_pval: str, sym: bool
    ):
        if correction_method is not None:
            pval = statistical_test_groups_comparison(
                    df[self.DIVERSITY_INDEXES_NAME], df[group_col], stats_test,
                    output='series', sym=False
                )

            if pval.dropna().shape[0] == 0:
                logger.warning(
                    "Only NaN in the p-value dataframe: Meaning there are too few samples in all groups."
                )
                return np.nan

            # dropping NaN (= comparison that couldn't have been generated due to too few samples in one or both groups)
            corrected_pval = pd.Series(multipletests(pval.dropna(), alpha=0.05, method=correction_method)[1])

            corrected_pval.index = pval.dropna().index   # postulate that the order hasn't changed
            if pval[pval.isnull()].size > 0:
                corrected_pval = corrected_pval.append(pval[pval.isnull()])

            # remodelling of p-values output
            corrected_pval = self._structure_remodelling(corrected_pval, structure=structure_pval, sym=sym)
            return corrected_pval
        else:
            pval = statistical_test_groups_comparison(
                df[self.DIVERSITY_INDEXES_NAME], df[group_col], stats_test,
                output=structure_pval, sym=sym
            )
            return pval

    def _visualize_pvalue_matrix(self, pval: pd.DataFrame, output_pval_file: str):
        graph = HeatmapGraph(pval)
        plotting_options = {
            'layout': {
                'title': 'Heatmap visualization of p-values',
            }
        }
        graph.plot_one_graph(
            colorscale=self.DEF_PVAL_COLORSCALE,
            plotting_options=plotting_options,
            output_file=output_pval_file
        )

    def _valid_pval_param(self, pval_to_compute):
        choices = [
            "all", "same group_col or group_col2 values", "same group_col values", None
        ]
        if pval_to_compute not in choices:
            logger.warning("pval_to_compute='%s' not valid, set to default (all).", pval_to_compute)
            pval_to_compute = "all"
        return pval_to_compute

    def _valid_correction_method_param(self, correction_method):
        if correction_method == "uncorrected":
            return None
        if correction_method not in [None, 'fdr_bh', 'bonferroni']:
            logger.warning("correction_method='%s' not valid, set to default (None).", correction_method)
            return None
        return correction_method

    def _compute_pval_inside_subgroups(
        self, diversity_index_dataframe: pd.DataFrame, group_col: str, final_group_col: str,
        stats_test: str, correction_method: str, structure_pval: str, sym: bool
    ):
        pval = pd.Series([])
        for g in diversity_index_dataframe[group_col].dropna().unique():
            df_gp = diversity_index_dataframe[diversity_index_dataframe[group_col] == g]
            if df_gp.shape[0] < 2:
                logger.warning(
                    f"Less than 2 samples in dataframe group {g} in data. P-val can't be computed."
                )
            else:
                pval = pval.append(self._run_statistical_test_groups(
                    df_gp, final_group_col, stats_test,
                    correction_method, structure_pval, sym
                ))
        pval.index = pd.MultiIndex.from_tuples(pval.index, names=('Group1', 'Group2'))
        return pval

    def analyse_groups(
        self, metadata_df: pd.DataFrame, group_col: str, group_col2: str = None,
        mode: str = 'boxplot',
        log_scale: bool = False, colors: dict = None,
        groups: list = None, groups2: list = None,
        show: bool = True, output_file: str = False, make_graph: bool = True,
        plotting_options: dict = None,
        stats_test: str = 'mann_whitney_u', correction_method: str = None,
        structure_pval: str = 'dataframe', sym: bool = True,
        pval_to_compute: bool = 'all',
        show_pval: bool = True, output_pval_file: str = False,
        **kwargs
    ) -> dict:
        """
        :param metadata_df: dataframe containing metadata and information to group the data
        :param group_col: column from metadata_df used to group the data
        :param group_col2: (optional) column from metadata_df used to further divide the data
        :param mode: how to display (boxplot, or violin)
        :param colors: overides color for groups. format {group_id: color}
        :param groups: specifically select groups to display among group_col
        :param groups2: specifically select groups to display among group_col2
        :param show: also visualize
        :param show_pval: visualize p-values
        :param output_file: file path to output your html graph
        :param make_graph: whether or not to make the graph
        :param plotting_options: plotly plotting_options
        :param stats_test: {'mann_whitney_u', 'ttest_independence', 'chi2_contingency'} statistical test
        used to calculate the p-values between each groups
        :param correction_method: {None, 'fdr_bh' (benjamini-hochberg), 'bonferroni'} method used (if any)
        to correct generated p-values
        :param structure_pval: {'series', 'dataframe'}
        :param sym: whether generated dataframe (or MultiIndexed series) is symetric or half-full
        :param pval_to_compute: if group_col2 used, problems of memory or in maximum recursion depth
        may occur. In this case, you may want to compute only p-values of specific comparisons.
        {"all" (default), None, "same group_col values", "same group_col or group_col2 values"}
        """
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)

        pval_to_compute = self._valid_pval_param(pval_to_compute)
        correction_method = self._valid_correction_method_param(correction_method)

        if group_col2:
            final_group_col = group_col+"_"+group_col2
            filtered_metadata_df[final_group_col] = np.where(
                np.logical_or(
                    filtered_metadata_df[group_col].isnull() is True, filtered_metadata_df[group_col2].isnull() is True
                ), np.nan,
                filtered_metadata_df[group_col].astype(str) + " - " +
                filtered_metadata_df[group_col2].astype(str)
            )
            df = self._get_grouped_df(filtered_metadata_df[[group_col, group_col2, final_group_col]])

            if pval_to_compute == "all":
                pval = self._run_statistical_test_groups(
                    df, final_group_col, stats_test, correction_method, structure_pval, sym
                )
            elif (pval_to_compute == "same group_col values" or
                  pval_to_compute == "same group_col or group_col2 values"):
                pval = self._compute_pval_inside_subgroups(
                    df, group_col, final_group_col, stats_test, correction_method, structure_pval, sym
                )
                if pval_to_compute == "same group_col or group_col2 values":
                    pval = pval.append(
                        self._compute_pval_inside_subgroups(
                            df, group_col2, final_group_col,
                            stats_test, correction_method, structure_pval, sym
                        )
                    )

        else:
            df = self._get_grouped_df(filtered_metadata_df[group_col])
            pval = self._run_statistical_test_groups(
                df, group_col, stats_test, correction_method, structure_pval, sym
                )
            # pval is in the right structure to be returned

        self.last_grouped_df = df
        report_dictionary = {
            'data': df,
            'pval': pval,
            'meta': {
                'pval_to_compute': pval_to_compute,
                'stats_test': stats_test,
                'correction_method': correction_method,
            }
        }

        if make_graph:
            fig = self._make_graph(
                df, mode, group_col, group_col2, plotting_options, log_scale, show, output_file,
                colors, groups, groups2, **kwargs
            )
            report_dictionary['fig'] = fig
            if show_pval:
                if structure_pval != 'dataframe' or not sym:
                    pval_for_visualization = self._structure_remodelling(pval, 'dataframe', sym=True)
                    self._visualize_pvalue_matrix(pval_for_visualization, output_pval_file)
                else:
                    self._visualize_pvalue_matrix(pval, output_pval_file)

        return report_dictionary


class PhylogeneticDiversityBase(DiversityBase):
    """
    Class for Phylogenetic Diversities that use a taxonomy tree to compute the diversity.
    """
    def __init__(
        self,
        taxonomy_dataframe: pd.DataFrame,
        taxonomy_tree: skbio.TreeNode,
        validate: bool = True,
        force_computation: bool = False
    ):
        """
        Args:
            validate: skbio argument. "If False, validation of the input won’t be performed.
            This step can be slow, so if validation is run elsewhere it can be disabled here.
            However, invalid input data can lead to invalid results or error messages that
            are hard to interpret, so this step should not be bypassed if you’re not certain
            that your input data are valid. See skbio.diversity for the description of what
            validation entails so you can determine if you can safely disable validation.
            force_computation: if True, doesn't raise error if OTU IDs are missing and compute
            the diversity with the OTU IDs that are present in the Tree
        """
        super().__init__(taxonomy_dataframe)
        if type(taxonomy_tree) == skbio.tree._tree.TreeNode:
            self.tree = taxonomy_tree
        else:
            raise RuntimeError("taxonomy_tree should be a skbio.TreeNode.")
        self.force_computation = force_computation
        self.validate = validate

    def _verification_otu_ids_in_tree(self, otu_ids):
        missing_ids = []
        for otu_id in otu_ids:
            try:
                self.tree.find(otu_id)
            except skbio.tree._exception.MissingNodeError:
                missing_ids += [otu_id]
        return missing_ids
