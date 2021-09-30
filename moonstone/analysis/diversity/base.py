import logging
import re
from abc import ABC, abstractmethod
import skbio
from string import capwords
from typing import Union

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

    def _get_grouped_df(self, metadata_series):
        return pd.concat([metadata_series, self.diversity_indexes], axis=1).dropna()

    def _get_filtered_df_from_metadata(self, metadata_df):
        return NamesFiltering(metadata_df, list(self.df.columns)).filtered_df

    def _make_graph(
        self, df, mode: str, group_col: str, plotting_options: dict, log_scale: bool,
        show: bool, output_file: str, colors: dict, groups: list, **kwargs
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
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            colors=colors,
            groups=groups,
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

    def analyse_groups(
        self, metadata_df: pd.DataFrame, group_col: str,  mode: str = 'boxplot',
        log_scale: bool = False, colors: dict = None, groups: list = None,
        show: bool = True, output_file: str = False, make_graph: bool = True,
        plotting_options: dict = None,
        stats_test: str = 'mann_whitney_u', correction_method: str = None,
        structure_pval: str = 'dataframe', sym: bool = True,
        show_pval: bool = True, output_pval_file: str = False,
        **kwargs
    ) -> dict:
        """
        :param metadata_df: dataframe containing metadata and information to group the data
        :param group_col: column from metadata_df used to group the data
        :param mode: how to display (boxplot, or violin)
        :param colors: overides color for groups. format {group_id: color}
        :param groups: specifically select groups to display among group_col
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
        """
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)
        df = self._get_grouped_df(filtered_metadata_df[group_col])

        pval = self._run_statistical_test_groups(df, group_col, stats_test, correction_method, structure_pval, sym)
        # pval is in the right structure to be returned

        if make_graph:
            self._make_graph(
                df, mode, group_col, plotting_options, log_scale, show, output_file,
                colors, groups, **kwargs
            )
            if show_pval:
                if structure_pval != 'dataframe' or not sym:
                    pval_for_visualization = self._structure_remodelling(pval, 'dataframe', sym=True)
                    self._visualize_pvalue_matrix(pval_for_visualization, output_pval_file)
                else:
                    self._visualize_pvalue_matrix(pval, output_pval_file)

        self.last_grouped_df = df
        return {
            'data': df,
            'pval': pval,
            'meta': {
                'stats_test': stats_test
            }
        }


class PhylogeneticDiversityBase(DiversityBase):
    """
    Class for Phylogenetic Diversities that use a taxonomy tree to compute the diversity.
    """
    def __init__(
        self,
        taxonomy_dataframe: pd.DataFrame,
        taxonomy_tree: skbio.TreeNode
    ):
        super().__init__(taxonomy_dataframe)
        if type(taxonomy_tree) == skbio.tree._tree.TreeNode:
            self.tree = taxonomy_tree
        else:
            raise RuntimeError("taxonomy_tree should be a skbio.TreeNode.")

    def _verification_otu_ids_in_tree(self, otu_ids):
        missing_ids = []
        for otu_id in otu_ids:
            try:
                self.tree.find(otu_id)
            except skbio.tree._exception.MissingNodeError:
                missing_ids += [otu_id]
        return missing_ids
