import logging
import re
from abc import ABC, abstractmethod
from string import capwords
from typing import Union

import pandas as pd

from moonstone.core.module_base import BaseModule, BaseDF
from moonstone.filtering.basics_filtering import NamesFiltering
from moonstone.plot.graphs.box import GroupBoxGraph, BoxGraph
from moonstone.plot.graphs.histogram import Histogram
from moonstone.plot.graphs.violin import GroupViolinGraph, ViolinGraph

logger = logging.getLogger(__name__)


class DiversityBase(BaseModule, BaseDF, ABC):
    DIVERSITY_INDEXES_NAME = "index"
    DEF_TITLE = "(index diversity) distribution across the samples"

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

    def analyse_groups(
        self, metadata_df: pd.DataFrame, group_col: str,  mode: str = 'boxplot',
        log_scale: bool = False, colors: dict = None, groups: list = None,
        show: bool = True, output_file: str = False, make_graph: bool = True,
        plotting_options: dict = None, **kwargs
    ) -> pd.DataFrame:
        """
        :param metadata_df: dataframe containing metadata and information to group the data
        :param group_col: column from metadata_df used to group the data
        :param mode: how to display (boxplot, or violin)
        :param colors: overides color for groups. format {group_id: color}
        :param groups: specifically select groups to display among group_col
        :param show: also visualize
        :param output_file: file path to output your html graph
        :param make_graph: whether or not to make the graph
        :param plotting_options: plotly plotting_options
        """
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)
        df = self._get_grouped_df(filtered_metadata_df[group_col])

        if make_graph:
            if mode not in ['violin', 'boxplot']:
                logger.warning("%s not a available mode, set to default (histogram)", mode)
                mode = "boxplot"

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

        self.last_grouped_df = df
        return df
