import logging
from abc import ABC, abstractmethod
from typing import Union

import pandas as pd
import skbio.diversity
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa

from moonstone.analysis.diversity.base import (
    DiversityBase, PhylogeneticDiversityBase
)
from moonstone.plot.graphs.scatter import GroupScatterGraph, GroupScatter3DGraph

logger = logging.getLogger(__name__)


class BetaDiversity(DiversityBase, ABC):
    DIVERSITY_INDEXES_NAME = "beta_index"
    DEF_TITLE = "(beta diversity) distribution across the samples"

    @abstractmethod
    def compute_beta_diversity(self, df) -> DistanceMatrix:
        """
        method that compute the beta diversity
        """
        pass

    def compute_diversity(self) -> pd.Series:
        series = self.beta_diversity.to_series()
        series.name = self.DIVERSITY_INDEXES_NAME
        return series

    @property
    def beta_diversity(self):
        """
        DistanceMatrix from skbio.
        """
        if getattr(self, '_beta_diversity', None) is None:
            self._beta_diversity = self.compute_beta_diversity(self.df)
        return self._beta_diversity

    @property
    def beta_diversity_series(self):
        return self.diversity_indexes

    @property
    def beta_diversity_df(self):
        return self.beta_diversity.to_data_frame()

    def _get_grouped_df_series(self, metadata_series):
        df_list = []
        for group in metadata_series.dropna().unique():
            group_df = self.df.loc[:, metadata_series[metadata_series == group].index]
            # computing beta diversities only between samples from the same group
            beta_div_multi_indexed_df = self.compute_beta_diversity(group_df).to_series().to_frame()
            if beta_div_multi_indexed_df.empty:  # Happens if only one item from the group
                continue
            # Make unique index from multi index
            beta_div_not_indexed_df = beta_div_multi_indexed_df.reset_index()
            index_col_names = ["level_0", "level_1"]
            beta_div_solo_indexed_df = beta_div_not_indexed_df.set_index(
                beta_div_not_indexed_df[index_col_names].agg('-'.join, axis=1)
            ).drop(index_col_names, axis=1)
            beta_div_solo_indexed_df.columns = [self.DIVERSITY_INDEXES_NAME]
            # Add corresponding group name
            beta_div_solo_indexed_df[metadata_series.name] = group
            df_list.append(beta_div_solo_indexed_df)
        return pd.concat(df_list).dropna()

    def _get_grouped_df_dataframe(self, metadata_dataframe):
        df_list = []
        final_group_col = metadata_dataframe.columns[-1]
        group_col = metadata_dataframe.columns[0]
        group_col2 = metadata_dataframe.columns[1]
        for group in metadata_dataframe[final_group_col].dropna().unique():
            group_df = self.df.loc[:, metadata_dataframe[metadata_dataframe[final_group_col] == group].index]
            # computing beta diversities only between samples from the same group
            beta_div_multi_indexed_df = self.compute_beta_diversity(group_df).to_series().to_frame()
            if beta_div_multi_indexed_df.empty:  # Happens if only one item from the group
                continue
            # Make unique index from multi index
            beta_div_not_indexed_df = beta_div_multi_indexed_df.reset_index()
            index_col_names = ["level_0", "level_1"]
            beta_div_solo_indexed_df = beta_div_not_indexed_df.set_index(
                beta_div_not_indexed_df[index_col_names].agg('-'.join, axis=1)
            ).drop(index_col_names, axis=1)
            beta_div_solo_indexed_df.columns = [self.DIVERSITY_INDEXES_NAME]
            # Add corresponding group name
            beta_div_solo_indexed_df[final_group_col] = group
            beta_div_solo_indexed_df[group_col] = group.split(" - ")[0]
            beta_div_solo_indexed_df[group_col2] = group.split(" - ")[1]
            df_list.append(beta_div_solo_indexed_df)
        return pd.concat(df_list).dropna()

    def _get_grouped_df(self, metadata_df):
        if type(metadata_df) == pd.core.series.Series:
            return self._get_grouped_df_series(metadata_df)
        else:
            return self._get_grouped_df_dataframe(metadata_df)

    @property
    def pcoa(self):
        if getattr(self, '_pcoa', None) is None:
            self._pcoa = pcoa(self.beta_diversity).samples
        return self._pcoa

    def visualize_pcoa(
        self, metadata_df: pd.DataFrame, group_col: str, mode: str = 'scatter',
        show: bool = True, output_file: Union[bool, str] = False,
        colors: dict = None, groups: list = None,
        plotting_options: dict = None,
    ):

        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)
        df = pd.concat([self.pcoa, filtered_metadata_df[group_col]], axis=1)

        if mode not in ['scatter', 'scatter3d']:
            logger.warning("%s not a available mode, set to default (scatter)", mode)
            mode = "scatter"
        title = f"PCOA of samples from {self.index_name} distance matrix"
        xlabel = "PC1"
        ylabel = "PC2"
        plotting_options = self._handle_plotting_options(
            plotting_options, title, xlabel, ylabel, False
        )

        if mode == "scatter":
            graph = GroupScatterGraph(df)
            args_for_plot = ["PC1", "PC2", group_col]
        elif mode == "scatter3d":
            graph = GroupScatter3DGraph(df)
            args_for_plot = ["PC1", "PC2", "PC3", group_col]
        graph.plot_one_graph(
            *args_for_plot,
            plotting_options=plotting_options,
            show=show,
            output_file=output_file,
            colors=colors,
            groups=groups,
        )


class BrayCurtis(BetaDiversity):
    """
    Perform calculation of the Bray Curtis for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df):    # compute_bray_curtis_diversity
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        # steps to compute the index
        return skbio.diversity.beta_diversity("braycurtis", df.transpose(), df.columns)


class WeightedUniFrac(BetaDiversity, PhylogeneticDiversityBase):
    """
    Perform calculation of the weighted UniFrac for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df) -> skbio.stats.distance._base.DistanceMatrix:
        # steps to compute the index
        otu_ids = df.index

        missing_ids = self._verification_otu_ids_in_tree(otu_ids)
        if len(missing_ids) > 0:
            if not self.force_computation:
                raise RuntimeError(f"INCOMPLETE TREE: missing {missing_ids}.")
            else:
                logger.warning(f"INCOMPLETE TREE: missing {missing_ids}.\n\
Computation of the Weighted UniFrac diversity using only the OTU IDs present in the Tree.")
                otu_ids = list(set(otu_ids) - set(missing_ids))

        return skbio.diversity.beta_diversity(
            "weighted_unifrac", df.loc[otu_ids].transpose(), df.columns,
            validate=self.validate, otu_ids=otu_ids, tree=self.tree,
            )


class UnweightedUniFrac(BetaDiversity, PhylogeneticDiversityBase):
    """
    Perform calculation of the unweighted UniFrac for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df) -> skbio.stats.distance._base.DistanceMatrix:
        # steps to compute the index
        otu_ids = df.index

        missing_ids = self._verification_otu_ids_in_tree(otu_ids)
        if len(missing_ids) > 0:
            if not self.force_computation:
                raise RuntimeError(f"INCOMPLETE TREE: missing {missing_ids}.")
            else:
                logger.warning(f"INCOMPLETE TREE: missing {missing_ids}.\n\
Computation of the Unweighted UniFrac diversity using only the OTU IDs present in the Tree.")
                otu_ids = list(set(otu_ids) - set(missing_ids))

        return skbio.diversity.beta_diversity(
            "unweighted_unifrac", df.loc[otu_ids].transpose(), df.columns,
            validate=self.validate, otu_ids=otu_ids, tree=self.tree,
            )
