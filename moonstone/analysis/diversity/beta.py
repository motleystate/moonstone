import logging
import re
from abc import ABC, abstractmethod

import pandas as pd
import plotly.graph_objects as go
import skbio.diversity
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa

from moonstone.analysis.diversity.base import DiversityBase

logger = logging.getLogger(__name__)


class BetaDiversity(DiversityBase, ABC):
    DIVERSITY_INDEXES_NAME = "beta_index"
    DEF_TITLE = "(beta diversity) distribution across the samples"

    def __init__(self, dataframe: pd.DataFrame):
        super().__init__(dataframe)
        self.index_name = " ".join(re.findall('[A-Z][^A-Z]*', self.__class__.__name__)).capitalize()

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

    def _get_grouped_df(self, metadata_series):
        df_list = []
        for group in metadata_series.unique():
            group_df = self.df.loc[:, metadata_series[metadata_series == group].index]
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

    @property
    def pcoa(self):
        if getattr(self, '_pcoa', None) is None:
            self._pcoa = pcoa(self.beta_diversity).samples
        return self._pcoa

    def visualize_pcoa(
        self, metadata_df: pd.DataFrame, group_col: str,
        mode: str = 'scatter', **kwargs
    ):
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)

        if mode not in ['scatter']:
            logger.warning("%s not a available mode, set to default (scatter)", mode)
            mode = "scatter"
        title = "PCOA"
        xlabel = "PC1"
        ylabel = "PC2"
        plotting_options = self._handle_plotting_options(
            {}, title, xlabel, ylabel, False
        )

        groups = list(filtered_metadata_df[group_col].unique())
        fig = go.Figure()
        for group in groups:
            df = self.pcoa.loc[filtered_metadata_df[filtered_metadata_df[group_col] == group].index,:]
            fig.add_trace(go.Scatter(
                x=df['PC1'],
                y=df['PC2'],
                text=df.index,
                mode='markers',
                name=str(group)
            ))

            fig.update_layout(
                title={'text': title},
                xaxis={'title': 'PC1'},
                yaxis={'title': 'PC2'}
            )
        fig.show()


class BrayCurtis(BetaDiversity):
    """
    Perform calculation of the Bray Curtis for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df):    # compute_shannon_diversity
        """
        :param base: logarithm base chosen (NB : for ln, base=math.exp(1))
        """
        # steps to compute the index
        return skbio.diversity.beta_diversity("braycurtis", df.transpose(), df.columns)
