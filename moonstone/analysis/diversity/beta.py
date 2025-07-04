import logging
from abc import ABC, abstractmethod
from typing import Union

import numpy as np
from packaging import version
import pandas as pd
import plotly.graph_objects as go
import skbio.diversity
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import pcoa, pcoa_biplot

from moonstone.analysis.diversity.base import (
    DiversityBase, PhylogeneticDiversityBase
)
from moonstone.plot.graphs.scatter import GroupScatterGraph, GroupScatter3DGraph
from moonstone.utils.log_messages import warn_once, reset_warnings_decorator

logger = logging.getLogger(__name__)


class BetaDiversity(DiversityBase, ABC):
    _DIVERSITY_INDEXES_NAME = "beta_index"
    _DEF_TITLE = "(beta diversity) distribution across the samples"

    @abstractmethod
    def compute_beta_diversity(self, df) -> DistanceMatrix:
        """
        method that compute the beta diversity
        """
        pass

    def compute_diversity(self) -> pd.Series:
        series = self.beta_diversity.to_series()
        series.name = self._DIVERSITY_INDEXES_NAME
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
            beta_div_solo_indexed_df.columns = [self._DIVERSITY_INDEXES_NAME]
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
            beta_div_solo_indexed_df.columns = [self._DIVERSITY_INDEXES_NAME]
            # Add corresponding group name
            beta_div_solo_indexed_df[final_group_col] = group
            beta_div_solo_indexed_df[group_col] = group.split(" - ")[0]
            beta_div_solo_indexed_df[group_col2] = group.split(" - ")[1]
            df_list.append(beta_div_solo_indexed_df)
        return pd.concat(df_list).dropna()

    def _get_grouped_df(self, metadata_df):
        if isinstance(metadata_df, pd.Series):
            return self._get_grouped_df_series(metadata_df)
        else:
            return self._get_grouped_df_dataframe(metadata_df)

    @property
    def pcoa(self):
        if getattr(self, '_pcoa', None) is None:
            self._pcoa = pcoa(self.beta_diversity)
        return self._pcoa

    def _label_with_proportion(self, var):
        string = "{} ({:2.2%})".format(var, self.pcoa.proportion_explained[var])
        return string

    def _unit_scale(self, x, y, scale=1):
        """
        Scales to a range of 1 to -1.
        Each value for a given component, x, is divided by the range of values of that same component, y.
        """
        z = x/y * scale
        return z

    def _scale_biplot(self, pcoa_biplot_features_df, pcoa_sample_coordinates_df) -> pd.DataFrame:
        """
        skbio PCoA defaults to a scale of 1 (e.g. -1 to 0, 0 to 1, -0.5 to 0.5, -0.3 to 0.7, etc...)
        This functions first calls the UNIT_SCALE function. It then centers the features to match the
        PCoA figure scaling so as to generate an accurate overlay of features and samples.

        Calculations and scaling is done individually for each axis.
        """
        d_pbf = pcoa_biplot_features_df

        for axis in d_pbf.columns:  # for each principle coordinate (PC)
            # Get the range of coordinates for the samples
            pcoa_range = abs(pcoa_sample_coordinates_df[axis].min()) + abs(pcoa_sample_coordinates_df[axis].max())

            # Feature have a range of 1.0, and so we scale this to sample range
            d_pbf[axis] = self._unit_scale(
                d_pbf[axis],
                d_pbf[axis].max() - d_pbf[axis].min(),
                scale=pcoa_range)

            # Center features with respect to samples so ploting fisuals make sense
            positive_feature_correction = pcoa_sample_coordinates_df[axis].max() - d_pbf[axis].max()
            d_pbf[axis] = d_pbf[axis] + positive_feature_correction
        return d_pbf

    def _map_loadings(self, pci, n):
        pci_pos = list(np.argsort(pci)[-n:])  # Sorts values but returns the argument position!
        pci_neg = list(np.argsort(pci)[:n])   # Allows alignment with feature name.
        return pci_pos, pci_neg

    def _bi_plotly(self, fig, pcx, pcy, features, n):
        """
        Args:
            loadings: is the vector for each indiviudal variable (columns) in as many dimentions as PCs (rows).
            n: Number of individual loadings to plot.
        """
        pcx_pos, pcx_neg = self._map_loadings(pcx, n)
        pcy_pos, pcy_neg = self._map_loadings(pcy, n)
        map_loadings = set(pcx_pos + pcx_neg + pcy_pos + pcy_neg)

        for m in map_loadings:
            fig.add_annotation(x=pcx[m], y=pcy[m], showarrow=True, arrowhead=2, arrowwidth=2,
                               axref='x', ayref='y', ax=0, ay=0)
            if type(features[m]) is str:
                formatted_text = f"{features[m]}"
            else:
                formatted_text = f"{features[m][-1]}"  # Handle multi-index parts
            fig.add_annotation(x=pcx[m]*1.1, y=pcy[m]*1.1, text=formatted_text,
                               ay=0, font={"color": "black", "size": 14})
        return fig

    def _bi_plotly3d(self, fig, pcx, pcy, pcz, features, n):
        """
        Args:
            loadings: is the vector for each indiviudal variable (columns) in as many dimentions as PCs (rows).
            n: Number of individual loadings to plot.
        """
        pcx_pos, pcx_neg = self._map_loadings(pcx, n)
        pcy_pos, pcy_neg = self._map_loadings(pcy, n)
        pcz_pos, pcz_neg = self._map_loadings(pcz, n)
        map_loadings = set(pcx_pos + pcx_neg + pcy_pos + pcy_neg + pcz_pos + pcz_neg)

        annotations = []
        for m in map_loadings:
            if type(features[m]) is str:
                formatted_text = f"{features[m]}"
            else:
                formatted_text = f"{features[m][-1]}"  # Handle multi-index parts

            fig.add_trace(go.Scatter3d(
                x=[0, pcx[m]],  # Start at (0, 0, 0) and go to (x, y, z)
                y=[0, pcy[m]],
                z=[0, pcz[m]],
                mode='lines', line=dict(color='black', width=2),
                showlegend=False
            ))

            # Add annotation near the arrowhead
            fig.add_trace(go.Scatter3d(
                x=[pcx[m]*1.1], y=[pcy[m]*1.1], z=[pcz[m]*1.1],
                mode='text', text=[formatted_text], textposition='top center',
                showlegend=False
            ))

            fig.add_trace(go.Cone(
                x=[pcx[m]],
                y=[pcy[m]],
                z=[pcz[m]],
                u=[pcx[m]],
                v=[pcy[m]],
                w=[pcz[m]],
                sizeref=0.05,  # Adjust size of the cone
                showlegend=False,
                showscale=False,
                colorscale=[[0, 'rgb(0,0,0)'], [1, 'rgb(0,0,0)']]
                ))

        fig.update_layout(scene={"annotations": annotations})
        return fig

    @reset_warnings_decorator
    def visualize_pcoa(
        self, metadata_df: pd.DataFrame, group_col: str, group_col2: str = None,
        mode: str = 'scatter',
        proportions: bool = False, x_pc: int = 1, y_pc: int = 2, z_pc: int = 3,
        show: bool = True, output_file: Union[bool, str] = False,
        colors: dict = None, symbols: dict = None,
        groups: list = None, groups2: list = None,
        n_biplot_features: int = 0, plotting_options: dict = None,
        q_confel: float = 0,
    ):
        """
        Args:
            metadata_df: dataframe containing metadata and information to group the data.
            group_col: column of the metadata_df to use for coloring.
            group_col2: (optional) column from metadata_df to use for.
            mode: type of graph to visualize the PCoA. { 'scatter' (default), 'scatter3d' }.
            proportions: write proportion explained for each PC in the x/y labels.
            colors: dictionary to impose specific colors to one or more groups of group_col.
            symbols: dictionary to impose specific symbols to one or more groups of group_col2.
            groups: select specific groups to display among group_col.
            groups2: select specific groups to display among group_col2.
            n_biplot_features: add arrows showing the n most explanatory features per axis direction.
            q_confel: proportion (between 0 and 1) of confidence intervals to represent with ellipses on the graph.
              NB: 0 = no ellipse. Only available with mode="scatter" (2D).
        """
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)
        if group_col2:
            df = pd.concat([self.pcoa.samples, filtered_metadata_df[[group_col, group_col2]]], axis=1)
        else:
            df = pd.concat([self.pcoa.samples, filtered_metadata_df[group_col]], axis=1)

        if mode not in ['scatter', 'scatter3d']:
            logger.warning("'%s' not a available mode, set to default (scatter).", mode)
            mode = "scatter"

        title = f"PCOA of samples from {self.index_name} distance matrix"
        xvar = "PC"+str(x_pc)
        yvar = "PC"+str(y_pc)
        if proportions:
            xlabel = self._label_with_proportion(xvar)
            ylabel = self._label_with_proportion(yvar)
        else:
            xlabel = xvar
            ylabel = yvar

        if n_biplot_features > 0:
            tmp_show = False
            tmp_output_file = False
        else:
            tmp_show = show
            tmp_output_file = output_file

        if mode == "scatter":
            plotting_options = self._handle_plotting_options(
                plotting_options, title, xlabel, ylabel, False
            )
            graph = GroupScatterGraph(df)
            args_for_plot = [xvar, yvar, group_col]
        else:                       # mode == "scatter3d"
            zvar = "PC"+str(z_pc)
            if proportions:
                zlabel = self._label_with_proportion(zvar)
            else:
                zlabel = zvar
            plotting_options = self._handle_plotting_options(
                plotting_options, title, xlabel, ylabel, False, zlabel=zlabel
            )
            graph = GroupScatter3DGraph(df)
            args_for_plot = [xvar, yvar, zvar, group_col]
            if q_confel != 0:
                logger.info("Ellipses are not available for mode='scatter3d' (yet).")
                # @TODO: eventually try to compute ellipse of confidence for 3D plot
        fig = graph.plot_one_graph(
            *args_for_plot,
            group_col2=group_col2,
            plotting_options=plotting_options,
            show=tmp_show,
            output_file=tmp_output_file,
            colors=colors,
            symbols=symbols,
            groups=groups,
            groups2=groups2,
            q_confel=q_confel
        )

        if n_biplot_features > 0:
            my_pcoa_biplot = pcoa_biplot(self.pcoa, self.df.transpose())
            d_bf = my_pcoa_biplot.features                           # DF of PCoA feature contributions
            d_sbf = self._scale_biplot(d_bf, self.pcoa.samples)      # DF of SCALLED PCoA features contribution
            loadings = d_sbf.transpose().values
            pcx = loadings[x_pc-1, :]
            pcy = loadings[y_pc-1, :]
            if mode == "scatter3d":
                pcz = loadings[z_pc-1, :]
                fig = self._bi_plotly3d(fig, pcx, pcy, pcz, d_bf.index, n_biplot_features)
            else:
                fig = self._bi_plotly(fig, pcx, pcy, d_bf.index, n_biplot_features)

            graph._handle_output_plotly(fig, show, output_file)
        return fig


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


class Jaccard(BetaDiversity):
    """
    Perform calculation of the Jaccard distance for each pairs of samples from the dataframe
    """
    def compute_beta_diversity(self, df):
        # steps to compute the index
        return skbio.diversity.beta_diversity("jaccard", df.transpose(), df.columns)


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
                warn_once(logger, f"INCOMPLETE TREE: missing {missing_ids}.\n\
Computation of the Weighted UniFrac diversity using only the OTU IDs present in the Tree.")
                otu_ids = list(set(otu_ids) - set(missing_ids))

        skbio_version = version.parse(skbio.__version__)
        if skbio_version >= version.parse("0.6.0"):
            return skbio.diversity.beta_diversity(
                "weighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
                validate=self.validate, taxa=otu_ids, tree=self.tree,
                )
        else:
            return skbio.diversity.beta_diversity(
                "weighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
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
                warn_once(logger, f"INCOMPLETE TREE: missing {missing_ids}.\n\
Computation of the Unweighted UniFrac diversity using only the OTU IDs present in the Tree.")
                otu_ids = list(set(otu_ids) - set(missing_ids))

        # df = counts dataframe
        skbio_version = version.parse(skbio.__version__)
        if skbio_version >= version.parse("0.6.0"):
            return skbio.diversity.beta_diversity(
                "unweighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
                validate=self.validate, taxa=otu_ids, tree=self.tree,
                )
        else:
            return skbio.diversity.beta_diversity(
                "unweighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
                validate=self.validate, otu_ids=otu_ids, tree=self.tree,
                )
