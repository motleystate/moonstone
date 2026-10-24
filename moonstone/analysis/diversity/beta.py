import logging
from abc import ABC, abstractmethod
from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import skbio.diversity
from skbio.stats.distance import DistanceMatrix, permanova, permdisp
from skbio.stats.ordination import pcoa
from statsmodels.stats.multitest import multipletests

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

    @staticmethod
    def _pcoa_biplot_features(ordination, features_df) -> pd.DataFrame:
        """
        Same computation as skbio.stats.ordination.pcoa_biplot (Legendre & Legendre 1998, chapter 9),
        returning the projected features only.

        skbio's version writes into a read-only view of the coordinates under pandas 3 (copy-on-write)
        and raises "output array is read-only", so the computation is done here on writable numpy copies.
        """
        coordinates = ordination.samples
        eigvals = np.array(ordination.eigvals, dtype=np.float64)
        n_samples = coordinates.shape[0]
        y = features_df.reindex(coordinates.index)

        coords = np.array(coordinates, dtype=np.float64)
        std = coords.std(axis=0, ddof=1)
        std[std == 0] = 1.0
        scaled_coords = (coords - coords.mean(axis=0)) / std

        spc = (1 / (n_samples - 1)) * np.asarray(y, dtype=np.float64).T.dot(scaled_coords)
        inv_sqrt_eigvals = np.zeros_like(eigvals)
        np.power(eigvals, -0.5, out=inv_sqrt_eigvals, where=eigvals > 0)
        uproj = np.sqrt(n_samples - 1) * spc.dot(np.diag(inv_sqrt_eigvals))
        return pd.DataFrame(uproj, index=y.columns.copy(), columns=coordinates.columns.copy()).fillna(0.0)

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
            d_bf = self._pcoa_biplot_features(self.pcoa, self.df.transpose())  # DF of PCoA feature contributions
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

    def _permanova_grouping(self, metadata_df: pd.DataFrame, group_col: str, min_group_size: int) -> tuple:
        """
        Return the group of each sample usable for PERMANOVA: samples present in both the counts and the metadata,
        with a value in group_col, and in a group of at least min_group_size samples. Empty or whitespace-only
        values count as missing.

        Returns:
            The grouping (a Series of strings indexed by sample) and the number of samples left out because their
            value is missing.
        """
        not_in_metadata = [s for s in self.df.columns if s not in metadata_df.index]
        if not_in_metadata:
            logger.warning(
                "%s: %i sample(s) have no metadata and are left out: %s",
                group_col, len(not_in_metadata), ", ".join(map(str, not_in_metadata))
            )
        grouping = metadata_df.loc[[s for s in self.df.columns if s in metadata_df.index], group_col]
        blank = grouping.notna() & grouping.astype(str).str.strip().eq('')
        grouping = grouping.mask(blank)
        n_missing = int(grouping.isna().sum())
        if n_missing:
            logger.warning("%s: %i sample(s) with a missing value are left out.", group_col, n_missing)
        # as strings before counting, so unused categories of a categorical column are not reported as groups
        grouping = grouping.dropna().astype(str)
        group_sizes = grouping.value_counts()
        too_small = group_sizes[group_sizes < min_group_size].index
        if len(too_small) > 0:
            logger.info(
                "%s: group(s) with fewer than %i samples are left out: %s",
                group_col, min_group_size, ", ".join(map(str, too_small))
            )
            grouping = grouping[~grouping.isin(too_small)]
        return grouping.rename(group_col), n_missing

    def permanova(
        self, metadata_df: pd.DataFrame, group_col: str,
        permutations: int = 999, seed: int = 0, min_group_size: int = 2,
        dispersion: bool = True, dispersion_test: str = 'median',
    ) -> pd.Series:
        """
        Test whether the groups of group_col differ in composition, with a Permutational Multivariate Analysis
        of Variance (PERMANOVA, Anderson 2001), on this object's distance matrix.

        PERMANOVA cannot tell a difference in composition (group centres) from a difference in spread (one group
        more variable than another). With dispersion=True, a Permutational Analysis of Multivariate Dispersions
        (PERMDISP, Anderson 2006; R vegan::betadisper) is run on the same samples. A small PERMDISP p-value means
        the groups differ in spread, so a significant PERMANOVA may reflect spread rather than composition.

        Samples without metadata or with a missing value in group_col (NaN, None, pd.NA, or an empty or
        whitespace-only string) are left out and counted, as are groups with fewer than min_group_size samples.
        PERMANOVA and PERMDISP use the same samples.

        Args:
            metadata_df: dataframe of metadata, indexed by sample.
            group_col: column of metadata_df defining the groups (treated as categories).
            permutations: number of label permutations used for the p-values.
            seed: random seed for the permutations, so results are reproducible. None for an unseeded run.
            min_group_size: smallest group kept in the test.
            dispersion: also run PERMDISP.
            dispersion_test: centre used by PERMDISP: 'median' (spatial median, as vegan::betadisper's default)
              or 'centroid'.

        Returns:
            Series named after group_col with: n_samples, n_missing (samples left out for a missing value),
            n_groups, pseudo_F, R2 (share of the variation explained by the groups), p_value, permdisp_F,
            permdisp_p_value (NaN when dispersion=False), and permutations.

        Note:
            PERMDISP works on Principal Coordinates Analysis (PCoA) axes. All axes with a positive eigenvalue are
            used. Non-Euclidean distances such as Bray-Curtis also give negative eigenvalues; scikit-bio leaves
            those axes out, whereas vegan::betadisper includes them, so PERMDISP values can differ slightly from R.
        """
        grouping, n_missing = self._permanova_grouping(metadata_df, group_col, min_group_size)
        return self._permanova_on_grouping(
            grouping, n_missing, min_group_size, permutations, seed, dispersion, dispersion_test
        )

    def _permanova_on_grouping(
        self, grouping: pd.Series, n_missing: int, min_group_size: int, permutations: int, seed: int,
        dispersion: bool, dispersion_test: str,
    ) -> pd.Series:
        group_col = grouping.name
        n_groups = grouping.nunique()
        if n_groups < 2:
            raise ValueError(
                f"{group_col}: fewer than 2 groups with at least {min_group_size} samples, PERMANOVA not possible."
            )
        if n_groups == grouping.shape[0]:
            raise ValueError(f"{group_col}: every sample is in its own group, PERMANOVA not possible.")
        distmat = self.beta_diversity.filter(grouping.index)
        grouping_df = grouping.to_frame()

        result = permanova(distmat, grouping_df, column=group_col, permutations=permutations, seed=seed)
        n_samples = int(result['sample size'])
        pseudo_f = float(result['test statistic'])
        # R2 = SS_between / SS_total, from pseudo-F = (SS_between / (g - 1)) / (SS_within / (n - g))
        r2 = pseudo_f * (n_groups - 1) / (pseudo_f * (n_groups - 1) + (n_samples - n_groups))

        permdisp_f = permdisp_p = np.nan
        if dispersion:
            disp = permdisp(
                distmat, grouping_df, column=group_col, test=dispersion_test, permutations=permutations,
                dimensions=distmat.shape[0], seed=seed, warn_neg_eigval=False,
            )
            permdisp_f = float(disp['test statistic'])
            permdisp_p = float(disp['p-value'])

        return pd.Series(
            {
                'n_samples': n_samples,
                'n_missing': n_missing,
                'n_groups': n_groups,
                'pseudo_F': pseudo_f,
                'R2': r2,
                'p_value': float(result['p-value']),
                'permdisp_F': permdisp_f,
                'permdisp_p_value': permdisp_p,
                'permutations': permutations,
            },
            name=group_col, dtype=object,
        )

    def permanova_screen(
        self, metadata_df: pd.DataFrame, columns: list = None, correction_method: str = 'fdr_bh',
        max_groups: int = 15, permutations: int = 999, seed: int = 0, min_group_size: int = 2,
        dispersion: bool = True, dispersion_test: str = 'median',
    ) -> pd.DataFrame:
        """
        Run permanova() for each column of metadata_df and correct the p-values for multiple testing.

        Testing many columns at p < 0.05 gives about one false positive per 20 columns with no real effect, so
        the corrected p-values (p_value_corrected, permdisp_p_value_corrected) are the ones to read.

        Args:
            metadata_df: dataframe of metadata, indexed by sample.
            columns: columns to test. Default: all columns of metadata_df.
            correction_method: 'fdr_bh' (Benjamini-Hochberg false discovery rate, default), 'bonferroni',
              or None for no correction (no corrected columns are added).
            max_groups: columns with more groups than this are skipped (usually continuous variables, which
              PERMANOVA on categories does not handle).
            permutations, seed, min_group_size, dispersion, dispersion_test: passed to permanova().

        Returns:
            Dataframe with one row per column: the fields of permanova(), the corrected p-values, and 'skipped'
            (why a column was not tested, empty otherwise). Correction is over the tested columns only.
        """
        correction_method = self._valid_correction_method_param(correction_method)
        if columns is None:
            columns = list(metadata_df.columns)
        rows = []
        for col in columns:
            grouping, n_missing = self._permanova_grouping(metadata_df, col, min_group_size)
            values = metadata_df.loc[metadata_df.index.isin(self.df.columns), col]
            n_values = values.mask(values.astype(str).str.strip().eq('')).nunique()
            n_groups = grouping.nunique()
            skipped = ''
            if n_values > max_groups:
                skipped = f'more than {max_groups} groups'
            elif n_groups < 2:
                skipped = f'fewer than 2 groups with at least {min_group_size} samples'
            elif n_groups == grouping.shape[0]:
                skipped = 'every sample in its own group'
            if skipped:
                logger.warning("%s: not tested, %s.", col, skipped)
                rows.append(pd.Series({'n_samples': grouping.shape[0], 'n_missing': n_missing, 'n_groups': n_groups,
                                       'skipped': skipped}, name=col))
                continue
            row = self._permanova_on_grouping(
                grouping, n_missing, min_group_size, permutations, seed, dispersion, dispersion_test
            )
            row['skipped'] = ''
            rows.append(row)

        columns_order = [
            'n_samples', 'n_missing', 'n_groups', 'pseudo_F', 'R2', 'p_value', 'permdisp_F', 'permdisp_p_value',
            'permutations', 'skipped'
        ]
        screen = pd.DataFrame(rows).reindex(columns=columns_order)
        for col in columns_order[:-1]:
            screen[col] = pd.to_numeric(screen[col])

        if correction_method is not None:
            for pval_col in ['p_value', 'permdisp_p_value']:
                corrected = pd.Series(np.nan, index=screen.index)
                tested = screen[pval_col].dropna()
                if not tested.empty:
                    corrected[tested.index] = multipletests(tested, method=correction_method)[1]
                screen.insert(screen.columns.get_loc(pval_col) + 1, f'{pval_col}_corrected', corrected)
        return screen


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

        return skbio.diversity.beta_diversity(
            "weighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
            validate=self.validate, taxa=otu_ids, tree=self.tree,
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
        return skbio.diversity.beta_diversity(
            "unweighted_unifrac", counts=df.loc[otu_ids].transpose(), ids=df.columns,
            validate=self.validate, taxa=otu_ids, tree=self.tree,
        )
