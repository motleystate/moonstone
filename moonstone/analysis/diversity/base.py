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
from moonstone.utils.dict_operations import merge_dict


logger = logging.getLogger(__name__)


class DiversityBase(BaseModule, BaseDF, ABC):
    _DIVERSITY_INDEXES_NAME = "index"
    _DEF_TITLE = "(index diversity) distribution across the samples"
    _AVAILABLE_GROUP_VIZ = ['violin', 'boxplot']
    _DEF_GROUP_VIZ = "boxplot"
    _DEF_PVAL_COLORSCALE = [
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
            self._diversity_indexes.name = self._DIVERSITY_INDEXES_NAME
        return self._diversity_indexes

    def _get_default_title(self) -> str:
        return f"{self.index_name} {self._DEF_TITLE}"

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
        if mode not in self._AVAILABLE_GROUP_VIZ:
            logger.warning("%s not a available mode, set to default (histogram)", mode)
            mode = self._DEF_GROUP_VIZ

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
            self._DIVERSITY_INDEXES_NAME, group_col,
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
            if isinstance(datastruct, pd.Series):
                datastruct = pd.concat([datastruct, datastruct.reorder_levels([1, 0])])
            else:  # ed. pd.DataFrame
                datastruct = datastruct.fillna(datastruct.transpose())
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
                    df[self._DIVERSITY_INDEXES_NAME], df[group_col], stats_test,
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
                # corrected_pval = corrected_pval.append(pval[pval.isnull()])
                corrected_pval = pd.concat([corrected_pval, pval[pval.isnull()]])
            # remodelling of p-values output
            corrected_pval = self._structure_remodelling(corrected_pval, structure=structure_pval, sym=sym)
            return corrected_pval
        else:
            pval = statistical_test_groups_comparison(
                df[self._DIVERSITY_INDEXES_NAME], df[group_col], stats_test,
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
            colorscale=self._DEF_PVAL_COLORSCALE,
            plotting_options=plotting_options,
            output_file=output_pval_file
        )

    def _valid_pval_param(self, pval_to_compute, pval_to_display):
        choices = [
            "all", "same group_col or group_col2 values", "same group_col values", None
        ]
        dicpval = {}
        for i in range(len(choices)):
            dicpval[choices[i]]=i
            
        if pval_to_compute not in choices:
            logger.warning("pval_to_compute='%s' not valid, set to default (all).", pval_to_compute)
            pval_to_compute = "all"

        if pval_to_display not in choices:
            logger.warning("pval_to_display='%s' not valid, set to default (None).", pval_to_display)
            pval_to_display = None
        elif dicpval[pval_to_display]<dicpval[pval_to_compute]:
            raise ValueError("pval_to_display='{}' not valid, when pval_to_compute='{}'. \
pval_to_display should be set to :{}".format(
                pval_to_display, pval_to_compute, choices[dicpval[pval_to_compute]:]
            ))
        
        return pval_to_compute, pval_to_display

    def _valid_correction_method_param(self, correction_method):
        if correction_method == "uncorrected":
            return None
        if correction_method not in [None, 'fdr_bh', 'bonferroni']:
            logger.warning("correction_method='%s' not valid, set to default (None).", correction_method)
            return None
        return correction_method

    def _pval_selection(
        self, pval_series, groups
    ):
        
        pval_series = pval_series[pval_series < 0.05]
        if groups is not None:
            pval_series = pval_series[
                (pval_series.index.get_level_values(0).isin(groups) & pval_series.index.get_level_values(1).isin(groups))
            ]
        return pval_series     
    
    def _pval_selection_with_group_col2(
        self, pval_series, final_groups, pval_to_compute, pval_to_display
    ):
        #Reminder: 
        #    1) This method called only if pval_to_display is not None. So pval_to_compute/pval_to_display = {"all", "same group_col or group_col2 values", "same group_col values"}
        #    2) Index values follow this pattern: "{group_col value} - {group_col2 value}"
        print("pval_to_compute:", pval_to_compute)
        print("pval_to_display:", pval_to_display)
        pval_series = self._pval_selection(pval_series, final_groups)
        if (pval_to_compute != "same group_col values" and
            pval_to_display == "same group_col values"):
            # we only have to check first part of index values
            pval_series = pval_series[
                (
                    pval_series.index.get_level_values(0).map(lambda x: x.split(" - ")[0]) == pval_series.index.get_level_values(1).map(lambda x: x.split(" - ")[0])
                )
            ]
        elif (pval_to_compute == "all" and
              pval_to_display == "same group_col or group_col2 values"):
            # we compare both part of the index values -> if first part is the same = same group_col value
                                                     # -> if second part is the same = same group_col2 value
            pval_series = pval_series[
                (
                    pval_series.index.get_level_values(0).map(lambda x: x.split(" - ")[0]) == pval_series.index.get_level_values(1).map(lambda x: x.split(" - ")[0])
                ) | (
                    pval_series.index.get_level_values(0).map(lambda x: x.split(" - ")[1]) == pval_series.index.get_level_values(1).map(lambda x: x.split(" - ")[1])    
                )]

        return pval_series

    def _order_pval_series(
        self, pval_series, groups, dic_gps
    ):
        # Reminder: pvalue series is a MultiIndex series -> 2 level of index are the 2 groups compared
        # to order p-value series in a specific order dictated in dic_gps
        # example: dic_gps = {"Group A": 0, "Group B": 1, "Group C": 3}
        names = pval_series.index.names     # default: ["Group1", "Group2"]

        # first we order p-value series index name to have the group that should come first as first member
        pval_series = pval_series.reset_index()
        for i in pval_series.index:
            level0 = pval_series.loc[i][names[0]]
            level1 = pval_series.loc[i][names[1]]
            if dic_gps[level0] > dic_gps[level1]:
                pval_series.loc[i, names[0]] = level1       # invert to have the Group that should be put first as first
                pval_series.loc[i, names[1]] = level0       # in example: "Group B - Group A" becomes "Group A - Group B"
        pval_series[names[0]] = pval_series[names[0]].astype("category")
        pval_series[names[0]].cat.set_categories(groups, inplace=True)
        pval_series[names[1]] = pval_series[names[1]].astype("category")
        pval_series[names[1]].cat.set_categories(groups, inplace=True)
        pval_series = pval_series.sort_values([names[0], names[1]]) # we sort by 1st member, and then 2nd member
        pval_series = pval_series.set_index([names[0], names[1]])
        return pval_series[0]

    def _generate_shapes_annotations_lists(
        self, pval_series, groups, hgt_min
    ):
        """
        To generate annotations to represent significant pvalues. Methods for group_col only (not group_col2)

        Args:
            pval_series: pd.Series of the pvalue to put on the graph (need to be filtered beforehand so that it only contains significant pvalues)
        """
        # Overview of this method: We're trying to generate the shapes (1 bracket = 3 lines =
        # 1 going from Group1 to Group2 and 2 small lines to form the edge of the bracket)
        # and the annotations.
        # To ensure that the brackets don't overlap, we create a table/list named `level`
        # with a number of columns equal to the number of groups and an expendable number of rows
        # When we fill the list, we replace 0 to 1 meaning that there is now a bracket at this
        # location and that we can't add another bracket on top of it, and that another level (row)
        # need to be added to `level`
        dic_gps = {}
        for i in range(len(groups)):
            dic_gps[groups[i]]=i        
        
        # the pvalues need to be ordered so that looking at the cell corresponding to the left edge 
        # of the bracket that needs to be added is enough to determine if there is already a
        # bracket there.
        pval_series = self._order_pval_series(pval_series, groups, dic_gps)
        
        det = (hgt_min/15)
        
        hgt_min += det/2
        fontsize = int(12+det)
        linewidth = 0.5+0.15*det
    
        level=[[0] * (len(groups))]
        list_shapes = []
        list_annotations = []
        
        for ind, val in pval_series.items():
            y = 0
            
            left_ind = dic_gps[ind[0]]  # -> could be directly ind for shapes but /!\ not for the annotations
            right_ind = dic_gps[ind[1]]
            
            for i in range(len(level)):
                if level[i][left_ind] == 0:
                    level[i][left_ind:right_ind+1]=[1]*(right_ind-left_ind+1)
                    list_shapes += [
                        {'x0':left_ind, 'y0':hgt_min+(i*det/2), 
                         'x1':right_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)},   # from Group1 to Group2
                        {'x0':left_ind, 'y0':hgt_min+(i*det/2)-0.15*det, 
                         'x1':left_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)},    # left edge of the bracket
                        {'x0':right_ind, 'y0':hgt_min+(i*det/2)-0.15*det, 
                         'x1':right_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)}    # right edge of the bracket
                    ]
                    if val < 0.01:
                        list_annotations += [{'text':'**', "font":dict(size=fontsize),
                                              'x':(left_ind+right_ind)/2, 
                                              'y':hgt_min+(i*det/2)+0.15*det, 'showarrow':False}]
                    else:
                        list_annotations += [{'text':'*', "font":dict(size=fontsize),
                                              'x':(left_ind+right_ind)/2, 
                                              'y':hgt_min+(i*det/2)+0.15*det, 'showarrow':False}]
                    y = 1
                    break
            if y == 0:
                # we need to add another level
                i += 1
                level += [[0] * (len(groups))]
                level[i][left_ind:right_ind+1]=[1]*(right_ind-left_ind+1)  # or ind
                list_shapes += [
                    {'x0':left_ind, 'y0':hgt_min+(i*det/2), 
                     'x1':right_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)},
                    {'x0':left_ind, 'y0':hgt_min+(i*det/2)-0.15*det, 
                     'x1':left_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)},
                    {'x0':right_ind, 'y0':hgt_min+(i*det/2)-0.15*det, 
                     'x1':right_ind, 'y1':hgt_min+(i*det/2), 'line':dict(width=linewidth)}]
                if val < 0.01:
                    list_annotations += [{'text':'**', "font":dict(size=fontsize),
                                          'x':(left_ind+right_ind)/2, 
                                          'y':hgt_min+(i*det/2)+0.15*det, 'showarrow':False}]
                else:
                    list_annotations += [{'text':'*',"font":dict(size=fontsize),
                                          'x':(left_ind+right_ind)/2,
                                          'y':hgt_min+(i*det/2)+0.15*det, 'showarrow':False}]
        
        return list_shapes, list_annotations, len(level)

    """
    def _generate_shapes_annotations_lists_supergroup(
        self, pval_series, groups, supergroups, hgt_min
    ):
    """
        #:param supergroups dictionary with supergroup edges
    """
        dic_gps = {}
        for i in range(len(groups)):
            dic_gps[groups[i]]=i
    
    
        list_shapes = []
        dic_middle = {}
        supergroups_to_display = set(list(pval_series.index.get_level_values(0)) + list(pval_series.index.get_level_values(1)))
        for i in supergroups_to_display:
            if isinstance(supergroups[i], list):
                list_shapes += [{'x0':supergroups[i][0], 'y0':hgt_min, 'x1':supergroups[i][1], 'y1':hgt_min, 'line':dict(width=1)}]                      
                dic_middle[i] = (dic_gps[supergroups[i][0]] + dic_gps[supergroups[i][1]])/2
            else:
                dic_middle[i] = dic_gps[supergroups[i]]

        dic_supergps = {}
        i = 0
        for k, v in sorted(dic_middle.items(), key=lambda item: item[1]):
            dic_supergps[k] = i
            i+=1
                
        level=[[0] * (len(dic_supergps))]
        list_annotations = []
        for ind, val in pval_series.items():
            y = 0
            for i in range(len(level)):
                if dic_middle[ind[0]] > dic_middle[ind[1]]:
                    ind = (ind[1], ind[0])
                if level[i][dic_supergps[ind[0]]] == 0:
                    level[i][dic_supergps[ind[0]]:dic_supergps[ind[1]]+1]=[1]*(dic_supergps[ind[1]]-dic_supergps[ind[0]]+1)  # or ind
                    list_shapes += [
                            {'x0':dic_middle[ind[0]], 'y0':hgt_min+i+0.5, 
                             'x1':dic_middle[ind[1]], 'y1':hgt_min+i+0.5, 'line':dict(width=1)},
                            {'x0':dic_middle[ind[0]], 'x1':dic_middle[ind[0]], 'y1':hgt_min+i+0.5,
                             #'y0':hgt_min, 'line':dict(width=1, dash="dot")},
                             'y0':hgt_min+i+0.35, 'line':dict(width=1)},
                            {'x0':dic_middle[ind[1]], 'x1':dic_middle[ind[1]], 'y1':hgt_min+i+0.5, 
                             #'y0':hgt_min, 'line':dict(width=1, dash="dot")}
                             'y0':hgt_min+i+0.35, 'line':dict(width=1)}
                                   ]
                    if val < 0.01:
                        list_annotations += [{'text':'**','x':(dic_middle[ind[0]]+dic_middle[ind[1]])/2, 
                                              'y':hgt_min+i+0.65, 'showarrow':False}]
                    else:
                        list_annotations += [{'text':'*','x':(dic_middle[ind[0]]+dic_middle[ind[1]])/2, 
                                              'y':hgt_min+i+0.65, 'showarrow':False}]
                    y = 1
                    break
            if y == 0:
                i += 1
                level += [[0] * (len(dic_supergps))]
                level[i][dic_supergps[ind[0]]:dic_supergps[ind[1]]+1]=[1]*(dic_supergps[ind[1]]-dic_supergps[ind[0]]+1)  # or ind
                list_shapes += [
                            {'x0':dic_middle[ind[0]], 'y0':hgt_min+i+0.5, 
                             'x1':dic_middle[ind[1]], 'y1':hgt_min+i+0.5, 'line':dict(width=1)},
                            {'x0':dic_middle[ind[0]], 'x1':dic_middle[ind[0]], 'y1':hgt_min+i+0.5, 
                             #'y0':hgt_min, 'line':dict(width=1, dash="dot")},
                             'y0':hgt_min+i+0.35, 'line':dict(width=1)},
                            {'x0':dic_middle[ind[1]], 'x1':dic_middle[ind[1]], 'y1':hgt_min+i+0.5, 
                             #'y0':hgt_min, 'line':dict(width=1, dash="dot")}
                             'y0':hgt_min+i+0.35, 'line':dict(width=1)}
                                ]
                if val < 0.01:
                    list_annotations += [{'text':'**','x':(dic_middle[ind[0]]+dic_middle[ind[1]])/2, 
                                              'y':hgt_min+i+0.65, 'showarrow':False}]
                else:
                    list_annotations += [{'text':'*','x':(dic_middle[ind[0]]+dic_middle[ind[1]])/2, 
                                              'y':hgt_min+i+0.65, 'showarrow':False}]
        
        return list_shapes, list_annotations, len(level)
    """
    # for now, method above and below in different methods
    #def _generate_dic_shapes_and_annotations(
    #    self, pval_series, dic_lev, hgt_min
    #):   
    #    list_shapes={}
    #    dic_annotations={}
    #    for ind, val in pval_series.items():

    def _compute_pval_inside_subgroups(
        self, diversity_index_dataframe: pd.DataFrame, group_col: str, final_group_col: str,
        stats_test: str, correction_method: str, structure_pval: str, sym: bool
    ):
        pval = pd.Series([], dtype='float64')
        for g in diversity_index_dataframe[group_col].dropna().unique():
            df_gp = diversity_index_dataframe[diversity_index_dataframe[group_col] == g]
            if df_gp.shape[0] < 2:
                logger.warning(
                    f"Less than 2 samples in dataframe group {g} in data. P-val can't be computed."
                )
            else:
                pval = pd.concat([
                    pval,
                    self._run_statistical_test_groups(
                        df_gp, final_group_col, stats_test,
                        correction_method, structure_pval, sym
                    )
                ])
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
        pval_to_compute: bool = 'all', pval_to_display: str = None,
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
        :param show_pval: visualize p-values's heatmap
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
        :param pval_to_display: whether you want the significant pvalues displayed on the graph ("all") or not (None)
          When group_col2 is used you may want to specify which type of comparisons you want to display the
          significant pvalues of. Otherwise the graph can appear crowded by pvalues lines.
          {None (default), "all", "same group_col values", "same group_col or group_col2 values"}
        """
        filtered_metadata_df = self._get_filtered_df_from_metadata(metadata_df)

        pval_to_compute, pval_to_display = self._valid_pval_param(pval_to_compute, pval_to_display)
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

            #if pval_to_display:
            if 1 == 1:
                # listing and sorting all final_groups possibles respecting the order given by
                # groups first and then by groups2
                t = df.drop_duplicates(subset=[final_group_col])
                t[group_col] = t[group_col].astype("category")
                t[group_col].cat.set_categories(groups, inplace=True)
                t[group_col2] = t[group_col2].astype("category")
                t[group_col2].cat.set_categories(groups2, inplace=True)
                t.dropna(how="any", subset=[group_col, group_col2], inplace=True)
                t = t.sort_values([group_col, group_col2])
                final_groups = list(t[final_group_col])

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
                    pval = pd.concat([
                        pval,
                        self._compute_pval_inside_subgroups(
                            df, group_col2, final_group_col,
                            stats_test, correction_method, structure_pval, sym
                        )
                    ])

            if pval_to_display:
                to_display = self._pval_selection_with_group_col2(pval, final_groups, pval_to_compute, pval_to_display)

        else:
            df = self._get_grouped_df(filtered_metadata_df[group_col])
            pval = self._run_statistical_test_groups(
                df, group_col, stats_test, correction_method, structure_pval, sym
                )
            # pval is in the right structure to be returned

            final_groups = groups   # to remove from here later on
            
            if pval_to_display:
                final_groups = groups
                to_display = self._pval_selection(pval, groups)
                
        if pval_to_display and to_display.empty:       # nothing to display
            pval_to_display = None

        if pval_to_display:
            hgt_min = df[self.DIVERSITY_INDEXES_NAME].max() #+ 0.5
            list_shapes, list_annotations, nb_lev = self._generate_shapes_annotations_lists(
                to_display, final_groups, hgt_min
            )
            if not plotting_options:
                plotting_options = {}
            det = (hgt_min/15)
            hgt_min += det/2
            nblevels=to_display.shape[0]
            plotting_options = merge_dict(plotting_options, {
                'layout': {
                    'shapes': list_shapes,            # we should had to previous list if there is one
                    'annotations': list_annotations,  # idem
                }
            })

        self.last_grouped_df = df
        self.report_data['analyse_groups'] = {
            'pval': pval,
            'param': {
                'pval_to_compute': pval_to_compute,
                'stats_test': stats_test,
                'correction_method': correction_method,
                'group_col': group_col,
                'group_col2': group_col2
            }
        }

        if make_graph:
            fig = self._make_graph(
                df, mode, group_col, group_col2, plotting_options, log_scale, show, output_file,
                colors, groups, groups2, **kwargs
            )
            if show_pval:
                if structure_pval != 'dataframe' or not sym:
                    pval_for_visualization = self._structure_remodelling(pval, 'dataframe', sym=True)
                    self._visualize_pvalue_matrix(pval_for_visualization, output_pval_file)
                else:
                    self._visualize_pvalue_matrix(pval, output_pval_file)

            return {**{'data': df, 'fig': fig}, **self.report_data['analyse_groups']}

        # 'data' different from 'diversity indexes' in the fact that it has been filtered on metadata, meaning that
        # samples without metadata for group_col (or group_col2) have been dropped
        return {**{'data': df}, **self.report_data['analyse_groups']}

    def generate_report_data(self) -> dict:
        """
        method that generates a report summurazing the diversity computed
        (parameters, results)
        """
        return {"title": self.index_name+" diversity", "diversity indexes": self.diversity_indexes}


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
              that your input data are valid. See
              `skbio.diversity <https://http://scikit-bio.org/docs/latest/diversity.html#module-skbio.diversity/>`_
              for the description of what validation entails so you can determine if you can safely disable validation."
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
