import pandas as pd
import typing
from typing import Optional, Union, List
import warnings

from moonstone.plot.plot_template import (
    BarGraph,
    Histogram
)

"""
plot that could be used on any dataframe no matter the module
Or should we put that in the analysis module?
"""


# What people might want to visualize?


def _check_list_of(listx, expectedtype):
    """
    check the type of what's in listx
    """
    print("called with ", listx, expectedtype)
    if all(isinstance(s, expectedtype) for s in listx) is True:
        return True
    else:
        return expectedtype._name+'['+expectedtype.__args__[0].__name__+']'


def _check_type(to_check, expectedtype: Union[type, typing._GenericAlias]):
    """
    check if type of 'to_check' is right
    """
    print("espectedtype =", expectedtype)
    if type(expectedtype) == list:      # if list of str/int etc. List[type]
        return _check_list_of(to_check, expectedtype)
    else:
        if type(to_check) == expectedtype:
            return True
        else:
            print(expectedtype)
            return expectedtype.__name__     # PROBLEM BC List[str]


def _check_type_s(to_check, expectedtype_s: Union[type, typing._GenericAlias]):
    """
    check if type of 'to_check' is right
    """
    if type(expectedtype_s) == list:    # if more than one type allowed
        type_in_str_for_warning = []
        for i in expectedtype_s:
            s = _check_type(to_check, i)
            if s is True:
                return True
            else:
                type_in_str_for_warning += s
        return " or ".join(type_in_str_for_warning)   # return what will be in the error raised
    else:                               # if only one type
        if _check_type(to_check, expectedtype_s):
            return True
        else:
            return expectedtype_s.__name__                         # return what will be in the error raised


def _check_types_in_plotting_options(plotting_options: dict):
    """
    check types of plotting options given by the user
    """
    expectedtype = {'log': bool, 'colorbar': [str, List[str]], 'tickangle': [int, float]}
    cleaned_plotting_options = {}
    for i in plotting_options.keys():
        print("hey ho", i, plotting_options[i], expectedtype[i])
        if i in expectedtype.keys():
            print("est-ce que tu m'entends")
            s = _check_type_s(plotting_options[i], expectedtype[i])
            if s is True:
                print("there")
                cleaned_plotting_options[i] = plotting_options[i]
            else:
                print("here")
                warnings.warn('Warning : %s value in plotting_options must be a %s, %s given. Value overidden' %
                              (i, s, type(plotting_options[i]).__name__))
    return cleaned_plotting_options


def _add_x_to_plotting_options(plotting_options: dict, x: str, defaultvalue):
    """
    don't override given plotting_options, meaning it only add the default value
    if value not already defined in plotting_options
    ~~~ MAYBE DELETE THIS METHOD ~~~
    """
    if x not in plotting_options.keys():    # if x not already specified (not given or not of the right type)
        plotting_options[x] = defaultvalue
        return plotting_options
    else:
        return plotting_options


class PlotStatsData():

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def plot_mean(self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options:
        """
        if plotting_options is None:
            plotting_options = {}
        # else:
        #    plotting_options = _check_types_in_plotting_options(plotting_options)

        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'log', True)
        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'tickangle', -60)

        df_mean = self.df.mean(axis=1)
        bar_fig = BarGraph(df_mean, plotting_options, show=show, output_file=output_file)

        bar_fig.compute_asymetric_bins()
        bar_fig.in_bins_and_count()    # normalize or not?
        bar_fig.plot_one_graph(
            "Distribution of %s mean" % self.items_name,
            "mean of the number of reads",
            "number of samples",
            )

    def plot_taxonomy_classification(self, level_of_interest, plotting_options: dict = None,
                                     show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the fractions of taxonomic levels (species/genus/...)
        :param level_of_interest: taxonomic classification level (species, genus etc.) to plot
        :param output_file: name of the output file
        :param plotting_options:
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        # mean by sample or nb of samples with ?

    # heatmap -> in analysis?


class PlotStatsMetadata():

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    def plot_sex(self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False):
        if plotting_options is None:
            plotting_options = {}
        # else:
        #    plotting_options = _check_types_in_plotting_options(plotting_options)

        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'colorbar', ['pink', 'blue'])

        bar_fig = BarGraph(self.metadata_df['sex'], plotting_options, show=show, output_file=output_file)
        bar_fig.count()    # normalize or not?
        bar_fig.reset_xnames({'F': 'Female', 'M': 'Male'})
        bar_fig.plot_one_graph(
            "Sex distribution in the samples",
            "sex",
            "number of samples",
            )

    def plot_age(self, step=None, plotting_options: dict = None,
                 show: Optional[bool] = True, output_file: Optional[str] = False):
        if plotting_options is None:
            plotting_options = {}
        # else:
        #    plotting_options = _check_types_in_plotting_options(plotting_options)

        hist_fig = Histogram(self.metadata_df['age'], plotting_options, show=show, output_file=output_file)

        if step is None:
            step = 1

        hist_fig.plot_one_graph(
            "Age distribution in the samples",
            "age",
            "number of samples",
            step
            )

    def plot_other(self, column_name, title=None, xlabel=None,
                   reset_xnames_dic: dict = None, plotting_options: dict = None,
                   show: Optional[bool] = True, output_file: Optional[str] = False):
        if plotting_options is None:
            plotting_options = {}
        # else:
        #    plotting_options = _check_types_in_plotting_options(plotting_options)

        bar_fig = BarGraph(self.metadata_df[column_name], plotting_options, show=show, output_file=output_file)
        bar_fig.count()    # normalize or not?
        if reset_xnames_dic is not None:
            bar_fig.reset_xnames(reset_xnames_dic)
        if title is None:
            title = column_name+" distribution in the samples"
        if xlabel is None:
            xlabel = column_name
        bar_fig.plot_one_graph(
            title,
            xlabel,
            "number of samples",
            )
