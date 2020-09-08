import pandas as pd
import sys
import typing
from typing import Optional, List
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
    return all(isinstance(s, expectedtype) for s in listx)


def _check_type(to_check, expectedtype: type):
    """
    check if type of 'to_check' is right
    """
    # type(List[str]) = typing._GenericAlias in python3.7 and typing.GenericMeta in python3.6
    if sys.version_info[1] >= 7:
        typeList = typing._GenericAlias
        # typeListName = expectedtype._name       # for now only typing type is List so not necessary
    else:
        typeList = typing.GenericMeta
        # typeListName = str(p._gorg).split('.')[1]
    if type(expectedtype) == typeList:      # if List[<type>] (<type> = str, int etc.)
        if type(to_check) == list:
            if _check_list_of(to_check, expectedtype.__args__[0]) is True:
                return True
        return 'List['+expectedtype.__args__[0].__name__+']'
    else:
        if type(to_check) == expectedtype:
            return True
        else:
            return expectedtype.__name__


def _check_type_s(to_check, expectedtype_s: type):
    """
    check if type of 'to_check' is right.
    Dispatcher if several types are allowed (then expectedtype_s is a list) or not
    """
    if type(expectedtype_s) == list:    # if more than one type allowed
        type_in_str_for_warning = []
        for i in expectedtype_s:
            s = _check_type(to_check, i)
            if s is True:
                return True
            else:
                type_in_str_for_warning += [s]
        return " or ".join(type_in_str_for_warning)   # return what will be in the error raised
    else:                               # if only one type
        if _check_type(to_check, expectedtype_s):
            return True
        else:
            return expectedtype_s.__name__                         # return what will be in the error raised


def _check_types_in_plotting_options(plotting_options: dict):
    """
    check types of plotting options given by the user.
    Create a new dictionary cleaned_plotting_options without the options that are of the wrong type
    """
    expectedtype = {'log': bool, 'colorbar': [str, List[str]], 'tickangle': [int, float]}
    cleaned_plotting_options = {}
    for i in plotting_options.keys():
        if i in expectedtype.keys():
            s = _check_type_s(plotting_options[i], expectedtype[i])
            if s is True:
                cleaned_plotting_options[i] = plotting_options[i]
            else:
                if type(plotting_options[i]) == list and 'List' in s:
                    # warning when the problem is not the overall type but the type of elements in list
                    warnings.warn(('Warning : %s value in plotting_options must be a %s. \
Please check the type of elements in the list given. Value overridden') % (i, s))
                else:
                    warnings.warn('Warning : %s value in plotting_options must be a %s, %s given. Value overridden' %
                                  (i, s, type(plotting_options[i]).__name__))
        else:
            cleaned_plotting_options[i] = plotting_options[i]
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


class PlotCountsStats():

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def plot_mean(self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the mean distribution of the number of reads by items

        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'log', True)
        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'tickangle', -60)

        df_mean = self.df.mean(axis=1)

        bar_fig = BarGraph(df_mean)
        bar_fig.compute_heterogeneous_bins()
        bar_fig.in_bins_and_count()    # normalize or not?
        bar_fig.plot_one_graph(
            "Distribution of %s mean" % self.items_name,
            "mean of the number of reads",
            "number of samples",
            plotting_options,
            show=show,
            output_file=output_file
            )


class PlotTaxonomyStats():

    def __init__(self, dataframe: pd.DataFrame, items_name: str = "items"):
        self.df = dataframe
        self.items_name = items_name

    def _plot_taxonomy_classification(self, level_of_interest, plotting_options: dict = None,
                                      show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        ~~~~ IN CONSTRUCTION (remove _ in title when finished) ~~~~
        method to visualize the fractions of taxonomic levels (species/genus/...)

        :param level_of_interest: taxonomic classification level (species, genus etc.) to plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        # mean by sample or nb of samples with ?

    # heatmap -> in analysis?


class PlotMetadataStats():

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    def plot_age(self, bins_size=None, plotting_options: dict = None,
                 show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the age distribution of patients (whose the samples are originated from)

        :param bins_size: size of the bins of the Histogram
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        hist_fig = Histogram(self.metadata_df['age'])

        if bins_size is None:
            bins_size = 1

        hist_fig.plot_one_graph(
            "Age distribution in the samples",
            "age",
            "number of samples",
            bins_size,
            plotting_options,
            show=show,
            output_file=output_file
            )

    def plot_sex(self, plotting_options: dict = None, show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        method to visualize the sex distribution of patients (whose the samples are originated from)

        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        plotting_options = _add_x_to_plotting_options(
            plotting_options, 'colorbar', ['pink', 'blue'])

        bar_fig = BarGraph(self.metadata_df['sex'])
        bar_fig.count()    # normalize or not?
        bar_fig.reset_xnames({'F': 'Female', 'M': 'Male'})
        bar_fig.plot_one_graph(
            "Sex distribution in the samples",
            "sex",
            "number of samples",
            plotting_options,
            show=show,
            output_file=output_file
            )

    def plot_other(self, column_name, title=None, xlabel=None,
                   reset_xnames_dic: dict = None, plotting_options: dict = None,
                   show: Optional[bool] = True, output_file: Optional[str] = False):
        """
        :param column_name: name of the column you wish to display into a barplot
        :param title: title of the graph
        :param xlabel: label of the axis
        :param reset_xnames_dic: to rename the names of the values in the x axis. \n
                                 Example for a plot of the distribution of smoking habits :
                                 reset_xnames_dic={'y': 'smoker', 'n': 'non smoker'}
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        if plotting_options is None:
            plotting_options = {}
        else:
            plotting_options = _check_types_in_plotting_options(plotting_options)

        bar_fig = BarGraph(self.metadata_df[column_name])
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
            plotting_options,
            show=show,
            output_file=output_file
            )
