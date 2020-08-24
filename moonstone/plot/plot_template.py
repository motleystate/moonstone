from abc import ABC, abstractmethod
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io
from typing import Union

from moonstone.plot.functions import (
    check_list_of
)


class BaseGraph(ABC):

    def __init__(self, dataframe: pd.DataFrame, plotting_options: dict = None,
                 show: bool = True, output_file: Union[bool, str] = False):
        """
        :param dataframe: pandas dataframe generated with moonstone
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        self.df = dataframe
        self.plotting_options = plotting_options

        if type(show) == bool:
            self.show = show
        else:
            raise ValueError('Error : show value must be a bool, %s given' % type(show).__name__)

        self.output_file = output_file
        if self.output_file is True:
            # if no name given for the output file, a generic name is generated
            if dataframe.name is not None:
                self.output_file = dataframe.name+"_"+self.__class__.__name__+".html"
            else:
                self.output_file = self.__class__.__name__+".html"

    @abstractmethod
    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass


class BarGraph(BaseGraph):

    def __init__(self, dataframe: pd.DataFrame, plotting_options: dict = None,
                 show: bool = True, output_file: Union[bool, str] = False):
        super().__init__(dataframe, plotting_options, show, output_file)

        # check type of data : numerical(int or float) or categorical(anything really)
        # -> np.nan considered as float so all good!

        if isinstance(dataframe, pd.Series):
            if dataframe.apply(np.isreal).all():
                self.datatype = 'numerical'
            else:
                self.datatype = 'categorical'
        # if self.dataframe.apply(np.isreal).all().all():   #only numerical values

    def compute_heterogeneous_bins(self):
        """Logish bins"""
        maximum = self.df.max()
        magnitude = int(math.log10(maximum))
        bval = [-0.1]                       # to have the 0 value
        i = 0
        while i < magnitude:
            bval += list(np.arange(2*10**i, 10**(i+1)+1, 10**i))
            i += 1
        # i=magnitude
        bval += list(np.arange(2*10**i, maximum+10**i, 10**i))        # up to maximum
        return bval

    @property
    def bins_values(self):
        """
        retrieves bins_values, and compute it if no values given
        """
        if getattr(self, '_bins_values', None) is None:
            self.bins_values = self.compute_heterogeneous_bins()
        return self._bins_values

    @bins_values.setter
    def bins_values(self, bins_values):
        if type(bins_values) == list and check_list_of(bins_values, (int, float, np.integer)):
            self._bins_values = bins_values
        else:
            raise ValueError("Error : expecting list of numerical values (int, float) in bins_values.")

    def in_bins_and_count(self, normalize: bool = False):
        binned_df = pd.cut(self.df, bins=self.bins_values)   # put every items in the appropriate bin
        counts = pd.value_counts(binned_df, normalize=normalize)
        counts = counts.reindex(binned_df.cat.categories)
        self.xnames = list(counts.index.astype(str))
        self.xnames[0] = self.xnames[0].replace('(-0.1', '[0.0')
        self.xnames = [self.xnames[i].replace('(', ']') for i in range(len(self.xnames))]
        self.yvalues = list(counts)

    def count(self, normalize: bool = False):
        counts = pd.value_counts(self.df, normalize=normalize)
        counts = counts.sort_index()
        self.xnames = list(counts.index)
        self.yvalues = list(counts)

    def reset_xnames(self, reset_xnames_dic: dict):
        for i in range(len(self.xnames)):
            self.xnames[i] = reset_xnames_dic[self.xnames[i]]

    def plot_one_graph(self, title: str = None, xlabel: str = None, ylabel: str = None,
                       reset_xnames_dic: dict = None, plotting_options: dict = None):
        """
        :param title: title of the graph
        :param xlabel: label of the x axis
        :param ylabel: label of the y axis
        :param reset_xnames_dic: to rename the names of the values in the x axis. \n
                                 Example for a plot of the distribution of smoking habits :
                                 reset_xnames_dic={'y': 'smoker', 'n': 'non smoker'}
        """
        if self.datatype == 'numerical':
            self.in_bins_and_count()
        else:
            self.count()

        if reset_xnames_dic is not None:
            self.reset_xnames(reset_xnames_dic)

        fig = go.Figure([go.Bar(x=self.xnames, y=self.yvalues)])

        if plotting_options is not None:
            if 'log' in self.plotting_options.keys() and self.plotting_options['log']:
                fig.update_layout(yaxis_type="log")
            if 'tickangle' in self.plotting_options.keys():
                fig.update_xaxes(tickangle=self.plotting_options['tickangle'])
            if 'colorbar' in self.plotting_options.keys():
                fig.update_traces(marker_color=self.plotting_options['colorbar'])

        if title is not None:
            fig.update_layout(title_text=title, title_x=0.5)
        if xlabel is not None:
            fig.update_xaxes(title_text=xlabel)
        if ylabel is not None:
            fig.update_yaxes(title_text=ylabel)

        if self.show is True:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)


class Histogram(BaseGraph):

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str,
                       step: Union[int, float]):   # see if we add nbinsx options
        minimum = int(self.df.min()/step) * step
        maximum = self.df.max()
        fig = go.Figure([go.Histogram(x=self.df, xbins=dict(start=minimum, end=maximum, size=step), autobinx=False)])

        if 'log' in self.plotting_options.keys() and self.plotting_options['log'] is True:
            fig.update_layout(yaxis_type="log")
        if 'tickangle' in self.plotting_options.keys():
            fig.update_xaxes(tickangle=self.plotting_options['tickangle'])
        if 'colorbar' in self.plotting_options.keys():
            fig.update_traces(marker_color=self.plotting_options['colorbar'])

        fig.update_layout(title_text=title, title_x=0.5)
        fig.update_xaxes(title_text=xlabel)
        fig.update_yaxes(title_text=ylabel)

        if self.show is True:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)
