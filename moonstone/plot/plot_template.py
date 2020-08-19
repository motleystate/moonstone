from abc import ABC
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io
from typing import Union


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

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        """
        method that plots the graph
        """
        # instantiate fig

        # treat every plotting options by updating layout

        # if self.show is True:
        #    fig.show()

        # if self.output_file:
        #    plotly.io.write_html(fig, self.output_file)
        pass


class BarGraph(BaseGraph):

    def compute_heterogeneous_bins(self):
        """Logish bins"""
        maximum = self.df.max()
        magnitude = int(math.log10(maximum))
        self.bins_values = [-0.1]                       # to have the 0 value
        i = 0
        while i < magnitude:
            self.bins_values += list(np.arange(2*10**i, 10**(i+1)+1, 10**i))
            i += 1
        # i=magnitude
        self.bins_values += list(np.arange(2*10**i, maximum+10**i, 10**i))        # up to maximum

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

    def reset_xnames(self, replace_dic: dict):
        for i in range(len(self.xnames)):
            self.xnames[i] = replace_dic[self.xnames[i]]

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        fig = go.Figure([go.Bar(x=self.xnames, y=self.yvalues)])

        if 'log' in self.plotting_options.keys() and self.plotting_options['log']:
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
