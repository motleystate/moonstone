from abc import ABC, abstractmethod
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io
from typing import Union


class BaseGraph(ABC):

    def __init__(self, dataframe: pd.DataFrame):
        """
        :param dataframe: pandas dataframe generated with moonstone
        """
        self.df = dataframe
        self.fig = go.Figure()

    @abstractmethod
    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass

    def _handle_colorbar_plotly(self, value):
        self.fig.update_traces(marker_color=value)

    def _handle_log_plotly(self, value):
        if "xaxis" in value:
            self.fig.update_layout(xaxis_type="log")
        elif "yaxis" in value or value is True:
            self.fig.update_layout(yaxis_type="log")

    def _handle_tickangle_plotly(self, value):
        self.fig.update_xaxes(tickangle=value)

    def _handle_shapes_plotly(self, value):
        """
        :param value: list of dictionnary for all the shapes to draw, with values for type,
                      x0, y0, x1, y1, xref [optional], yref [optional],
                      line dictionnary [optional]
                      -> line dictionnary contains lines descriptor like width, color,
                         or dash (dash, dashdot, dot etc.)
        """
        self.fig.update_layout(shapes=value)

    def _handle_plotting_options_plotly(self, plotting_options):
        for option in plotting_options.keys():
            handler = f"_handle_{option}_plotly"
            getattr(self, handler)(plotting_options[option])

    def _display_titles_plotly(self, title, xlabel, ylabel):
        self.fig.update_layout(title_text=title, title_x=0.5)
        self.fig.update_xaxes(title_text=xlabel)
        self.fig.update_yaxes(title_text=ylabel)

    def _handle_output_plotly(self, show, output_file):
        if show is True:
            self.fig.show()

        if output_file:
            if output_file is True:
                # if no name given for the output file, a generic name is generated
                if self.df.name is not None:
                    output_file = self.df.name+"_"+self.__class__.__name__+".html"
                else:
                    output_file = self.__class__.__name__+".html"
            plotly.io.write_html(self.fig, output_file)


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

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str, plotting_options: dict,
                       show: bool = True, output_file: Union[bool, str] = False):
        self.fig = go.Figure([go.Bar(x=self.xnames, y=self.yvalues)])

        self._handle_plotting_options_plotly(plotting_options)
        self._display_titles_plotly(title, xlabel, ylabel)
        self._handle_output_plotly(show, output_file)


class Histogram(BaseGraph):

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str,
                       bins_size: Union[int, float], plotting_options: dict = None,
                       show: bool = True, output_file: Union[bool, str] = False):
        # see if we add nbinsx options
        minimum = int(self.df.min()/bins_size) * bins_size
        maximum = self.df.max()
        self.fig = go.Figure([go.Histogram(x=self.df, xbins=dict(start=minimum, end=maximum, size=bins_size),
                                           autobinx=False)])

        self._handle_plotting_options_plotly(plotting_options)
        self._display_titles_plotly(title, xlabel, ylabel)
        self._handle_output_plotly(show, output_file)
