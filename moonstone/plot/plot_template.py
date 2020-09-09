from abc import ABC, abstractmethod
import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io
from typing import Union

from moonstone.plot.functions import (
    check_list_type
)


class BaseGraph(ABC):

    def __init__(self, data: Union[pd.Series, pd.DataFrame], plotting_options: dict = None,
                 show: bool = True, output_file: Union[bool, str] = False):
        """
        :param data: data to plot
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        self.data = data
        self.plotting_options = plotting_options

        if type(show) == bool:
            self.show = show
        else:
            raise ValueError('Error : show value must be a bool, %s given' % type(show).__name__)

        self.output_file = output_file
        if self.output_file is True:
            # if no name given for the output file, a generic name is generated
            if data.name is not None:
                self.output_file = data.name+"_"+self.__class__.__name__+".html"
            else:
                self.output_file = self.__class__.__name__+".html"

    @abstractmethod
    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass


class BaseBarGraph(BaseGraph, ABC):

    def _set_plotting_options(self, fig, plotting_options):
        if 'log' in self.plotting_options.keys() and self.plotting_options['log']:
            fig.update_layout(yaxis_type="log")
        if 'tickangle' in self.plotting_options.keys():
            fig.update_xaxes(tickangle=self.plotting_options['tickangle'])
        if 'colorbar' in plotting_options.keys():
            fig.update_traces(marker_color=plotting_options['colorbar'])
        return fig

    def _set_labels(self, fig, title, xlabel, ylabel):
        if title is not None:
            fig.update_layout(title_text=title, title_x=0.5)
        if xlabel is not None:
            fig.update_xaxes(title_text=xlabel)
        if ylabel is not None:
            fig.update_yaxes(title_text=ylabel)
        return fig

    def plot_one_graph(self, orientation: str = "v", ascending: bool = None,
                       title: str = None, xlabel: str = None, ylabel: str = None,
                       plotting_options: dict = None):
        """
        :param title: title of the graph
        :param xlabel: label of the x axis
        :param ylabel: label of the y axis
        :param reset_xnames_dic: to rename the names of the values in the x axis. \n
                                 Example for a plot of the distribution of smoking habits :
                                 reset_xnames_dic={'y': 'smoker', 'n': 'non smoker'}
        """
        fig = go.Figure(self._get_chart(orientation=orientation, ascending=ascending))

        if plotting_options is not None:
            fig = self._set_plotting_options(fig, plotting_options)

        fig = self._set_labels(fig, title, xlabel, ylabel)

        if self.show is True:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)


class CategoryBarGraph(BaseBarGraph):

    def _get_chart(self, orientation: str = "v", ascending: bool = None) -> go.Bar:
        if ascending is not None:
            data = self.data.sort_values(ascending=ascending)
        else:
            data = self.data
        x = list(data.index)
        y = list(data)
        if orientation == "v":
            return go.Bar(x=x, y=y, orientation=orientation)
        return go.Bar(x=y, y=x, orientation=orientation)


class DistributionBarGraph(BaseBarGraph):

    def compute_heterogeneous_bins(self):
        """Logish bins"""
        maximum = self.data.max()
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
        if type(bins_values) == list and check_list_type(bins_values, (int, float, np.integer)):
            self._bins_values = bins_values
        else:
            raise ValueError("Error : expecting list of numerical values (int, float) in bins_values.")

    def in_bins_and_count(self, normalize: bool = False):
        binned_df = pd.cut(self.data, bins=self.bins_values)   # put every items in the appropriate bin
        data = pd.value_counts(binned_df, normalize=normalize)
        data = data.reindex(binned_df.cat.categories)
        new_xnames = list(data.index.astype(str))
        new_xnames[0] = new_xnames[0].replace('(-0.1', '[0.0')
        new_xnames = [new_xnames[i].replace('(', ']') for i in range(len(new_xnames))]
        data.index = new_xnames
        return data

    def _get_chart(self, orientation: str = "v", ascending: bool = None) -> go.Bar:
        data = self.in_bins_and_count()
        if ascending is not None:
            data = data.sort_values(ascending=ascending)
        x = list(data.index)
        y = list(data)
        if orientation == "v":
            return go.Bar(x=x, y=y, orientation=orientation)
        return go.Bar(x=y, y=x, orientation=orientation)


class Histogram(BaseGraph):

    def plot_one_graph(self, title: str, xlabel: str, ylabel: str,
                       step: Union[int, float]):   # see if we add nbinsx options
        minimum = int(self.data.min()/step) * step
        maximum = self.data.max()
        fig = go.Figure([go.Histogram(x=self.data, xbins=dict(start=minimum, end=maximum, size=step), autobinx=False)])

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
