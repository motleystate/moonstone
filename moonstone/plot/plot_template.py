import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io


class BaseGraph():

    def __init__(self, dataframe, plotting_options, **kwargs):
        self.df = dataframe
        self.plotting_options = plotting_options
        if 'show' in kwargs.keys():
            if type(kwargs['show']) == bool:
                self.show = kwargs['show']
            else:
                raise ValueError('Error : show value must be a bool, %s given' % type(kwargs['show']).__name__)
        if 'output_file' in kwargs.keys() and kwargs['output_file']:
            self.output_file = kwargs['output_file']
            if self.output_file is True:
                # no name given for the output file so a generic name will be given
                self.output_file = "graph.html"
        else:
            self.output_file = False


class BarGraph(BaseGraph):

    def compute_symetric_bins(self, step):
        maximum = self.df.max()
        minimum = int(self.df.min()/step) * step
        self.bins_values = [minimum]
        i = 1
        while self.bins_values[i-1] < maximum:
            self.bins_values += [self.bins_values[i-1]+step]
            i += 1

    def compute_asymetric_bins(self):
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

    def in_bins_and_count(self, normalize=False):
        binned_df = pd.cut(self.df, bins=self.bins_values)   # put every items in the appropriate bin
        counts = pd.value_counts(binned_df, normalize=normalize)
        counts = counts.reindex(binned_df.cat.categories)

        self.xnames = list(counts.index.astype(str))
        self.xnames[0] = self.xnames[0].replace('(-0.1', '[0.0')
        self.xnames = [self.xnames[i].replace('(', ']') for i in range(len(self.xnames))]
        self.yvalues = list(counts)

    def count(self, normalize=False):
        counts = pd.value_counts(self.df, normalize=normalize)
        counts = counts.sort_index()

        self.xnames = list(counts.index)
        self.yvalues = list(counts)

    def reset_xnames(self, replace_dic: dict):
        for i in range(len(self.xnames)):
            self.xnames[i] = replace_dic[self.xnames[i]]

    def plot_one_graph(self, title, xlabel, ylabel):
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

        if self.show:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)


class Histogram(BaseGraph):

    def plot_one_graph(self, title, xlabel, ylabel, step):   # see if we add nbinsx options
        minimum = int(self.df.min()/step) * step
        maximum = self.df.max()
        fig = go.Figure([go.Histogram(x=self.df, xbins=dict(start=minimum, end=maximum, size=5), autobinx=False)])

        if 'log' in self.plotting_options.keys() and self.plotting_options['log'] is True:
            print("ici")
            fig.update_layout(yaxis_type="log")
        if 'tickangle' in self.plotting_options.keys():
            fig.update_xaxes(tickangle=self.plotting_options['tickangle'])
        if 'colorbar' in self.plotting_options.keys():
            fig.update_traces(marker_color=self.plotting_options['colorbar'])

        fig.update_layout(title_text=title, title_x=0.5)
        fig.update_xaxes(title_text=xlabel)
        fig.update_yaxes(title_text=ylabel)

        if self.show:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)
