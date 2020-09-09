from typing import Union

import plotly.graph_objects as go
import plotly.io

from .base import BaseGraph


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
