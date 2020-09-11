from typing import Union

import plotly.graph_objects as go

from .base import BaseGraph


class Histogram(BaseGraph):

    def plot_one_graph(self,
                       bins_size: Union[int, float], plotting_options: dict = None,
                       show: bool = True, output_file: Union[bool, str] = False):
        # see if we add nbinsx options
        minimum = int(self.data.min()/bins_size) * bins_size
        maximum = self.data.max()
        fig = go.Figure([go.Histogram(x=self.data,
                                      xbins=dict(start=minimum, end=maximum, size=bins_size),
                                      autobinx=False)])

        fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)
