from typing import Union

import plotly.graph_objects as go

from .base import BaseGraph
from moonstone.utils.plot import (
    add_x_to_plotting_options
)


class Histogram(BaseGraph):

    def _get_chart(self, bins_size: Union[int, float] = None) -> go.Histogram:
        if bins_size is not None:
            return go.Histogram(x=self.data,
                                xbins=dict(size=bins_size),
                                autobinx=False)
        return go.Histogram(x=self.data)

    def plot_one_graph(self, bins_size: Union[int, float] = None,
                       xmin: Union[int, float, str] = None, plotting_options: dict = None,
                       show: bool = True, output_file: Union[bool, str] = False):
        # see if we add nbinsx options

        fig = go.Figure(self._get_chart(bins_size=bins_size))

        if plotting_options is None and (bins_size is not None or xmin is not None):
            plotting_options = {}
        if plotting_options is not None:
            if bins_size is not None:
                plotting_options = add_x_to_plotting_options(plotting_options, 'xaxes', 'dtick', bins_size)
            if xmin is not None:
                xmax = (int(self.data.max() / bins_size) + 1) * bins_size
                plotting_options = add_x_to_plotting_options(plotting_options, 'xaxes', 'range', [xmin, xmax])
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
