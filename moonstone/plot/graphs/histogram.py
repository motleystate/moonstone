from typing import Union

import plotly.graph_objects as go

from .base import BaseGraph


class Histogram(BaseGraph):

    def plot_one_graph(
        self,
        bins_size: Union[int, float],
        xmin: Union[int, float,str] = 'auto', xmax: Union[int, float,str] = 'auto',
        plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False):
        # see if we add nbinsx options
        if xmin == 'auto':
            minimum = int(self.data.min()/bins_size) * bins_size
        else:
            minimum = xmin
        if xmax == 'auto':
            maximum = self.data.max()
        else:
            maximum = xmax
        fig = go.Figure(
            [
                go.Histogram(
                    x=self.data,
                    xbins=dict(start=minimum, end=maximum, size=bins_size),
                    autobinx=False,
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
