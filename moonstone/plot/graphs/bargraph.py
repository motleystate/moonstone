from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph


class BarGraph(BaseGraph):

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

    def plot_one_graph(self, plotting_options: dict = None,
                       orientation: str = "v", ascending: bool = None,
                       show: bool = True, output_file: Union[bool, str] = False):
        fig = go.Figure(self._get_chart(orientation=orientation, ascending=ascending))

        fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)
