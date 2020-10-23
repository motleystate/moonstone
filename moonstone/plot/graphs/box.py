from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph


class BoxGraph(BaseGraph):

    def plot_one_graph(
        self, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
    ):
        fig = go.Figure(
            [
                go.Box(
                    y=self.data,
                    boxpoints='all',
                    text=self.data.index,
                    name="All",
                    marker_color=self.DEFAULT_COLOR,
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)


class GroupBoxGraph(GroupBaseGraph):

    def _gen_fig_trace(self, x: list, y: list, name: str, text: list, color: str):
        return go.Box(
            x=x,
            y=y,
            name=name,
            text=text,
            marker_color=color,
            boxpoints='all',
        )
