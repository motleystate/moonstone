from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph


class ViolinGraph(BaseGraph):

    def plot_one_graph(
        self, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        log_scale: bool = False
    ):
        fig = go.Figure(
            [
                go.Violin(
                    y=self.data,
                    points='all',
                    meanline_visible=True,
                    text=self.data.index,
                    name="All"
                )
            ]
        )

        self._handle_output_plotly(fig, show, output_file, log_scale)


class GroupViolinGraph(GroupBaseGraph):

    def _gen_fig_trace(self, x: list, y: list, name: str, text: list, color: str):
        return go.Violin(
            x=x,
            y=y,
            name=name,
            text=text,
            line_color=color,
            box_visible=True,
            points='all',
            meanline_visible=True,
        )
