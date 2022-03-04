from typing import Union, List

import plotly.graph_objects as go

from .base import BaseGraph


class HeatmapGraph(BaseGraph):
    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colorscale: List[list] = None,
    ) -> go.Figure:
        fig = go.Figure(
            [
                go.Heatmap(
                    x=self.data.columns,
                    y=self.data.index,
                    z=self.data,
                    text=self.data,
                    colorscale=colorscale,
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
