import logging

from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph

logger = logging.getLogger(__name__)


class BaseBoxGraph:
    def _valid_boxpoints_param(self, boxpoints: str) -> str:
        if boxpoints not in ["all", "outliers", "suspectedoutliers", False]:
            logger.warning(
                "boxpoints=%s not valid, set to default (suspectedoutliers)", boxpoints
            )
            return "suspectedoutliers"
        return boxpoints


class BoxGraph(BaseGraph, BaseBoxGraph):
    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        boxpoints: Union[bool, str] = "suspectedoutliers",
        marker_color: str = None,
    ) -> go.Figure:
        """
        Plot box graph.

        Args:
            plotting_options: options for the layout of the graph
            show: whether or not to show the graph
            output_file: path for output file
            boxpoints: boxpoints param from plotly
            marker_color: color for the boxplot
        """
        marker_color = self.DEFAULT_COLOR if marker_color is None else marker_color
        fig = go.Figure(
            [
                go.Box(
                    y=self.data,
                    boxpoints=self._valid_boxpoints_param(boxpoints),
                    text=self.data.index,
                    name="All",
                    marker_color=marker_color,
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


class GroupBoxGraph(GroupBaseGraph, BaseBoxGraph):
    def _gen_fig_trace(
        self,
        x: list,
        y: list,
        name: str,
        text: list,
        color: str,
        boxpoints: Union[bool, str] = "suspectedoutliers",
    ) -> go.Box:
        return go.Box(
            x=x,
            y=y,
            name=name,
            text=text,
            marker_color=color,
            boxpoints=self._valid_boxpoints_param(boxpoints),
        )
