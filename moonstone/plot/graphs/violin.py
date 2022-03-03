import logging

from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph

logger = logging.getLogger(__name__)


class BaseViolinGraph:
    def _valid_points_param(self, points: str) -> str:
        if points not in ["all", "outliers", "suspectedoutliers", False]:
            logger.warning(
                "points=%s not valid, set to default (suspectedoutliers)", points
            )
            return "suspectedoutliers"
        return points


class ViolinGraph(BaseGraph, BaseViolinGraph):
    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        points: Union[bool, str] = "suspectedoutliers",
        box: bool = False,
        color: str = None,
        line_color: str = None,
    ) -> go.Figure:
        color = self.DEFAULT_COLOR if color is None else color
        line_color = self.DEFAULT_COLOR if line_color is None else line_color
        fig = go.Figure(
            [
                go.Violin(
                    y=self.data,
                    points=self._valid_points_param(points),
                    meanline_visible=True,
                    text=self.data.index,
                    box_visible=box,
                    name="All",
                    fillcolor=color,
                    line_color=line_color,
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


class GroupViolinGraph(GroupBaseGraph, BaseViolinGraph):
    def _gen_fig_trace(
        self,
        x: list,
        y: list,
        name: str,
        text: list,
        color: str,
        points: Union[bool, str] = "suspectedoutliers",
        box: bool = False,
    ) -> go.Violin:
        return go.Violin(
            x=x,
            y=y,
            name=name,
            text=text,
            line_color=color,
            points=self._valid_points_param(points),
            meanline_visible=True,
            box_visible=box,
        )
