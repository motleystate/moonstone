import logging

from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph

logger = logging.getLogger(__name__)


class BaseViolinGraph:

    def _valid_points_param(self, points):
        if points not in ['all', 'outliers', 'suspectedoutliers', False]:
            logger.warning("points=%s not valid, set to default (suspectedoutliers)", points)
            return "suspectedoutliers"
        return points


class ViolinGraph(BaseGraph, BaseViolinGraph):

    def plot_one_graph(
        self, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        points: Union[bool, str] = "suspectedoutliers",
    ):
        fig = go.Figure(
            [
                go.Violin(
                    y=self.data,
                    points=self._valid_points_param(points),
                    meanline_visible=True,
                    text=self.data.index,
                    name="All"
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)


class GroupViolinGraph(GroupBaseGraph, BaseViolinGraph):

    def _gen_fig_trace(
        self, x: list, y: list, name: str, text: list, color: str,
        points: Union[bool, str] = "suspectedoutliers",
    ):
        return go.Violin(
            x=x,
            y=y,
            name=name,
            text=text,
            line_color=color,
            box_visible=True,
            points=self._valid_points_param(points),
            meanline_visible=True,
        )
