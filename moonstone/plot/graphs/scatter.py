import logging
from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph, GroupBaseGraph

logger = logging.getLogger(__name__)


class ScatterGraph(BaseGraph):
    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
    ) -> go.Figure:
        """
        :param first_col: col name for x data
        :param second_col: col name for y data
        """
        fig = go.Fig(
            [
                go.Scatter(
                    x=self.data[first_col],
                    y=self.data[second_col],
                    text=self.data.index,
                    mode="markers",
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


class GroupScatterGraph(GroupBaseGraph):
    def _gen_fig_trace(
        self,
        x: list,
        y: list,
        name: str,
        text: list,
        color: str,
    ) -> go.Scatter:
        return go.Scatter(
            x=x,
            y=y,
            name=name,
            text=text,
            marker_color=color,
            mode="markers",
        )

    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        group_col: str,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        groups: list = None,
        **kwargs
    ) -> go.Figure:
        """
        :param first_col: col name for x data
        :param second_col: col name for y data
        :param group_col: column used to group data
        """
        if groups is None:
            groups = list(self.data[group_col].unique())
        if colors is None:
            colors = self._gen_default_color_dict(groups)
        fig = go.Figure()

        for group in groups:
            filtered_df = self.data[self.data[group_col] == group]
            fig.add_trace(
                self._gen_fig_trace(
                    filtered_df[first_col],
                    filtered_df[second_col],
                    str(group),
                    filtered_df.index,
                    self._get_group_color(group, colors),
                    **kwargs,
                )
            )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


class Scatter3DGraph(BaseGraph):
    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        third_col: str,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
    ) -> go.Figure:
        """
        :param first_col: col name for x data
        :param second_col: col name for y data
        :param third_col: col name for z data
        """
        fig = go.Fig(
            [
                go.Scatter3d(
                    x=self.data[first_col],
                    y=self.data[second_col],
                    z=self.data[third_col],
                    text=self.data.index,
                    mode="markers",
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


class GroupScatter3DGraph(GroupBaseGraph):
    def _gen_fig_trace(
        self,
        x: list,
        y: list,
        z: list,
        name: str,
        text: list,
        color: str,
    ) -> go.Scatter3d:
        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            name=name,
            text=text,
            marker_color=color,
            mode="markers",
        )

    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        third_col: str,
        group_col: str,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        groups: list = None,
        plotting_options: dict = None,
        **kwargs
    ) -> go.Figure:
        """
        :param first_col: col name for x data
        :param second_col: col name for y data
        :param third_col: col name for z data
        :param group_col: column used to group data
        """
        if groups is None:
            groups = list(self.data[group_col].unique())
        if colors is None:
            colors = self._gen_default_color_dict(groups)
        fig = go.Figure()

        for group in groups:
            filtered_df = self.data[self.data[group_col] == group]
            fig.add_trace(
                self._gen_fig_trace(
                    filtered_df[first_col],
                    filtered_df[second_col],
                    filtered_df[third_col],
                    str(group),
                    filtered_df.index,
                    self._get_group_color(group, colors),
                    **kwargs,
                )
            )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
