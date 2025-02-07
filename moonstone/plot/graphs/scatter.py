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
        fig = go.Figure(
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
        marker: dict = None
    ) -> go.Scatter:
        return go.Scatter(
            x=x,
            y=y,
            name=name,
            text=text,
            marker_color=color,
            mode="markers",
            marker=marker,
        )

    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        group_col: str,
        group_col2: str = None,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        symbols: dict = None,
        groups: list = None,
        groups2: list = None,
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
        if group_col2 and groups2 is None:
            groups2 = list(self.data[group_col2].unique())

        fig = go.Figure()

        if group_col2:
            dic_symbol = self._symbol_scheme(groups2, symbols)
            for gp2 in groups2:
                gp2_filtered_df = self.data[self.data[group_col2] == gp2]
                for gp1 in groups:
                    filtered_df = gp2_filtered_df[gp2_filtered_df[group_col] == gp1]
                    if not filtered_df.empty:
                        fig.add_trace(
                            self._gen_fig_trace(
                                filtered_df[first_col],
                                filtered_df[second_col],
                                str(gp1)+" - "+str(gp2),
                                filtered_df.index,
                                self._get_group_color(gp1, colors),
                                marker=dict(symbol=dic_symbol[gp2]),
                            )
                        )
        else:       # = if group_col2 is None
            for group in groups:
                filtered_df = self.data[self.data[group_col] == group]
                fig.add_trace(
                    self._gen_fig_trace(
                        filtered_df[first_col],
                        filtered_df[second_col],
                        str(group),
                        filtered_df.index,
                        self._get_group_color(group, colors),
                        # **kwargs,  # removing kwargs because I think plotting_options
                        # is already supposed to get the stuff done
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
        fig = go.Figure(
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
        marker: dict = None
    ) -> go.Scatter3d:
        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            name=name,
            text=text,
            marker_color=color,
            mode="markers",
            marker=marker
        )

    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        third_col: str,
        group_col: str,
        group_col2: str = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        symbols: dict = None,
        groups: list = None,
        groups2: list = None,
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
        if group_col2 and groups2 is None:
            groups2 = list(self.data[group_col2].unique())

        fig = go.Figure()

        if group_col2:
            dic_symbol = self._symbol_scheme3d(groups2, symbols)
            for gp2 in groups2:
                gp2_filtered_df = self.data[self.data[group_col2] == gp2]
                for gp1 in groups:
                    filtered_df = gp2_filtered_df[gp2_filtered_df[group_col] == gp1]
                    if not filtered_df.empty:
                        fig.add_trace(
                            self._gen_fig_trace(
                                filtered_df[first_col],
                                filtered_df[second_col],
                                filtered_df[third_col],
                                str(gp1)+" - "+str(gp2),
                                filtered_df.index,
                                self._get_group_color(gp1, colors),
                                marker=dict(symbol=dic_symbol[gp2]),
                            )
                        )
        else:       # = if group_col2 is None
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
