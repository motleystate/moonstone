from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph


class BoxGraph(BaseGraph):

    def plot_one_graph(
        self, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        log_scale: bool = False
    ):
        fig = go.Figure(
            [
                go.Box(
                    y=self.data,
                    boxpoints='all',
                    text=self.data.index,
                    name="All"
                )
            ]
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file, log_scale)


class GroupBoxGraph(BaseGraph):

    def plot_one_graph(
        self, data_col: str, group_col: str, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        log_scale: bool = False
    ):
        """
        :param data_col: column with data to visualize
        :param group_col: column used to group data
        """
        groups = list(self.data[group_col].unique())
        groups.sort()
        fig = go.Figure()
        for group in groups:
            fig.add_trace(go.Box(x=self.data[group_col][self.data[group_col] == group],
                                 y=self.data[data_col][self.data[group_col] == group],
                                 name=str(group),
                                 boxpoints='all',
                                 text=self.data.index))

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file, log_scale)
