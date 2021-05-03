from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph
from moonstone.utils.colors import generate_color_code


class BarGraph(BaseGraph):
    def _get_chart(
        self,
        orientation: str = "v",
        ascending: bool = None,
        marker_color: str = "crimson",
        colors_from_string: bool = False,
        **kwargs
    ) -> go.Bar:
        if ascending is not None:
            data = self.data.sort_values(ascending=ascending)
        else:
            data = self.data
        x = list(data.index)
        y = list(data)
        if colors_from_string:
            marker_color = [generate_color_code(name) for name in data.index]
        if orientation == "v":
            return go.Bar(
                x=x, y=y, orientation=orientation, marker_color=marker_color, **kwargs
            )
        return go.Bar(
            x=y, y=x, orientation=orientation, marker_color=marker_color, **kwargs
        )

    def plot_one_graph(
        self,
        plotting_options: dict = None,
        orientation: str = "v",
        ascending: bool = None,
        marker_color: str = "crimson",
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors_from_string: bool = False,
        **kwargs
    ):
        fig = go.Figure(
            self._get_chart(
                orientation=orientation,
                ascending=ascending,
                marker_color=marker_color,
                colors_from_string=colors_from_string,
                **kwargs
            )
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)


class MatrixBarGraph(BaseGraph):
    """
    Represent a matrix using stacked Bar from plotly.
    """

    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None
    ):
        """
        Args:
            colors: Selected colors for a group
        """
        final_colors = {name: generate_color_code(name) for name in self.data.index}
        if colors is not None:
            final_colors.update(**colors)
        fig = go.Figure()
        for col_name in self.data.index:
            fig.add_trace(
                go.Bar(name=col_name, x=self.data.columns, y=self.data.loc[col_name], marker_color=final_colors.get(col_name, None))
            )

        fig.update_layout(barmode="stack")
        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)
