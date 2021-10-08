from typing import Union

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
    def _add_chart(
        self, fig,
        final_colors,
        **kwargs
    ) -> go.Bar:
        for col_name in self.data.index:
            fig.add_trace(
                go.Bar(
                    name=col_name,
                    x=self.data.columns,
                    y=self.data.loc[col_name],
                    marker_color=final_colors.get(col_name, None)
                )
            )
        return fig

    def _color_scheme(self, colors: dict = None):
        final_colors = {name: generate_color_code(name) for name in self.data.index}
        if colors is not None:
            final_colors.update(**colors)
        return final_colors

    def _color_scheme_metadata(self, metadata, colors: dict = None):
        final_colors = {}
        all_gp = []
        for cc in metadata.columns:
            all_gp += list(metadata[cc].unique())
        all_gp = list(set(all_gp))  # doesn't remove all different nan
# we can't do it manually because then some nan don't correspond to the one manually added
        if len(all_gp) <= 10:
            final_colors = dict(zip(all_gp, px.colors.qualitative.Plotly))
        elif len(all_gp) <= 26:
            final_colors = dict(zip(all_gp, px.colors.qualitative.Alphabet))
        else:
            c = px.colors.qualitative.Alphabet+px.colors.qualitative.Set3
            final_colors = dict(zip(all_gp, c*(int(len(all_gp)/len(c))+1)))
        if colors is not None:
            final_colors.update(**colors)
        return final_colors

    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
    ):
        """
        Args:
            colors: Selected colors for a group
        """
        final_colors = self._color_scheme(colors)

        fig = go.Figure()
        fig = self._add_chart(fig, final_colors)

        fig.update_layout(barmode="stack")
        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig

    def plot_complex_graph(
        self,
        metadata,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        colors_metadata: dict = None
    ):
        # metadata = samples (row) * metadata (col)
        # data = species * samples

        final_colors = self._color_scheme(colors)
        final_colors_metadata = self._color_scheme_metadata(metadata, colors_metadata)

        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.02,
            row_width=[0.02 * len(metadata.columns), 1 - (0.02 * len(metadata.columns))]
        )

        # main graph
        fig = self._add_chart(fig, final_colors)

        # metadata "legends" subplot.s
        for cc in metadata.columns:
            for lbl in metadata[cc].unique():
                if type(lbl) != str and np.isnan(lbl):
                    dfp = metadata.loc[self.data.columns][metadata.loc[self.data.columns][cc].isna()]
                else:
                    dfp = metadata.loc[self.data.columns][metadata.loc[self.data.columns][cc] == lbl]
                dfp['y'] = 1
                fig.add_trace(
                    go.Bar(x=dfp.index,
                           y=dfp['y'],
                           name=lbl,
                           marker=dict(
                               color=final_colors_metadata[lbl]
                           ),
                           legendgroup=cc,
                           legendgrouptitle_text=cc,
                           ),
                    row=2, col=1)

        if 'layout' in plotting_options.keys():
            xaxis_title = plotting_options['layout'].pop("xaxis_title", "Samples")
            if 'legend' in plotting_options['layout'].keys():
                plotting_options['layout']['legend'].pop("traceorder", None)

        fig.update_layout(
            xaxis2=dict(                             # xaxis of the 2nd subplot (to not have "samples" * 2)
                title_text=xaxis_title,
            ),
            yaxis2=dict(                             # yaxis of the 2nd subplot
                showticklabels=False
            ),
            barmode="stack"
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
