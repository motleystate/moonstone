from typing import Union

import numpy as np
import pandas as pd
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
    ) -> go.Figure:
        fig = go.Figure(
            self._get_chart(
                orientation=self._valid_orientation_param(orientation),
                ascending=ascending,
                marker_color=marker_color,
                colors_from_string=colors_from_string,
                **kwargs
            )
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)

        return fig


class MatrixBarGraph(BaseGraph):
    """
    Represent a matrix using stacked Bar from plotly.
    """

    def _add_chart(
        self, fig, 
        final_colors, orientation,
        **kwargs
    ) -> go.Figure:
        for col_name in self.data.index:
            if orientation == "v":
                fig.add_trace(
                    go.Bar(
                        name=col_name,
                        x=self.data.columns,
                        y=self.data.loc[col_name],
                        orientation="v",
                        marker_color=final_colors.get(col_name, None),
                    )
                )
            else:
                # horizontal
                fig.add_trace(
                    go.Bar(
                        name=col_name,
                        x=self.data.loc[col_name],
                        y=self.data.columns,
                        orientation="h",
                        marker_color=final_colors.get(col_name, None),
                    ),
                    **kwargs
                )
        return fig

    def _gen_traces_metadata_legends_subplot(
        self,
        fig: go.Figure,
        metadata_ser: pd.Series,
        name: str,
        final_colors_metadata: dict,
        orientation: str
    ) -> go.Figure:
        lbls = list(metadata_ser.unique())
        lbls.sort(reverse=True)
        for lbl in lbls:
            if type(lbl) is not str and np.isnan(lbl):
                dfp = pd.DataFrame(
                    metadata_ser.loc[self.data.columns][
                        metadata_ser.loc[self.data.columns].isna()
                    ]
                )
            else:
                dfp = pd.DataFrame(
                    metadata_ser.loc[self.data.columns][
                        metadata_ser.loc[self.data.columns] == lbl
                    ]
                )
            dfp["y"] = 1
            if orientation == "v":
                fig.add_trace(
                    go.Bar(
                        x=dfp.index,
                        y=dfp["y"],
                        name=lbl,
                        marker=dict(color=final_colors_metadata[lbl]),
                        legendgroup=name,
                        legendgrouptitle_text=name,
                    ),
                    row=2,
                    col=1,
                )
            else:
                # horizontal
                if orientation == "h-l":
                    fig.add_trace(
                        go.Bar(
                            x=dfp["y"],
                            y=dfp.index,
                            orientation="h",
                            name=lbl,
                            marker=dict(color=final_colors_metadata[lbl]),
                            legendgroup=name,
                            legendgrouptitle_text=name,
                        ),
                        row=1,
                        col=1,
                    )
                else:
                    fig.add_trace(
                        go.Bar(
                            x=dfp["y"],
                            y=dfp.index,
                            orientation="h",
                            name=lbl,
                            marker=dict(color=final_colors_metadata[lbl]),
                            legendgroup=name,
                            legendgrouptitle_text=name,
                        ),
                        row=1,
                        col=2,
                    )
        return fig

    def _gen_traces_metadata_legends_subplots(
        self, fig: go.Figure, metadata_df: pd.DataFrame, final_colors_metadata: dict, orientation: str
    ) -> go.Figure:
        for cc in metadata_df.columns:
            fig = self._gen_traces_metadata_legends_subplot(
                fig, metadata_df[cc], cc, final_colors_metadata, orientation
            )
        return fig

    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        orientation: str = "v",
    ) -> go.Figure:
        """
        Args:
            colors: Selected colors for a group.
            orientation: orientation of the graph. {"v" (or "vertical")(default), "h-l" (or "horizontal-left"),
              "h-r" (or "horizontal-right"}.
        """
        orientation = self._valid_orientation_param(orientation, hplus=True)
        
        final_colors = self._color_scheme_species(colors)

        fig = go.Figure()
        fig = self._add_chart(fig, final_colors, orientation)
        if orientation == "h-r":
            fig.update_xaxes(autorange="reversed")

        fig.update_layout(barmode="stack", legend_traceorder="reversed")
        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig

    def plot_complex_graph(
        self,
        metadata: Union[pd.DataFrame, pd.Series],
        plotting_options: dict = {},
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        colors_metadata: dict = None,
        orientation: str = "v",
    ) -> go.Figure:
        """
        Args:
            metadata: pandas dataframe or series with the metadata relevant to show below/side-to-side to the bar graph.
            colors: Selected colors for a group in the bar graph part of the graph.
            orientation: orientation of the graph. {"v" (or "vertical")(default), "h-l" (or "horizontal-left"),
              "h-r" (or "horizontal-right"}.
        """
        # metadata = samples (row) * metadata (col)
        # data = species * samples

        orientation = self._valid_orientation_param(orientation, hplus=True)

        final_colors = self._color_scheme_species(colors)  # attribute a color to each species
        final_colors_metadata = self._color_scheme_metadata(metadata, colors_metadata)

        if isinstance(metadata, pd.Series):
            nrows = 1
        else:
            nrows = len(metadata.columns)

        if orientation == "v":
            fig = make_subplots(
                rows=2,
                cols=1,
                shared_xaxes=True,
                vertical_spacing=0.02,
                start_cell='top-left',  # default
                row_heights=[1 - (0.02 * nrows), 0.02 * nrows],
            )
        else:
            # horizontal
            if orientation == "h-l":
                fig = make_subplots(
                    rows=1,
                    cols=2,
                    shared_yaxes=True,
                    horizontal_spacing=0.02,
                    column_widths=[0.02 * nrows, 1 - (0.02 * nrows)],  # start from left
                )
            else:
                fig = make_subplots(
                    rows=1,
                    cols=2,
                    shared_yaxes=True,
                    horizontal_spacing=0.02,
                    column_widths=[1 - (0.02 * nrows), 0.02 * nrows],  # start from left
                )

        # main graph
        if orientation == "h-l":
            fig = self._add_chart(fig, final_colors, orientation, row=1, col=2)
        else:
            fig = self._add_chart(fig, final_colors, orientation)
            if orientation == "h-r":
                fig.update_xaxes(autorange="reversed")

        # metadata "legends" subplot.s
        if isinstance(metadata, pd.Series):
            fig = self._gen_traces_metadata_legends_subplot(
                fig, metadata, metadata.name, final_colors_metadata, orientation
            )
        else:
            fig = self._gen_traces_metadata_legends_subplots(
                fig, metadata, final_colors_metadata, orientation
            )
        
        if "layout" in plotting_options.keys():
            xaxis_title = plotting_options["layout"].pop("xaxis_title", None)
            fig.update_layout(
                xaxis2=dict(  # xaxis of the 2nd subplot (to not have "samples" * 2)
                    title_text=xaxis_title,
                )
            )

        if orientation == "v":
            fig.update_layout(
                yaxis2=dict(showticklabels=False),  # yaxis of the 2nd subplot
            )
        else:
            # horizontal
            if orientation == "h-l":
                fig.update_layout(
                    xaxis1=dict(showticklabels=False),  # xaxis of the 1st subplot
                )
            else:
                # h-r
                fig.update_layout(
                    xaxis2=dict(showticklabels=False),  # xaxis of the 2nd subplot
                )

        if "layout" in plotting_options.keys() and "legend" in plotting_options["layout"].keys():
            plotting_options["layout"]["legend"].pop("traceorder", None)
            # traceorder: "normal" given by `plot_sample_composition_most_abundant_taxa`
        fig.update_layout(barmode="stack")

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
