from typing import Union

import plotly.graph_objects as go

import logging

from .base import BaseGraph, GroupBaseGraph

logger = logging.getLogger(__name__)


class Histogram(BaseGraph):

    def _get_chart(self, xbins: dict = None) -> go.Histogram:
        return go.Histogram(x=self.data, xbins=xbins)

    def plot_one_graph(self, xbins: dict = None,
                       plotting_options: dict = None,
                       show: bool = True, output_file: Union[bool, str] = False):
        """
        Args:
            xbins: `plotly.graph_objects.histogram.XBins 
            <https://plotly.github.io/plotly.py-docs/generated/plotly.graph_objects.histogram.html#plotly.graph_objects.histogram.XBins>`_
            instance or dict with compatible properties.
        """  # noqa
        fig = go.Figure(self._get_chart(xbins=xbins))
        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)
        return fig


class GroupHistogram(GroupBaseGraph, Histogram):
    def _gen_fig_trace(
        self,
        x: list,
        name: str,
        color: str,
        xbins: dict = None,
        **kwargs,
    ) -> go.Histogram:
        return go.Histogram(
            x=x,
            name=name,
            marker_color=color,
            xbins=xbins,
            **kwargs,
        )

    def plot_one_graph(
        self,
        data_col: str,
        group_col: str,
        groups: list = None,
        sort_groups: bool = False,
        colors: dict = None,
        bins_size: Union[int, float] = None,
        plotting_options: dict = None,
        show_counts: bool = False,
        show: bool = True,
        output_file: Union[bool, str] = False,
        **kwargs,
    ) -> go.Figure:
        """
        Args:
            data_col: column with data to visualize
            group_col: column used to group data
            groups: specifically select groups to display among group_col
            sort_groups: whether to sort groups to display or not
            colors: dictionnary with group_col2 (or group_col if no group_col2) values as keys
            bins_size: size of the histogram bins
            show_count: write the number of samples per group in the groups' name (True) or not (False; default)
        """
        fig = go.Figure()

        groups, colors, names = self._prep_for_plot_one_graph(
            group_col, groups, colors, show_counts, sort_groups
        )

        # bin_size = kwargs.get('bin_size', kwargs.get('size', None))
        # Putting it as an hidden argument for now because I don't want to take time to make it robust yet
        if 'xbins' in kwargs.keys():
            logger.warning("using the xbins found in the kwargs.")
            xbins = kwargs['xbins']
        else:
            xbins = dict(start=self.data[data_col].min(), end=self.data[data_col].max(), size=bins_size)
        for group in groups:
            filtered_df = self.data[self.data[group_col] == group]
            fig = fig.add_trace(self._gen_fig_trace(
                x=filtered_df[data_col],
                name=str(names[group]),
                color=self._get_group_color(group, colors),
                xbins=xbins,
                opacity=0.6, hoverinfo="x+y",
                **kwargs,
            ))

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)

        return fig
