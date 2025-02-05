import logging
from typing import Union

import plotly.graph_objects as go
from plotly.validators.scatter.marker import SymbolValidator

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

    @property
    def dictionary_marker_symbols(self):
        """
        All the symbols accepted by scatter as integer.
        """
        if getattr(self, '_dictionary_marker_symbols', None) is None:
            raw_symbols = SymbolValidator().values
            self._dictionary_marker_symbols = {}
            for i in range(0, len(raw_symbols), 3):
                self._dictionary_marker_symbols[raw_symbols[i+2]] = raw_symbols[i]
        return self._dictionary_marker_symbols

    def _translating_forced_symbols(self, forced_symbols: Union[list, set]) -> set:
        translated_symbols = []
        txt_symbols = []
        for x in forced_symbols:
            if type(x) is int:
                translated_symbols += [x]
            elif x.isnumeric():
                translated_symbols += [int(x)]
            else:
                txt_symbols += [x]

        if len(txt_symbols) > 0:
            for i in txt_symbols:
                if i in self.dictionary_marker_symbols.keys():
                    translated_symbols += [self.dictionary_marker_symbols[i]]
                    # txt_symbols.remove(i)
                else:
                    # raise error
                    error_message = f"Invalid symbol values: {txt_symbols}.\n\
Accepted values: {SymbolValidator().values}"
                    raise ValueError(error_message)
        return set(translated_symbols)

    def _symbol_scheme(self, groups: list, symbols: dict = None):
        """
        Generate a dictionary with the different groups as key and an int corresponding to a symbol (following
        go.Marker's definition) as value.
        If number of different groups > 24, then symbols are repeated.

        Args:
            metadata: Can take DataFrame or Series, then all different values are listed and assigned a color.
                Or can take a list/array, then all elements of the list are assigned a color.
                Finally, metadata can be an integer which act like a list ranging from 0 to that integer (excluded),
                and each integer is assigned a color.
            symbols: Dictionary that dictate a particular symbol to a particular group.
        """
        if len(groups) > 162:
            list_symbols = list(self.dictionary_marker_symbols.values())
            dic_symbols = dict(zip(groups, list_symbols * (int(len(groups) / len(list_symbols)) + 1)))
        else:
            if len(groups) <= 72:
                list_fill_pattern = ["", "1", "3"]     # 2 is filled with a dot inside ("-dot")
                # but it is not very clear
                list_shape = [
                    "00", "01", "02", "03", "04", "05", "06", "17", "13", "18", "07", "08", "14", "23", "24", "19",
                    "20", "21", "22", "16", "09", "10", "11", "12", "15",
                ]
                n = -(-len(groups) // 3)
                if n <= 8:
                    n = 8

                list_symbols = [int(prefix + suffix) for prefix in list_fill_pattern for suffix in list_shape[:n]]
            else:
                list_symbols = self.dictionary_marker_symbols.values()

            if symbols is not None:
                to_remove = self._translating_forced_symbols(set(symbols.values()))

                # for reusability, we can't use set().difference()
                list_symbols = [x for x in list_symbols if x not in to_remove]
                skeys = set(symbols.keys())
                groups = [x for x in groups if x not in skeys]
            dic_symbols = dict(zip(groups, list_symbols))

        if symbols is not None:
            dic_symbols.update(**symbols)

        return dic_symbols

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
            dic_symbol = self._symbol_scheme(groups2)
            for gp2 in groups2:
                gp2_filtered_df = self.data[self.data[group_col2] == gp2]
                for gp1 in groups:
                    filtered_df = gp2_filtered_df[gp2_filtered_df[group_col] == gp1]
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
