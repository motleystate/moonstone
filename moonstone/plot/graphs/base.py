from abc import ABC, abstractmethod
from typing import Union, Tuple

import copy
import pandas as pd
import plotly.io
import plotly.express as px
import plotly.graph_objects as go
from plotly.validators.scatter.marker import SymbolValidator

import logging

from moonstone.utils.colors import generate_color_code

logger = logging.getLogger(__name__)


class BaseGraph(ABC):
    DEFAULT_COLOR = "#666666"

    def __init__(
        self,
        data: Union[pd.Series, pd.DataFrame],
    ):
        """
        :param data: data to plot
        """
        self.data = data

    @abstractmethod
    def plot_one_graph(self):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass

    def _color_scheme_species(self, colors: dict = None) -> dict:
        final_colors = {name: generate_color_code(name) for name in self.data.index}
        if colors is not None:
            final_colors.update(**colors)
        return final_colors

    def _color_scheme_metadata(
        self, metadata: Union[pd.DataFrame, pd.Series, list, int], colors: dict = None
    ) -> dict:
        """
        Generate a dictionary with the different groups as key and a color as value.
        If number of different groups > 38, then colors are repeated.

        Args:
            metadata: Can take DataFrame or Series, then all different values are listed and assigned a color.
                Or can take a list/array, then all elements of the list are assigned a color.
                Finally, metadata can be an integer which act like a list ranging from 0 to that integer (excluded),
                and each integer is assigned a color.
            colors: Dictionary that dictate a particular color to a particular group.
        """
        final_colors = {}
        if isinstance(metadata, pd.Series):
            all_gp = list(metadata.unique())
        elif isinstance(metadata, pd.DataFrame):
            all_gp = []
            for cc in metadata.columns:
                all_gp += list(metadata[cc].unique())
            all_gp = list(set(all_gp))  # doesn't remove all different nan
        # we can't do it manually because then some nan don't correspond to the one manually added
        elif type(metadata) is int:
            all_gp = list(range(0, metadata))
        else:  # array
            all_gp = metadata

        if len(all_gp) <= 10:
            final_colors = dict(zip(all_gp, px.colors.qualitative.Plotly))
        elif len(all_gp) <= 26:
            final_colors = dict(zip(all_gp, px.colors.qualitative.Alphabet))
        else:
            c = px.colors.qualitative.Alphabet + px.colors.qualitative.Set3
            final_colors = dict(zip(all_gp, c * (int(len(all_gp) / len(c)) + 1)))
        if colors is not None:
            final_colors.update(**colors)
        return final_colors

    def _valid_orientation_param(self, orientation: str) -> str:
        if orientation == "v" or orientation == "vertical":
            return "v"
        elif orientation == "h" or orientation == "horizontal":
            return "h"
        logger.warning("orientation=%s not valid, set to default (v).", orientation)
        return "v"

    def _handle_plotting_options_plotly(
        self, fig: go.Figure, plotting_options: dict
    ) -> go.Figure:
        """
        :param plotting_options: dictionary of dictionaries where the keys are the level to update
        ( "layout" | "traces" | "xaxes" | "yaxes" ...)

        Examples of keys that can be specified in layout dictionary :

        :param yaxis_type xaxis_type: Sets the axis type. By default, plotly attempts to determined the axis type
        by looking into the data of the traces that referenced the axis in question.
        ("-" | "linear" | "log" | "date" | "category" | "multicategory")
        :param shapes: list of dictionary for all the shapes to draw, with values for type,
                      x0, y0, x1, y1, xref [optional], yref [optional],
                      line dictionary [optional]
                      -> line dictionary contains lines descriptor like width, color,
                         or dash (dash, dashdot, dot etc.)
        :param title_text: title of the graph
        :param title_y title_x: coordinates of the title

        See `plotly documentation
        <https://plotly.com/python-api-reference/generated/plotly.graph_objects.Layout.html>`_.


        Examples of keys that can be specified in traces dictionary :
        :param marker_color

        Examples of keys that can be specified in xaxes dictionary or yaxes dictionary :
        :param tickangle: degree of orientation of tick labels
        :param title_text: label of the axis
        """
        for option in plotting_options.keys():
            updater = f"update_{option}"
            getattr(fig, updater)(plotting_options[option])
        return fig

    def _handle_output_plotly(self, fig: go.Figure, show: bool, output_file: str):
        if show is True:
            fig.show()

        if output_file:
            if output_file is True:
                # if no name given for the output file, a generic name is generated
                if self.data.name is not None:
                    output_file = (
                        self.data.name + "_" + self.__class__.__name__ + ".html"
                    )
                else:
                    output_file = self.__class__.__name__ + ".html"
            if output_file.split(".")[-1] == "html":
                plotly.io.write_html(fig, output_file)
            else:
                plotly.io.write_image(fig, output_file)


class GroupBaseGraph(BaseGraph):

    DEFAULT_COLORS = [
        "#A63A50",
        "#FFBF00",
        "#68ace8",
        "#97bf8f",
        "#28464B",
        "#6D5A72",
        "#FF8A5B",
        "#7C9EB2",
        "#F4F1DE",
        "#9CD08F",
    ]

    def __init__(self, *args, **kwargs):
        self.color_counter = 0
        super().__init__(*args, **kwargs)

    def _gen_default_color_dict(self, groups: list):
        return {
            groups[i]: self.DEFAULT_COLORS[i % len(self.DEFAULT_COLORS)]
            for i in range(0, len(groups))
        }

    def _get_group_color(self, group: str, group_color: dict):
        return group_color.get(group, self.DEFAULT_COLOR)

    @property
    def dictionary_marker_symbols(self):
        """
        Dictionary of all the symbols accepted by a scatter (or other) plot.
        key is the name of the symbol and value is its number
        ex: {'circle': 0, 'circle-open': 100, ...}
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
        If number of different groups > 162, then symbols are repeated.

        Args:
            metadata: Can take DataFrame or Series, then all different values are listed and assigned a color.
                Or can take a list/array, then all elements of the list are assigned a color.
                Finally, metadata can be an integer which act like a list ranging from 0 to that integer (excluded),
                and each integer is assigned a color.
            symbols: Dictionary that dictate a particular symbol to a particular group.
        """
        if symbols is not None and len(set(groups).difference(symbols.keys())) == 0:
            return symbols

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

    def _symbol_scheme3d(self, groups: list, symbols: dict = None):
        """
        Generate a dictionary with the different groups as key and the name of a symbol as value.
        If number of different groups > 8, then symbols are repeated.
        NB: Unlike for the 2D plots that accept 162 different symbols, 3D plots accept only 8 different symbols
        and only their name (and not their number).

        Args:
            metadata: Can take DataFrame or Series, then all different values are listed and assigned a color.
                Or can take a list/array, then all elements of the list are assigned a color.
                Finally, metadata can be an integer which act like a list ranging from 0 to that integer (excluded),
                and each integer is assigned a color.
            symbols: Dictionary that dictate a particular symbol to a particular group.
        """
        if symbols is not None and len(set(groups).difference(symbols.keys())) == 0:
            return symbols

        accepted_symbols = [
            'circle', 'cross', 'diamond', 'square',
            'circle-open', 'diamond-open', 'square-open',
            'x'  # this one is crazy big so let's put it as the last resort
        ]
        if len(groups) > 8:
            dic_symbols = dict(zip(groups, accepted_symbols * (int(len(groups) / len(accepted_symbols)) + 1)))
        else:
            if symbols is not None:
                # for reusability, we can't use set().difference()
                accepted_symbols = [x for x in accepted_symbols if x not in symbols.values()]
                skeys = set(symbols.keys())
                groups = [x for x in groups if x not in skeys]
            dic_symbols = dict(zip(groups, accepted_symbols))

        if symbols is not None:
            dic_symbols.update(**symbols)

        return dic_symbols

    @abstractmethod
    def _gen_fig_trace(self, x: list, y: list, name: str, text: list, color: str):
        pass

    def _gen_oriented_fig_trace(
        self,
        fig: go.Figure,
        axis1: pd.Series,
        axis2: pd.Series,
        name: str,
        text: str,
        color: str,
        orientation: str,
        **kwargs,
    ) -> go.Figure:
        if orientation == "h":
            fig.add_trace(
                self._gen_fig_trace(
                    axis2,
                    axis1,
                    name,
                    text,
                    color,
                    **kwargs,
                )
            )
        else:  # default "v"
            fig.add_trace(
                self._gen_fig_trace(
                    axis1,
                    axis2,
                    name,
                    text,
                    color,
                    **kwargs,
                )
            )
        return fig

    def _prep_for_plot_one_graph(
        self,
        group_col: str,
        groups: list,
        colors: dict,
        show_counts: bool,
        sort_groups: bool,
    ) -> Tuple[list, dict, dict]:
        if groups is None:
            groups = list(self.data[group_col].unique())
        if sort_groups:
            groups.sort()
        if colors is None:
            colors = self._gen_default_color_dict(groups)

        if show_counts:
            counts = self.data[group_col].value_counts().to_dict()
            names = {group: f"{group} (n={counts[group]})" for group in groups}
        else:
            names = {group: group for group in groups}
        return groups, colors, names

    def plot_one_graph(
        self,
        data_col: str,
        group_col: str,
        group_col2: str = None,
        groups: list = None,
        groups2: list = None,
        sort_groups: bool = False,
        colors: dict = None,
        plotting_options: dict = None,
        orientation: str = "v",
        show_counts: bool = False,
        show: bool = True,
        output_file: Union[bool, str] = False,
        **kwargs,
    ) -> go.Figure:
        """
        :param data_col: column with data to visualize
        :param group_col: column used to group data
        :param group_col2: (optional) second column used to group data
        :param groups: specifically select groups to display among group_col
        :param groups2: specifically select groups to display among group_col2
        :param sort_groups: whether to sort groups to display or not
        :param colors: dictionnary with group_col2 (or group_col if no group_col2) values as keys
        :param orientation: orientation of the graph. {"v" (or "vertical")(default), "h" (or "horizontal")}
        and their associated color as values
        """

        orientation = self._valid_orientation_param(orientation)

        fig = go.Figure()

        if group_col2:
            groups2, colors, names = self._prep_for_plot_one_graph(
                group_col2, groups2, colors, show_counts, sort_groups
            )

            if sort_groups and not groups:
                groups = list(self.data[group_col].unique())
                groups.sort()

            if groups:
                filtered_df = self.data[self.data[group_col].isin(groups)]
                filtered_df[group_col] = filtered_df[group_col].astype("category")
                filtered_df[group_col].cat = filtered_df[group_col].cat.set_categories(groups)
                filtered_df = filtered_df.sort_values([group_col])
            else:
                filtered_df = copy.deepcopy(self.data)

            for group in groups2:
                filtered_df2 = filtered_df[filtered_df[group_col2] == group]
                fig = self._gen_oriented_fig_trace(
                    fig,
                    filtered_df2[group_col],
                    filtered_df2[data_col],
                    str(names[group]),
                    filtered_df.index,
                    self._get_group_color(group, colors),
                    orientation,
                    **kwargs,
                )

            fig.update_layout(boxmode="group", legend_title_text=group_col2)
        else:
            groups, colors, names = self._prep_for_plot_one_graph(
                group_col, groups, colors, show_counts, sort_groups
            )
            for group in groups:
                filtered_df = self.data[self.data[group_col] == group]
                fig = self._gen_oriented_fig_trace(
                    fig,
                    filtered_df[group_col],
                    filtered_df[data_col],
                    str(names[group]),
                    filtered_df.index,
                    self._get_group_color(group, colors),
                    orientation,
                    **kwargs,
                )

        if orientation == "h":
            fig.update_traces(orientation="h")

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig
