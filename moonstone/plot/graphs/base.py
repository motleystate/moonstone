from abc import ABC, abstractmethod
from typing import Union

import copy
import pandas as pd
import plotly.io
import plotly.graph_objects as go


class BaseGraph(ABC):
    DEFAULT_COLOR = "#666666"

    def __init__(self, data: Union[pd.Series, pd.DataFrame], plotting_options: dict = None,
                 show: bool = True, output_file: Union[bool, str] = False):
        """
        :param data: data to plot
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        self.data = data

    @abstractmethod
    def plot_one_graph(self):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass

    def _valid_orientation_param(self, orientation):
        if orientation == "v" or orientation == "vertical":
            return "v"
        elif orientation == "h" or orientation == "horizontal":
            return "h"
        logger.warning("orientation=%s not valid, set to default (v)", orientation)
        return "v"

    def _handle_plotting_options_plotly(self, fig, plotting_options: dict):
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

    def _handle_output_plotly(self, fig, show: bool, output_file: str):
        if show is True:
            fig.show()

        if output_file:
            if output_file is True:
                # if no name given for the output file, a generic name is generated
                if self.data.name is not None:
                    output_file = self.data.name+"_"+self.__class__.__name__+".html"
                else:
                    output_file = self.__class__.__name__+".html"
            if output_file.split(".")[-1] == "html":
                plotly.io.write_html(fig, output_file)
            else:
                plotly.io.write_image(fig, output_file)


class GroupBaseGraph(BaseGraph):

    DEFAULT_COLORS = [
        "#A63A50", "#FFBF00", "#68ace8", "#97bf8f", "#28464B",
        "#6D5A72", "#FF8A5B", "#7C9EB2", "#F4F1DE", "#9CD08F",
    ]

    def __init__(self, *args, **kwargs):
        self.color_counter = 0
        super().__init__(*args, **kwargs)

    def _gen_default_color_dict(self, groups: list):
        return {groups[i]: self.DEFAULT_COLORS[i % len(self.DEFAULT_COLORS)] for i in range(0, len(groups))}

    def _get_group_color(self, group: str, group_color: dict):
        return group_color.get(group, self.DEFAULT_COLOR)

    @abstractmethod
    def _gen_fig_trace(self, x: list, y: list, name: str, text: list, color: str):
        pass

    def _gen_fig_traces(
        self, fig, dataframe: pd.DataFrame, 
        data_col: str, group_col: str, groups: list,
        names: dict, orientation: str, colors: dict,
        **kwargs
        ):
        for group in groups:
            filtered_df = dataframe[dataframe[group_col] == group]
            if orientation == "h":
                fig.add_trace(self._gen_fig_trace(
                    filtered_df[data_col], filtered_df[group_col],
                    str(names[group]), filtered_df.index, self._get_group_color(group, colors),
                    **kwargs,
                ))
            else:    # default "v"
                fig.add_trace(self._gen_fig_trace(
                    filtered_df[group_col], filtered_df[data_col],
                    str(names[group]), filtered_df.index, self._get_group_color(group, colors),
                    **kwargs,
                ))
        return fig
        

    def _prep_for_plot_one_graph(
        self, group_col: str, groups: list, colors: dict,
        show_counts: bool, sort_groups: bool
    ):
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
        self, data_col: str, group_col: str, group_col2: str = None,
        groups: list = None, groups2: list = None,
        sort_groups: bool = False, colors: dict = None,
        plotting_options: dict = None, orientation: str = "v",
        show_counts: bool = False,
        show: bool = True, output_file: Union[bool, str] = False,
        **kwargs
    ):
        """
        :param data_col: column with data to visualize
        :param group_col: column used to group data
        :param colors: dictionnary with group_col2 (or group_col if no group_col2) values as keys
        and their associated color as values
        """

        orientation=self._valid_orientation_param(orientation)
        
        fig = go.Figure()

        if group_col2:
            groups2, colors, names = self._prep_for_plot_one_graph(
                group_col2, groups2, colors, show_counts, sort_groups
            )
            if sort_groups and groups:
                groups.sort()

            fig = go.Figure()

            if groups:
                filtered_df = self.data[self.data[group_col].isin(groups)]
                filtered_df[group_col] = filtered_df[group_col].astype("category")
                filtered_df[group_col].cat.set_categories(groups, inplace=True)
                filtered_df = filtered_df.sort_values([group_col])
            else:
                filtered_df = copy.deepcopy(self.data)

            fig = self._gen_fig_traces(
                fig, filtered_df, data_col, group_col2, groups2, names, orientation, colors, 
                **kwargs
                )
            for group in groups2:
                filtered_df2 = filtered_df[filtered_df[group_col2] == group]

                if orientation == "h":
                    fig.add_trace(self._gen_fig_trace(
                        filtered_df2[data_col], filtered_df2[group_col],
                        str(names[group]), filtered_df.index, self._get_group_color(group, colors),
                        **kwargs,
                    ))
                else:    # default "v"
                    fig.add_trace(self._gen_fig_trace(
                        filtered_df2[group_col], filtered_df2[data_col],
                        str(names[group]), filtered_df.index, self._get_group_color(group, colors),
                        **kwargs,
                    ))
            fig.update_layout(
                boxmode='group', legend_title_text=group_col2
            )
        else:
            groups, colors, names = self._prep_for_plot_one_graph(
                group_col, groups, colors, show_counts, sort_groups
            )
            fig = self._gen_fig_traces(
                fig, self.data, data_col, group_col, groups, names, orientation, colors,
                **kwargs
                )

        if orientation == "h":
            fig.update_traces(orientation="h")

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)
