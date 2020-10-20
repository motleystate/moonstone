from abc import ABC, abstractmethod
from typing import Union

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

    def _handle_output_plotly(self, fig, show: bool, output_file: str, log_scale: bool = False):
        if log_scale:
            fig.update_yaxes(type="log", title=f"{fig.layout.yaxis.title.text} (log)")

        if show is True:
            fig.show()

        if output_file:
            if output_file is True:
                # if no name given for the output file, a generic name is generated
                if self.data.name is not None:
                    output_file = self.data.name+"_"+self.__class__.__name__+".html"
                else:
                    output_file = self.__class__.__name__+".html"
            plotly.io.write_html(fig, output_file)


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

    def plot_one_graph(
        self, data_col: str, group_col: str, plotting_options: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        log_scale: bool = False, colors: dict = None, sort_groups: bool = False,
        groups: list = None,
    ):
        """
        :param data_col: column with data to visualize
        :param group_col: column used to group data
        """
        if groups is None:
            groups = list(self.data[group_col].unique())
        if sort_groups:
            groups.sort()
        if colors is None:
            colors = self._gen_default_color_dict(groups)
        fig = go.Figure()

        for group in groups:
            filtered_df = self.data[self.data[group_col] == group]
            fig.add_trace(self._gen_fig_trace(
                filtered_df[group_col], filtered_df[data_col],
                str(group), filtered_df.index, self._get_group_color(group, colors)
            ))

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file, log_scale)
