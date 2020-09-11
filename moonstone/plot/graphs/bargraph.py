import plotly.graph_objects as go
import plotly.io

from moonstone.plot.graphs.base import BaseGraph


class BarGraph(BaseGraph):

    def _set_plotting_options(self, fig, plotting_options):
        if 'log' in self.plotting_options.keys() and self.plotting_options['log']:
            fig.update_layout(yaxis_type="log")
        if 'tickangle' in self.plotting_options.keys():
            fig.update_xaxes(tickangle=self.plotting_options['tickangle'])
        if 'colorbar' in plotting_options.keys():
            fig.update_traces(marker_color=plotting_options['colorbar'])
        return fig

    def _set_labels(self, fig, title, xlabel, ylabel):
        if title is not None:
            fig.update_layout(title_text=title, title_x=0.5)
        if xlabel is not None:
            fig.update_xaxes(title_text=xlabel)
        if ylabel is not None:
            fig.update_yaxes(title_text=ylabel)
        return fig

    def _get_chart(self, orientation: str = "v", ascending: bool = None) -> go.Bar:
        if ascending is not None:
            data = self.data.sort_values(ascending=ascending)
        else:
            data = self.data
        x = list(data.index)
        y = list(data)
        if orientation == "v":
            return go.Bar(x=x, y=y, orientation=orientation)
        return go.Bar(x=y, y=x, orientation=orientation)

    def plot_one_graph(self, orientation: str = "v", ascending: bool = None,
                       title: str = None, xlabel: str = None, ylabel: str = None,
                       plotting_options: dict = None):
        """
        :param title: title of the graph
        :param xlabel: label of the x axis
        :param ylabel: label of the y axis
        :param reset_xnames_dic: to rename the names of the values in the x axis. \n
                                 Example for a plot of the distribution of smoking habits :
                                 reset_xnames_dic={'y': 'smoker', 'n': 'non smoker'}
        """
        fig = go.Figure(self._get_chart(orientation=orientation, ascending=ascending))

        if plotting_options is not None:
            fig = self._set_plotting_options(fig, plotting_options)

        fig = self._set_labels(fig, title, xlabel, ylabel)

        if self.show is True:
            fig.show()

        if self.output_file:
            plotly.io.write_html(fig, self.output_file)
