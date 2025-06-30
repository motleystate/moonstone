import itertools
import logging
import numpy as np
from numpy.linalg import eig
import plotly.graph_objects as go
from scipy.stats import chi2
from typing import Union

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

    def _confidence_ellipse(
        self, x, y, q_confidence: float = 0.95, n_points: int = 100
    ):
        """
        Args:
            x: array of the x data.
            y: array of the y data.
            q_confidence: proportion (between 0 and 1) of confidence of the desired ellipses.
            n_points: number of points to trace the ellipse. It will change the smoothness of the curve,
              but not it's shape or size.
        """
        cov = np.cov(x, y)
        mean = np.mean(x), np.mean(y)
        vals, vecs = eig(cov)
        order = vals.argsort()[::-1]
        vals, vecs = vals[order], vecs[:, order]

        theta = np.linspace(0, 2*np.pi, n_points)
        unit_circle = np.array([np.cos(theta), np.sin(theta)])
        scale = np.sqrt(chi2.ppf(q_confidence, df=2))
        ellipse = vecs @ np.diag(np.sqrt(vals)) @ unit_circle * scale
        ellipse[0, :] += mean[0]
        ellipse[1, :] += mean[1]
        return ellipse

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
        symbols: dict = None,
        groups: list = None,
        groups2: list = None,
        q_confel: float = 0,
        **kwargs
    ) -> go.Figure:
        """
        Args:
            first_col: col name for x data.
            second_col: col name for y data.
            group_col: column used to group data (differentiated by color).
            group_col2: second columns used to group data (differentiated by symbols).
            colors: dictionary to impose specific colors to one or more groups of group_col.
            symbols: dictionary to impose specific symbols to one or more groups of group_col2.
            groups: select specific groups to display among group_col.
            groups2: select specific groups to display among group_col2.
            q_confel: proportion (between 0 and 1) of confidence intervals to represent with ellipses on the graph.
              NB: 0 = no ellipse.
        """
        if groups is None:
            groups = list(self.data[group_col].unique())
        if colors is None:
            colors = self._gen_default_color_dict(groups)
        if group_col2 and groups2 is None:
            groups2 = list(self.data[group_col2].unique())

        fig = go.Figure()

        if group_col2:
            dic_symbol = self._symbol_scheme(groups2, symbols)
            for gp2, linetype in zip(
                groups2, itertools.cycle(["dot", "dash", "dashdot", "longdash", "longdashdot", "solid"])
            ):
                gp2_filtered_df = self.data[self.data[group_col2] == gp2]
                for gp1 in groups:
                    filtered_df = gp2_filtered_df[gp2_filtered_df[group_col] == gp1]
                    if not filtered_df.empty:
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
                        if bool(q_confel):
                            if filtered_df.shape[0] < 2:
                                logger.warning(f"Cannot compute an ellipse of confidence for group {gp1} - {gp2}, \
because only one sample in that group")
                            else:
                                # computing the data to trace the ellipse of confidence
                                ell = self._confidence_ellipse(
                                    filtered_df[first_col], filtered_df[second_col],
                                    q_confidence=q_confel, n_points=kwargs.get("n_points", 100)
                                )

                                # adding the ellipse of confidence to the figure
                                fig.add_trace(go.Scatter(
                                    x=ell[0], y=ell[1],
                                    mode='lines',
                                    name=f'Ellipse {gp1} - {gp2}',
                                    line=dict(color=self._get_group_color(gp1, colors), dash=linetype),
                                    showlegend=False,
                                    hovertext=f'Ellipse {gp1} - {gp2} [{int(q_confel*100)}% confidence]',
                                    hoverinfo="text",
                                ))
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
                if bool(q_confel):
                    if filtered_df.shape[0] < 2:
                        logger.warning(f"Cannot compute an ellipse of confidence for group {group}, \
because only one sample in that group")
                    else:
                        ell = self._confidence_ellipse(
                            filtered_df[first_col], filtered_df[second_col],
                            q_confidence=q_confel, n_points=kwargs.get("n_points", 100)
                        )
                        fig.add_trace(go.Scatter(
                            x=ell[0], y=ell[1],
                            mode='lines',
                            name=f'Ellipse {group}',
                            line=dict(color=self._get_group_color(group, colors), dash="dot"),
                            showlegend=False,
                            hovertext=f'Ellipse {group} [{int(q_confel*100)}% confidence]',
                            hoverinfo="text",
                        ))

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
        marker: dict = None
    ) -> go.Scatter3d:
        return go.Scatter3d(
            x=x,
            y=y,
            z=z,
            name=name,
            text=text,
            marker_color=color,
            mode="markers",
            marker=marker
        )

    def plot_one_graph(
        self,
        first_col: str,
        second_col: str,
        third_col: str,
        group_col: str,
        group_col2: str = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
        colors: dict = None,
        symbols: dict = None,
        groups: list = None,
        groups2: list = None,
        plotting_options: dict = None,
        **kwargs
    ) -> go.Figure:
        """
        Args:
            first_col: col name for x data
            second_col: col name for y data
            third_col: col name for z data
            group_col: column used to group data (differentiated by color).
            group_col2: second columns used to group data (differentiated by symbols).
            colors: dictionary to impose specific colors to one or more groups of group_col.
            symbols: dictionary to impose specific symbols to one or more groups of group_col2.
            groups: select specific groups to display among group_col.
            groups2: select specific groups to display among group_col2.
        """
        kwargs.pop("q_confel", None)
        # @TODO: eventually try to compute ellipse of confidence for 3D plot

        if groups is None:
            groups = list(self.data[group_col].unique())
        if colors is None:
            colors = self._gen_default_color_dict(groups)
        if group_col2 and groups2 is None:
            groups2 = list(self.data[group_col2].unique())

        fig = go.Figure()

        if group_col2:
            dic_symbol = self._symbol_scheme3d(groups2, symbols)
            for gp2 in groups2:
                gp2_filtered_df = self.data[self.data[group_col2] == gp2]
                for gp1 in groups:
                    filtered_df = gp2_filtered_df[gp2_filtered_df[group_col] == gp1]
                    if not filtered_df.empty:
                        fig.add_trace(
                            self._gen_fig_trace(
                                filtered_df[first_col],
                                filtered_df[second_col],
                                filtered_df[third_col],
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
