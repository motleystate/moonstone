import logging
from typing import Union

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats
import statsmodels.api as sm


from moonstone.plot.graphs.base import BaseGraph
from moonstone.utils.dict_operations import merge_dict

logger = logging.getLogger(__name__)


class ScatterTrendlines(BaseGraph):

    def _generate_default_plotting_options(self, log_x, log_y):
        if log_x:
            plotting_options = {"xaxes": {"type": "log"}}
        else:
            plotting_options = {}
        if log_y:
            plotting_options = merge_dict(
                plotting_options, {"yaxes": {"type": "log"}}
            )
        return plotting_options

    def _valid_outliers_param(self, outliers):
        if outliers == 'trim' or outliers == 'trimmed':
            return 'trim'
        elif outliers != 'keep' and outliers != 'kept':
            logger.warning("outliers='%s' not valid, set to default (keep).", outliers)
        return 'keep'

    def _color_palette(self, n):
        if n <= 10:
            return px.colors.qualitative.Plotly[:n]
        elif n <= 26:
            return px.colors.qualitative.Alphabet[:n]
        else:
            return px.colors.qualitative.Alphabet * int(n/26) + px.colors.qualitative.Alphabet[:n % 26]

    def Scatter_plot(
        self,
        x_ser: pd.Series, y_ser: pd.Series,
        trendlines: list,
        color_ser: pd.Series = None,
        forced_color_dic: dict = None,
        show: bool = True, output_file: Union[bool, str] = False,
        plotting_options: dict = {},
    ):
        """
        To plot
        Args:
            - x_ser: pandas Series containing the data to plot along the x axis.
            - y_ser: pandas Series containing the data to plot along the y axis.
            - forced_color_dic: directory
            - output_file: str = None,
            - plotting_options: dict = {},
            - trendlines: list of tuples describing the trendlines.
              Tuples contain a dataframe of x and y coordinates tracing the line, and the number of the group
            - show: {True (default), False}. Set to False if you do not wish to generate the graph.
        """
        if color_ser is not None and not forced_color_dic:
            u = list(color_ser.value_counts().index)
            u.sort()
            color_dic = {}
            if type(u[0]) != float and type(u[0]) != int:
                for i in range(0, len(u)):
                    color_dic[u[i]] = i
                if color_ser.dropna().shape[0] != color_ser.shape[0]:
                    color_dic[pd.NA] = i+1
                color_ser.replace(color_dic, inplace=True)

        fig = go.Figure(data=go.Scatter(
            x=x_ser, y=y_ser,
            mode='markers',
            marker_color=color_ser,
            text=x_ser.index            # hover text goes here
        ))

        for tl in trendlines:
            if forced_color_dic:
                marker_color = {'color': forced_color_dic[tl[1]]}
            else:
                marker_color = None
            fig.add_trace(go.Scatter(
                x=tl[0].x, y=tl[0].y,
                mode='lines',
                marker=marker_color
            ))

        title = f'Correlation between {y_ser.name} and {x_ser.name}'
        default_plotting_options = {
            "layout": {
                "title": title,
                "xaxis_title": x_ser.name,
                "yaxis_title": y_ser.name,
            }}
        plotting_options = merge_dict(
            plotting_options, default_plotting_options
        )

        self._handle_plotting_options_plotly(fig, plotting_options)
        self._handle_output_plotly(fig, show, output_file)

        return fig

    def _compute_trendline(
        self,
        X: pd.Series, Y: pd.Series,
        log_x: bool = False, log_y: bool = False,
        y_trendline: bool = False
    ):
        """
        Args:
            - X: series of x coordinates.
            - Y: series of y coordinates.
            - log_x, log_y: set to True if the x (or y) axis is logarithmic.
            - y_trendline: set to True to generate the array of y values when 'projecting' x on the trendline
        """
        if log_y:
            if np.any(Y <= 0):
                raise ValueError(
                    "Can't do OLS trendline with `log_y=True` when `y` contains non-positive values."
                )
            Y = np.log10(Y)

        if log_x:
            if np.any(X <= 0):
                raise ValueError(
                    "Can't do OLS trendline with `log_x=True` when `x`  contains non-positive values."
                )
            X = np.log10(X)

        X = sm.add_constant(X)
        model = sm.OLS(Y, X, missing="drop")
        results = model.fit()
        if y_trendline:
            y_tl = results.predict()

            if log_y:
                y_tl = np.power(10, y_tl)
            return results.params, y_tl

        else:
            return results.params

    def _distance_point_to_line(self, x, y, a, cst):
        d = abs(-a*x+y-cst)/np.sqrt(np.power(-a, 2)+1)
        return d

    # NB: should work if both axis are linear or both axis are log
    # -> to think about cases where one is linear and the other log
    def _assign_sample_to_a_linegroup(self, x, y, list_params):
        best_d = np.inf
        params_best_d = None
        for params in list_params:
            d = self._distance_point_to_line(x, y, params[1], params[0])
            if d < best_d:
                best_d = d
                params_best_d = params
        return params_best_d, best_d

    def _assign_samples_to_linegroups(
        self,
        dataframe: pd.DataFrame, x_col: str, y_col: str,
        list_params: list, log_x: bool, log_y: bool
    ):
        if log_x:
            X = np.log10(dataframe[x_col])
        else:
            X = dataframe[x_col]

        if log_y:
            Y = np.log10(dataframe[y_col])
        else:
            Y = dataframe[y_col]

        dic_l = {}
        dic_d = {}
        for s in dataframe.index:
            lp, d = self._assign_sample_to_a_linegroup(X.loc[s], Y.loc[s], list_params)
            dic_l[s] = lp[2]
            dic_d[s] = d
        return pd.Series(dic_l, name="group"), pd.DataFrame.from_dict(dic_d, orient="index", columns=["distance"])

    def define_n_trendlines(
        self,
        x_col: str, y_col: str, n: int,
        log_x: bool = False, log_y: bool = False,
        outliers: str = "keep", nb_iter: int = 100,
        random_state: int = None,
        show: str = "last", output_file: Union[bool, str] = False,
        v: bool = True,
    ) -> dict:
        """
        Args:
            - x_col: column of the dataframe to plot on the x axis of the graph.
            - y_col: column of the dataframe to plot on the y axis of the graph.
            - n: number of different trendlines to trace. >1 (otherwise use the plotly function).
            - log_x, log_y: set to True if the x (or y) axis is logarithmic.
            - outliers: {'keep' (default), 'trim'}. Set to 'trim' to remove outliers (|z-score| > 3) between iteration.
            - nb_iter: maximum number of iteration to run the research of n stable trendlines.
              NB: If stability is reached before, then function ends before this maximum number.
            - show: {'last' (default), 'all', False}. Whether to show the graph for every iterations ('all')
              or just the graph corresponding to the 'last' iteration or not (False).
            - output_file: file path to output the generated graph.
            - v: verbosity.
            - random_state: Set to an integer for reproductibility purpose.
        """
        to_return = {}
        try:
            self.outliers
        except AttributeError:
            self.outliers = self._valid_outliers_param(outliers)

        plotting_options = self._generate_default_plotting_options(log_x, log_y)
        forced_color_dic = dict(zip(range(0, n), self._color_palette(n)))

        # random first split of the samples into n groups
        dataframe = self.data.sample(frac=1, random_state=random_state)  # shuffles dataframe's rows
        chunks_list = np.array_split(dataframe, n)
        group_ser_prec = pd.Series()
        for c in range(len(chunks_list)):
            group_ser_prec = pd.concat(
                [group_ser_prec, pd.Series([c]*len(chunks_list[c].index), index=chunks_list[c].index)]
            )
        i = 1
        stable = False
        while i <= nb_iter:
            list_params = []
            trendlines = []
            for c in range(len(chunks_list)):
                if chunks_list[c].shape[0] >= 2:        # cannot compute a line from less than 2 points
                    p, y_tl = self._compute_trendline(
                        chunks_list[c][x_col], chunks_list[c][y_col], log_x=log_x, log_y=log_y, y_trendline=True
                    )  # if outliers == 'trim', outliers not used to compute trendline... [1/2]
                    list_params += [(p[0], p[1], c)]
                    trendlines += [
                        (pd.DataFrame({"x": chunks_list[c][x_col], "y": y_tl}), c)
                        # x = same x as dot/sample, y = projection of x on trendline
                    ]
            group_ser, distance_df = self._assign_samples_to_linegroups(
                dataframe, x_col, y_col, list_params, log_x, log_y
            )
            # ... but outliers still put inside a group and still get a distance to the trendline [2/2]

            if group_ser.equals(group_ser_prec):
                if v:
                    print(f"Stable after {i} iterations.")
                stable = True
                break

            if show == 'all':
                self.Scatter_plot(
                    dataframe[x_col], dataframe[y_col],
                    color_ser=group_ser.apply(lambda x: forced_color_dic[x]),
                    forced_color_dic=forced_color_dic,
                    show=True, output_file=None,
                    plotting_options=plotting_options,
                    trendlines=trendlines
                )

            chunks_list = []
            if self.outliers == 'trim':
                distance_df['z_score'] = stats.zscore(distance_df['distance'])

                # 1) outliers defined as the data point with a distance to trendline diverging from the mean distance
                #    by more than 3 standard deviations
                to_return['list_outliers'] = distance_df.loc[distance_df['z_score'].abs() > 3].index
                # 2) we only keep non outliers to create the new chunks_list to compute the new trendlines
                smaller_group_ser = group_ser.loc[~group_ser.index.isin(to_return['list_outliers'])]
                for j in group_ser.unique():
                    chunks_list += [dataframe.loc[smaller_group_ser[smaller_group_ser == j].index]]
            else:
                for j in group_ser.unique():
                    chunks_list += [dataframe[group_ser == j]]

            group_ser_prec = group_ser
            i += 1

        if not stable:
            logger.warning(f"Not stable after {i-1} iterations.")

        to_return["figure"] = self.Scatter_plot(
            dataframe[x_col], dataframe[y_col],
            color_ser=group_ser.apply(lambda x: forced_color_dic[x]),
            forced_color_dic=forced_color_dic,
            show=bool(show), output_file=output_file,
            plotting_options=plotting_options,
            trendlines=trendlines
        )

        # outputs: dictionary containing
        # 1) 'trendlines' = list of n tuples containing:
        #    a) DataFrame with samples as rows, and x (= sample's x) and y (= projection of x on trendline) as columns
        #    b) group number
        # 2) 'group_series' = Series with samples as keys and group number of sample as values
        # 3) 'distance_dataframe' = DataFrame with samples as rows, and distance (of sample to projection on trendline)
        #    and z-score as columns
        # 4) 'figure' = scatter plot with trendlines
        # 5) (optional) 'list_outliers' = list of the datapoint that were considered outliers in the last iteration

        # [self.data.index] and .loc[self.data.index] -> to put samples back in the original order
        to_return["trendlines"] = trendlines
        to_return["group_series"] = group_ser[self.data.index]
        to_return["distance_dataframe"] = distance_df.loc[self.data.index]
        return to_return

    def bootstraps_define_n_trendlines(
        self,
        x_col: str, y_col: str, n: int,
        log_x: bool = False, log_y: bool = False,
        outliers: str = "keep",
        nb_iter: int = 100, nb_bootstraps: int = 25,
        random_state: int = None,
        show: bool = True, output_file: Union[bool, str] = False
    ) -> dict:
        """
        Since there's a part of randomness in generating the first separation of points into n parts, more than one
        stable state can be achieved (local min). So function called many time to find the state with the lowest mean
        distance between points and their trendline.

        Args:
            - x_col: column of the dataframe to plot on the x axis of the graph.
            - y_col: column of the dataframe to plot on the y axis of the graph.
            - n: number of different trendlines to trace. >1 (otherwise use the plotly function).
            - log_x, log_y: set to True if the x (or y) axis is logarithmic.
            - outliers: {'keep' (default), 'trim'}. Set to 'trim' to remove outliers (|z-score| > 3) between iteration.
            - nb_iter: maximum number of iteration to run the research of n stable trendlines.
              NB: If stability is reached before, then function ends before this maximum number.
            - nb_bootstraps: number of times define_n_trendlines is called. The best results, that is the ones were the
              mean distance of data point to their associated trendline is the smallest, are returned.
            - show: {True (default), False}. Set to False if you do not wish to generate the graph.
            - output_file: file path to output the generated graph.
            - random_state: Set to an integer for reproductibility purpose; increased by 1 for each iteration.
        """
        self.outliers = self._valid_outliers_param(outliers)

        best_mean = np.inf
        for i in range(nb_bootstraps):
            dfnt_output = self.define_n_trendlines(
                x_col, y_col, n,
                log_x=log_x, log_y=log_y,
                outliers=self.outliers, nb_iter=nb_iter,
                show=False, output_file=False, v=False,
                random_state=random_state,
            )
            if dfnt_output["distance_dataframe"]['distance'].mean() < best_mean:
                best_mean = dfnt_output["distance_dataframe"]['distance'].mean()
                best_dfnt_output = dfnt_output
                best_rs = random_state
            if random_state:
                random_state += 1

        if show or output_file:
            self._handle_output_plotly(best_dfnt_output["figure"], show, output_file)
        print(best_rs)
        # output: cf define_n_trendlines' output
        return best_dfnt_output

    def plot_one_graph(
        self,
        x_col: str, y_col: str, n: int,
        log_x: bool = False, log_y: bool = False,
        outliers: str = "keep",
        nb_iter: int = 100, nb_bootstraps: int = 1,
        random_state: int = None,
        show: bool = True, output_file: Union[bool, str] = False
    ) -> go.Figure:
        """
        Args:
            - x_col: column of the dataframe to plot on the x axis of the graph.
            - y_col: column of the dataframe to plot on the y axis of the graph.
            - n: number of different trendlines to trace. >1 (otherwise use the plotly function).
            - log_x, log_y: set to True if the x (or y) axis is logarithmic.
            - outliers: {'keep' (default), 'trim'}. Set to 'trim' to remove outliers (|z-score| > 3) between iteration.
            - nb_iter: maximum number of iteration to run the research of n stable trendlines.
              NB: If stability is reached before, then function ends before this maximum number.
            - nb_bootstraps: number of times define_n_trendlines is called. The best results, that is the ones were the
              mean distance of data point to their associated trendline is the smallest, are returned.
            - show: {True (default), False}. Set to False if you do not wish to generate the graph.
            - output_file: file path to output the generated graph.
            - random_state: Set to an integer for reproductibility purpose.
        """
        if nb_bootstraps > 1:
            fig = self.bootstraps_define_n_trendlines(
                x_col, y_col, n,
                log_x=log_x, log_y=log_y,
                outliers=outliers,
                nb_iter=nb_iter, nb_bootstraps=nb_bootstraps,
                random_state=random_state,
                show=show, output_file=output_file,
            )["figure"]
        else:
            if show:
                show = 'last'

            fig = self.define_n_trendlines(
                x_col, y_col, n,
                log_x=log_x, log_y=log_y,
                outliers=outliers, nb_iter=nb_iter,
                random_state=random_state,
                show=show, output_file=output_file, v=False
            )["figure"]

        return fig
