from typing import Dict, Optional, Union

import pandas as pd
import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph
from moonstone.plot.graphs.category_flows import prepare_categories


class SankeyDiagram(BaseGraph):
    """
    Sankey diagram: one column of boxes per categorical variable, with ribbons between neighbouring columns whose
    width is the number of samples (or the weight) going from one category to the next.

    Ribbons only join neighbouring columns, so a sample cannot be followed across three columns or more; use
    :class:`ParallelCategoriesGraph` for that.
    """

    def __init__(self, data: pd.DataFrame, weights: Optional[pd.Series] = None):
        """
        :param data: one row per sample, one column per categorical variable, plotted left to right
        :param weights: weight of each row (default 1 per row), with the same index as data
        """
        super().__init__(data)
        self.weights = weights

    def plot_one_graph(
        self,
        order: Optional[Dict[str, list]] = None,
        column_titles: Optional[Dict[str, str]] = None,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
    ) -> go.Figure:
        """
        :param order: for some columns, the order wanted for the categories, from top to bottom. By default, from the
            largest to the smallest.
        :param column_titles: for some columns, the title to display above the column instead of the column name
        :param plotting_options: options of plotting that will override the default setup
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file (.html or an image format)
        """
        str_df, weights, cat_order, titles = prepare_categories(self.data, self.weights, order, column_titles)
        columns = list(str_df.columns)
        n_cols = len(columns)

        # one node per (column, category), numbered column by column in the order wanted
        node_id = {}
        node_label, node_x, node_y, node_total = [], [], [], []
        for i, col in enumerate(columns):
            totals = weights.groupby(str_df[col]).sum()
            gap = 0.03
            usable = 1 - gap * (len(totals) - 1)
            position = 0.0
            for cat in cat_order[col]:
                share = totals[cat] / totals.sum() * usable if totals.sum() > 0 else 0.0
                node_id[(col, cat)] = len(node_label)
                node_label.append(f"{cat}<br>{totals[cat]:.4g}")
                node_x.append(0.001 + 0.998 * i / (n_cols - 1))
                node_y.append(min(max(position + share / 2, 0.001), 0.999))
                node_total.append(totals[cat])
                position += share + gap

        source, target, values = [], [], []
        for left, right in zip(columns[:-1], columns[1:]):
            flows = weights.groupby([str_df[left], str_df[right]]).sum()
            for (cat_left, cat_right), value in flows.items():
                if value == 0:
                    continue
                source.append(node_id[(left, cat_left)])
                target.append(node_id[(right, cat_right)])
                values.append(value)

        fig = go.Figure(
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=15, thickness=20, line=dict(color="black", width=0.5),
                    label=node_label, x=node_x, y=node_y,
                ),
                link=dict(source=source, target=target, value=values, color="rgba(203, 204, 204, 0.4)"),
            )
        )
        for i, title in enumerate(titles):
            fig.add_annotation(
                x=i / (n_cols - 1), y=1.06, xref="paper", yref="paper", text=f"<b>{title}</b>",
                showarrow=False, xanchor="left" if i == 0 else "right" if i == n_cols - 1 else "center",
            )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)
        return fig
