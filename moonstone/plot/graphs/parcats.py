from typing import Dict, Optional, Union

import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

from moonstone.plot.graphs.base import BaseGraph
from moonstone.plot.graphs.category_flows import MISSING_LABEL, prepare_categories


class ParallelCategoriesGraph(BaseGraph):
    """
    Parallel categories diagram: one column per categorical variable and one ribbon per combination of categories,
    followed through every column. Unlike :class:`SankeyDiagram`, a sample can be followed across any number of
    columns, and the ribbons can be coloured by one variable.
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
        color_by: Union[str, pd.Series, None] = None,
        order: Optional[Dict[str, list]] = None,
        column_titles: Optional[Dict[str, str]] = None,
        sort_categories: str = "size",
        colors: Optional[dict] = None,
        colorscale: Union[str, list, None] = None,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
    ) -> go.Figure:
        """
        :param color_by: name of a column of data, or a series with the same index as data (so a variable that is
            not one of the plotted columns can be used), to colour the ribbons by. Default: one colour.
        :param order: for some columns, the order wanted for the categories, from top to bottom. By default, from the
            largest to the smallest.
        :param column_titles: for some columns, the title to display above the column instead of the column name
        :param sort_categories: order of the categories not listed in ``order``: ``"size"`` (largest first,
            default) or ``"alphabetical"``
        :param colors: dictionary giving a colour to some categories of color_by
        :param colorscale: for an ordered color_by (age groups, for example): a Plotly colorscale (name such as
            ``"Viridis"``, or a list of colours or of ``[position, colour]`` pairs) along which the categories are
            spread evenly, in the order given for that variable in ``order`` (alphabetical by default). Each category
            still gets exactly one colour. ``colors`` overrides it for some categories.
        :param plotting_options: options of plotting that will override the default setup
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file (.html or an image format)
        """
        order = order or {}
        color_name = color_by if isinstance(color_by, str) else getattr(color_by, "name", None)
        # color_by may be a variable that is not plotted, whose order is then only used for the colours
        plotted_order = {col: cats for col, cats in order.items() if col in self.data.columns or col != color_name}
        str_df, weights, cat_order, titles = prepare_categories(
            self.data, self.weights, plotted_order, column_titles, sort_categories)
        dimensions = [
            go.parcats.Dimension(
                values=str_df[col], label=title, categoryorder="array", categoryarray=cat_order[col],
            )
            for col, title in zip(str_df.columns, titles)
        ]

        line = {"color": self.DEFAULT_COLOR}
        if color_by is not None:
            if isinstance(color_by, str):
                color_ser = self.data[color_by]
            else:
                if not color_by.index.equals(self.data.index):
                    raise ValueError("color_by must have the same index as data (same rows, same order).")
                color_ser = color_by
            color_ser = color_ser.astype(object).where(color_ser.notna(), MISSING_LABEL).astype(str)
            categories = sorted(color_ser.unique())
            wanted = [str(cat) for cat in order.get(color_name, [])]
            categories = [cat for cat in wanted if cat in categories] + [cat for cat in categories if cat not in wanted]
            if colorscale is None:
                color_scheme = self._color_scheme_metadata(categories, colors)
            else:
                n_cat = len(categories)
                positions = [0.5] if n_cat == 1 else [i / (n_cat - 1) for i in range(n_cat)]
                sampled = sample_colorscale(colorscale, positions, colortype="rgb")
                color_scheme = {**dict(zip(categories, sampled)), **(colors or {})}
            codes = color_ser.map({cat: i for i, cat in enumerate(categories)})
            # a stepped colourscale, so that each category keeps exactly its own colour
            n = len(categories)
            colorscale = []
            for i, cat in enumerate(categories):
                colorscale += [[i / n, color_scheme[cat]], [(i + 1) / n, color_scheme[cat]]]
            # + 0.5 puts each category in the middle of its step
            line = {"color": codes + 0.5, "colorscale": colorscale, "cmin": 0, "cmax": n}

        fig = go.Figure(
            go.Parcats(
                dimensions=dimensions,
                counts=weights,
                line=line,
                hoverinfo="count",
                arrangement="freeform",
            )
        )

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)
        return fig
