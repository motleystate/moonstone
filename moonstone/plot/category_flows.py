"""Plots showing how samples (or taxa) flow across categorical levels."""
import logging
from typing import Dict, List, Optional, Union

import pandas as pd
import plotly.graph_objects as go

from moonstone.plot.graphs.parcats import ParallelCategoriesGraph
from moonstone.plot.graphs.sankey import SankeyDiagram

logger = logging.getLogger(__name__)

KINDS = {"sankey": SankeyDiagram, "parcats": ParallelCategoriesGraph}
OTHER_LABEL = "Other"


def _graph_class(kind: str):
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {sorted(KINDS)}, not {kind!r}.")
    return KINDS[kind]


def plot_category_flows(
    data: pd.DataFrame,
    columns: Optional[List[str]] = None,
    kind: str = "sankey",
    order: Optional[Dict[str, list]] = None,
    column_titles: Optional[Dict[str, str]] = None,
    color_by: Union[str, pd.Series, None] = None,
    weights: Optional[pd.Series] = None,
    plotting_options: dict = None,
    show: bool = True,
    output_file: Union[bool, str] = False,
) -> go.Figure:
    """
    Show how samples flow across several categorical columns, for instance of the metadata.

    :param data: one row per sample, one column per variable
    :param columns: the columns to plot, from left to right (default: all of them)
    :param kind: ``"sankey"`` (boxes with counts, ribbons between neighbouring columns; best for two or three columns)
        or ``"parcats"`` (parallel categories: one ribbon per sample through every column, which can be coloured
        with ``color_by``; best for following samples across many columns)
    :param order: for some columns, the order wanted for the categories, from top to bottom. By default, from the
        largest to the smallest.
    :param column_titles: for some columns, the title to display instead of the column name
    :param color_by: (``"parcats"`` only) a column of data, or a series with the same index, to colour the ribbons by
    :param weights: weight of each row (default 1, so the widths are numbers of samples)
    :param plotting_options: options of plotting that will override the default setup
    :param show: set to False if you don't want to show the plot
    :param output_file: name of the output file (.html or an image format)

    Missing values are not filled in or dropped: they are shown as their own category, "(missing)".
    """
    graph_class = _graph_class(kind)
    plotted = data if columns is None else data[columns]
    graph = graph_class(plotted, weights=weights)
    kwargs = dict(
        order=order, column_titles=column_titles, plotting_options=plotting_options,
        show=show, output_file=output_file,
    )
    if kind == "parcats":
        if isinstance(color_by, str) and color_by not in plotted.columns:
            color_by = data[color_by]
        kwargs["color_by"] = color_by
    elif color_by is not None:
        raise ValueError('color_by is only available with kind="parcats".')
    return graph.plot_one_graph(**kwargs)


def _lump_into_other(flows_df: pd.DataFrame, weights: pd.Series, top_n: int) -> pd.DataFrame:
    """Keep the top_n categories of each column (by total weight) and call all the others "Other"."""
    flows_df = flows_df.copy()
    for col in flows_df.columns:
        totals = weights.groupby(flows_df[col]).sum()  # missing values are left out of the ranking
        kept = set(sorted(totals.index, key=lambda cat: (-totals[cat], str(cat)))[:top_n])
        lumped = flows_df[col].notna() & ~flows_df[col].isin(kept)
        flows_df[col] = flows_df[col].where(~lumped, OTHER_LABEL)
    return flows_df


def plot_taxonomy_flows(
    counts: pd.DataFrame,
    levels: List[str],
    kind: str = "sankey",
    top_n: int = 10,
    weight: str = "mean_relative_abundance",
    order: Optional[Dict[str, list]] = None,
    column_titles: Optional[Dict[str, str]] = None,
    plotting_options: dict = None,
    show: bool = True,
    output_file: Union[bool, str] = False,
) -> go.Figure:
    """
    Show how the taxa of a count table collapse into higher taxa, for instance species -> genus -> family.

    :param counts: count table, one row per taxon (index with one level per taxonomic rank), one column per sample
    :param levels: the taxonomic levels to plot, from left to right, as named in the index of counts
    :param kind: ``"sankey"`` (default) or ``"parcats"``. See :func:`plot_category_flows`.
    :param top_n: number of taxa to keep in each level; the others are lumped into "Other". This only affects what
        is drawn: the weights are computed beforehand from all the taxa. Use ``None`` to keep every taxon.
    :param weight: what the width of a ribbon stands for. ``"mean_relative_abundance"`` (default): the relative
        abundance (%), averaged across samples. ``"total_reads"``: the sum of the reads across samples.
    :param order: for some levels, the order wanted for the taxa, from top to bottom. By default, from the largest to
        the smallest.
    :param column_titles: for some levels, the title to display instead of the level name
    :param plotting_options: options of plotting that will override the default setup
    :param show: set to False if you don't want to show the plot
    :param output_file: name of the output file (.html or an image format)

    A taxon not classified down to one of the levels (missing in the index) is shown as "(missing)" in that level.
    """
    graph_class = _graph_class(kind)
    if weight == "mean_relative_abundance":
        row_weights = (counts / counts.sum() * 100).mean(axis=1)
    elif weight == "total_reads":
        row_weights = counts.sum(axis=1)
    else:
        raise ValueError('weight must be "mean_relative_abundance" or "total_reads".')
    missing_levels = [lvl for lvl in levels if lvl not in counts.index.names]
    if missing_levels:
        raise ValueError(f"Level(s) {missing_levels} not found in the index of counts ({list(counts.index.names)}).")

    # one row per combination of levels (a genus can hold many species, and many taxa can share the same levels)
    row_weights = row_weights.reset_index(drop=True)
    flows_df = counts.index.to_frame(index=False)[levels]
    combos = flows_df.assign(_weight=row_weights).groupby(levels, dropna=False, sort=False)["_weight"].sum()
    flows_df = combos.index.to_frame(index=False)
    weights = pd.Series(combos.to_numpy(), index=flows_df.index)
    if top_n is not None:
        flows_df = _lump_into_other(flows_df, weights, top_n)
    return graph_class(flows_df, weights=weights).plot_one_graph(
        order=order, column_titles=column_titles, plotting_options=plotting_options,
        show=show, output_file=output_file,
    )
