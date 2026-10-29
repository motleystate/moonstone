"""Shared input handling for the plots that show how rows flow across categorical columns."""
import logging
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

MISSING_LABEL = "(missing)"


def prepare_categories(
    data: pd.DataFrame,
    weights: Optional[pd.Series] = None,
    order: Optional[Dict[str, list]] = None,
    column_titles: Optional[Dict[str, str]] = None,
) -> Tuple[pd.DataFrame, pd.Series, Dict[str, List[str]], List[str]]:
    """
    Turn the columns of ``data`` into strings and work out the order of the categories of each column.

    Missing values are not filled in or dropped: they are shown as their own category, "(missing)", and a warning
    gives how many there are per column.

    :param data: one row per sample (or per taxon), one column per categorical variable. Plotted left to right.
    :param weights: weight of each row (default 1). Must have the same index as ``data``.
    :param order: for some columns, the categories in the order wanted. Categories left out come after, from the
        largest to the smallest.
    :param column_titles: for some columns, the title to display instead of the column name.
    :return: the string dataframe, the weights, the category order of each column, and the column titles.
    """
    if isinstance(data, pd.Series):
        data = data.to_frame()
    if data.shape[1] < 2:
        raise ValueError("At least two columns are needed to show flows between categories.")
    if data.columns.duplicated().any():
        raise ValueError("Column names must be unique.")
    if weights is None:
        weights = pd.Series(1.0, index=data.index)
    else:
        if not weights.index.equals(data.index):
            raise ValueError("weights must have the same index as data (same rows, same order).")
        if weights.isna().any() or (weights < 0).any():
            raise ValueError("weights must be non-negative numbers without missing values.")
    order = order or {}
    column_titles = column_titles or {}
    unknown = (set(order) | set(column_titles)) - set(data.columns)
    if unknown:
        raise ValueError(f"Unknown column(s) in order or column_titles: {sorted(unknown, key=str)}")

    str_df = pd.DataFrame(index=data.index)
    for col in data.columns:
        n_missing = int(data[col].isna().sum())
        if n_missing:
            logger.warning("Column %r: %d missing value(s), shown as %r.", col, n_missing, MISSING_LABEL)
        str_df[col] = data[col].astype(object).where(data[col].notna(), MISSING_LABEL).astype(str)

    category_order = {}
    for col in data.columns:
        totals = weights.groupby(str_df[col]).sum()
        # largest first, ties broken by name so that the order is reproducible
        by_size = sorted(totals.index, key=lambda cat: (-totals[cat], cat))
        wanted = [str(cat) for cat in order.get(col, [])]
        bad = [cat for cat in wanted if cat not in totals.index]
        if bad:
            logger.warning("Column %r: categories %s in order not found in the data.", col, bad)
        wanted = [cat for cat in wanted if cat in totals.index]
        category_order[col] = wanted + [cat for cat in by_size if cat not in wanted]

    titles = [str(column_titles.get(col, col)) for col in data.columns]
    return str_df, weights, category_order, titles
