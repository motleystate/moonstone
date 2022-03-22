import numpy as np
from packaging import version
import pandas as pd
import scipy
from typing import List, Union, Tuple, Optional

from moonstone.utils.dict_operations import filter_dict

import logging

logger = logging.getLogger(__name__)
DEFAULT_STATS_TEST = "mann_whitney_u"
"""
TESTS_FUNCTIONS_USED = {
    'wilcoxon_rank_test': st.ranksums,
    'one_way_anova': st.f_oneway,
    'kruskal_test': st.kruskal
    }
"""


def _preprocess_groups_comparison(
    series: pd.Series, group_series: pd.Series, stat_test: str, force_computation: bool
):
    # If samples in group_series/metadata but not in series/count_dataframe
    # then we need to remove them from the group_series/metadata
    # to not get an error like "None of [Index(['sample7'], dtype='object')] are in the [index]"
    group_series_index_to_keep = group_series.index.intersection(series.index)
    if len(group_series_index_to_keep) != len(group_series.index):
        logger.info(
            "Some index values in group_series aren't found in the series. Dropping those rows."
        )
        group_series = group_series.loc[group_series_index_to_keep]
    groups = list(group_series.unique())
    groups.sort()

    if stat_test == "mann_whitney_u":
        threshold = 20
    elif stat_test == "chi2_contingency":
        threshold = 10
    else:
        threshold = 0

    list_of_series = []
    new_groups = []
    for i in groups:
        groupi_index = group_series[group_series == i].index
        groupi = series[series.index.intersection(groupi_index)].dropna()
        groupi.name = str(i)
        if groupi.size < 1:
            logger.warning(f"No observations for group {i} in data. Group dropped.")
            continue
        elif groupi.size < threshold:
            if not force_computation:
                logger.warning(
                    f"Less than {threshold} observations for group {i} in data. Group dropped"
                )
                continue  # we don't add the group to the list of series = group is dropped
            else:
                logger.warning(
                    f"Less than {threshold} observations for group {i} in data."
                )

        list_of_series += [groupi]
        new_groups += [i]

    if len(new_groups) == 0:
        raise RuntimeError(
            "All groups have been dropped: not enough observations by group."
        )

    return new_groups, list_of_series


def statistical_test_groups_comparison(
    series: pd.Series,
    group_series: pd.Series,
    stat_test: str,
    output: str = "dataframe",
    sym: str = True,
    **kwargs,
):
    """
    Args:
        output: {'series', 'dataframe'}
        sym: whether generated dataframe (or MultiIndexed series) is symetric or half-full

    In kwargs, you can pass argument for statistical test, like:
        equal_var: For ttest_ind, set to True if your samples have the same variance and
    you wish to perform a Student's t-test rather than a Welch's t-test (here default is False)
        alternative: {None, 'two-sided', 'less', 'greater'} For Mann - Whitney-U, you can define
    the alternative hypothesis.
        bins: For chi2_contingency, the criteria to bin by.
          * int : Defines the number of bins in the range of x. The range of x is extended by .1% on each side to
            include the minimum and maximum values of x.
          * sequence of scalars : Defines the bin edges. No extension of the range of x is done. So don't forget to make
            the first bin edge a little smaller than the minimum value you want included.
          * "best pvalue": Returns the contingency table associated with the best chi2 contingency pvalue
          * "max nbins": Returns the contingency table with the maximum number of bins with at least 5 observations by
            cell
    """

    method = ["mann_whitney_u", "ttest_independence", "chi2_contingency"]
    if stat_test not in method:
        logger.warning(
            "%s not a available mode, set to default (%s)",
            stat_test,
            DEFAULT_STATS_TEST,
        )
        stat_test = DEFAULT_STATS_TEST
        # raise NotImplementedError("Method %s not implemented" % stat_test)

    if kwargs.get("force_computation", None) is None:
        if stat_test == "chi2_contingency" and kwargs.get("bins", None) is None:
            kwargs["force_computation"] = False
        else:
            kwargs["force_computation"] = True

    # split dataframe by group + warn and/or drop groups not respecting minimum number of observations
    groups, list_of_series = _preprocess_groups_comparison(
        series, group_series, stat_test, force_computation=kwargs["force_computation"]
    )

    if output == "series":
        dic_df = {}
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                pval = eval(stat_test)(list_of_series[i], list_of_series[j], **kwargs)[1]
                dic_df[(groups[i], groups[j])] = pval
                if sym:
                    dic_df[(groups[j], groups[i])] = pval
        pvalue_df = pd.Series(dic_df)
        pvalue_df.index = pd.MultiIndex.from_tuples(
            pvalue_df.index, names=["Group1", "Group2"]
        )
        return pvalue_df

    tab = [[np.nan] * len(groups) for _ in range(len(groups))]
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            pval = eval(stat_test)(list_of_series[i], list_of_series[j], **kwargs)[1]
            tab[i][j] = pval
            if sym:
                tab[j][i] = pval
    return pd.DataFrame(tab, index=groups, columns=groups)


def mann_whitney_u(
    series1: pd.Series, series2: pd.Series, alternative: str = "two-sided", **kwargs
):
    # alternative = None is deprecated
    authorized_keys = ["use_continuity"]
    scipy_version = version.parse(scipy.__version__)
    if scipy_version >= version.parse("1.7.0"):
        authorized_keys += ["method", "axis"]
    if scipy_version >= version.parse("1.8.0"):
        authorized_keys += ["nan_policy"]
    kwargs, raisewarning = filter_dict(
        kwargs, authorized_keys, ["method", "axis", "nan_policy"]
    )
    if raisewarning:
        logger.warning(
            f"{raisewarning} not available with version {scipy_version} of scipy."
        )
    return scipy.stats.mannwhitneyu(series1, series2, alternative=alternative, **kwargs)


def ttest_independence(
    series1: pd.Series, series2: pd.Series, equal_var: bool = False, **kwargs
):
    # equal_var = False is Welch's t-test, which does not assume equal population variance
    authorized_keys = [
        "nan_policy",
        "axis",
    ]
    scipy_version = version.parse(scipy.__version__)
    if scipy_version >= version.parse("1.6.0"):
        authorized_keys += ["alternative"]
    if scipy_version >= version.parse("1.7.0"):
        authorized_keys += ["permutations", "random_state", "trim"]
    kwargs, raisewarning = filter_dict(
        kwargs, authorized_keys, ["alternative", "permutations", "random_state", "trim"]
    )
    if raisewarning:
        logger.warning(
            f"{raisewarning} not available with version {scipy_version} of scipy."
        )
    return scipy.stats.ttest_ind(series1, series2, equal_var=equal_var, **kwargs)


def _compute_contingency_table(
    numerical_series: pd.Series,
    categorical_series: pd.Series,
    bins: Union[List[Union[int, float]], int],
    cut_type: str = "equal-width",
    na: bool = False,
    force_computation: bool = False,
    warn: bool = True
) -> pd.DataFrame:
    """
    binning numerical series and computing contingency table
    """
    # binning numerical series
    if cut_type == "equal-size":
        binned_series = pd.qcut(numerical_series, q=bins)
    else:  # cut_type == "equal-width" (default)
        binned_series = pd.cut(numerical_series, bins=bins)
    if na:
        # crosstab doesn't count np.nan as a value
        # so we need to change it to a string for those values to appear in the contingency table
        binned_series.cat.add_categories("NaN", inplace=True)
        binned_series.fillna("NaN", inplace=True)
    # creation of the contingency table
    tab = pd.crosstab(binned_series, categorical_series)  # rows -> numerical_series; columns -> categorical_series
    if (
        tab[tab >= 5].dropna().size != tab.size
    ):  # which means that every cells have at least 5 occurences
        if warn:
            logger.warning(
                "Some cells have less than 5 observations. \
Another statistical test would be more appropriate to compare these groups."
            )
        if not force_computation:
            return np.nan
    return tab


def _compute_max_nbins_contingency_table(
    numerical_series: pd.Series,
    categorical_series: pd.Series,
    cut_type: str = "equal-width",
    na: bool = False,
    force_computation: bool = False,
) -> Tuple[pd.DataFrame, int]:
    max_nbins = int(
        categorical_series.value_counts().iloc[-1] / 5
    )  # .iloc[:,1] is the number column
    if max_nbins > 100:  # more is too much anyway
        max_nbins = 100

    nbins = max_nbins
    while nbins > 1:
        tab = _compute_contingency_table(
            numerical_series, categorical_series, nbins,
            cut_type=cut_type, na=na, force_computation=False, warn=False
        )
        if type(tab) != float:  # type(np.nan) == float
            return tab, nbins
        nbins -= 1

    # nbins == 1
    if force_computation:
        logger.warning(
            f"moonstone wasn't able to compute a contingency table with at least 5 observations per cell. \
Another statistical test would be more appropriate to compare these groups.\n\
force_computation return contingency table of 2 x {len(categorical_series.unique())}."
        )
        return _compute_contingency_table(
            numerical_series, categorical_series, 2,
            cut_type=cut_type, na=na, force_computation=True, warn=False
        ), 2
    else:
        logger.warning(
            "moonstone wasn't able to compute a contingency table with at least 5 observations per cell. \
Another statistical test would be more appropriate to compare these groups."
        )
        return np.nan, None


def _comparison_chi2_pval_depending_on_nbins(
    numerical_series: pd.Series,
    categorical_series: pd.Series,
    cut_type: str = "equal-width",
    na: bool = False,
    force_computation: bool = False,
):
    df_dic = {}

    tab, nbins = _compute_max_nbins_contingency_table(
        numerical_series, categorical_series, cut_type=cut_type, na=na, force_computation=force_computation
    )
    df_dic[nbins] = list(
        scipy.stats.chi2_contingency(tab)
        ) + [tab]
    nbins -= 1
    while nbins > 1:
        tab = _compute_contingency_table(
            numerical_series, categorical_series, nbins,
            cut_type=cut_type, na=na, force_computation=force_computation,
            warn=False
        )
        df_dic[nbins] = list(
            scipy.stats.chi2_contingency(tab)
            ) + [tab]
        nbins -= 1
    return pd.DataFrame.from_dict(
        df_dic, orient='index',
        columns=['chi2', 'pval', 'dof', 'expected table', 'observed table']
        )


def compute_contingency_table(
    numerical_series: pd.Series,
    categorical_series: pd.Series,
    bins: Union[List[Union[int, float]], int, str],
    cut_type: str = "equal-width",
    na: bool = False,
    force_computation: bool = False,
    retpvalcomparisondf: bool = False
):
    """
    binning numerical series and computing contingency table

    Try to find the maximum number of bins, where every bins have at least 5 occurrences

    Args:
        numerical_series: Series that needs to be put into bins
        categorical_series: Other series used to make the contingency table
        bins: The criteria to bin by.
          * int : Defines the number of bins in the range of x. The range of x is extended by .1% on each side to
            include the minimum and maximum values of x.
          * sequence of scalars : Defines the bin edges. No extension of the range of x is done. So don't forget to make
            the first bin edge a little smaller than the minimum value you want included.
          * "best pvalue": Returns the contingency table associated with the best chi2 contingency pvalue
          * "max nbins": Returns the contingency table with the maximum number of bins with at least 5 observations by
            cell
        cut_type: {"equal-width" (default), "equal-size"} Defines how the bins are cut when a integer is given to bins.
          Note: they are defined on numerical_series only.
          So all cells of the contingency table might not all have the same weight at the end.
        na: Defines if NaN should be considered as a value
        force_computation: Returns contingency table even if not every cell has at least 5 observations
        retpvalcomparisondf: Returns the comparison dataframe used to evaluate the number of bins that gives the lowest
          chi2 contingency pvalues. (used with bins set to "best pvalue")
    """
    if na:
        categorical_series = categorical_series.replace(np.nan, "NaN")
    samples = numerical_series.index.intersection(categorical_series.index)
    if len(samples) == 0:
        raise ValueError("Index of numerical_series and categorical_series don't match")

    categorical_series = categorical_series.loc[samples]
    numerical_series = numerical_series.loc[samples]

    if numerical_series.name is None:
        numerical_series.name = "numerical_series"
    if categorical_series.name is None:
        categorical_series.name = "categorical_series"

    if bins == "max nbins":
        return _compute_max_nbins_contingency_table(
            numerical_series,
            categorical_series,
            cut_type=cut_type,
            na=na,
            force_computation=force_computation,
        )[0]
    elif bins == "best pvalue":
        comparison_df = _comparison_chi2_pval_depending_on_nbins(
            numerical_series,
            categorical_series,
            cut_type=cut_type,
            na=na,
            force_computation=force_computation,
        )
        nbins_best_pvalue = comparison_df["pval"].iloc[::-1].idxmin()
        if retpvalcomparisondf:
            return comparison_df.loc[nbins_best_pvalue]["observed table"], comparison_df
        return comparison_df.loc[nbins_best_pvalue]["observed table"]
    else:
        return _compute_contingency_table(
            numerical_series,
            categorical_series,
            bins,
            na=na,
            cut_type=cut_type,
            force_computation=force_computation
        )


def _add_category_column(series, defaultname: str):
    df = pd.DataFrame(series)
    if series.name is not None:
        name = series.name
    else:
        name = defaultname
    df.columns = ["number"]
    df["category"] = name
    return df


def chi2_contingency(
    series1: pd.Series,
    series2: pd.Series,
    bins: Union[List[Union[int, float]], int, str] = "best pvalue",
    cut_type: str = "equal-width",
    na: bool = False,
    force_computation: bool = False,
    rettab: bool = False,
    **kwargs,
) -> Tuple[float, float, int, np.ndarray, Optional[pd.DataFrame]]:
    """
    NB : Cells with 0 raise an error in the scipy.stats.chi2_contingency test. Furthermore, they recommand to use the
    test only if the observed and expected frequencies in each cell are at least 5.

    :param rettab: Whether to return the Chi2 contingency table.
    """
    if (series1.size < 10 or series2.size < 10) and not force_computation:
        logger.warning(
            "Data have less than 10 observations by groups. \
Another statistical test would be more appropriate to compare those 2 groups."
        )
        to_return = [np.nan] * 4
        if rettab:
            to_return += [np.nan]
        return tuple(to_return)

    df1 = _add_category_column(series1, defaultname="series1")
    df2 = _add_category_column(series2, defaultname="series2")
    df = df1.append(df2)

    if bins == "best pvalue":
        tab, comparison_df = compute_contingency_table(
            df["number"],
            df["category"],
            bins,
            cut_type=cut_type,
            na=na,
            force_computation=force_computation,
            retpvalcomparisondf=True
        )
        nbins_best_pvalue = tab.shape[0]
        if rettab:
            return comparison_df.loc[nbins_best_pvalue].to_list()
        else:
            return comparison_df.loc[nbins_best_pvalue].to_list()[:-1]
    else:
        tab = compute_contingency_table(
            df["number"],
            df["category"],
            bins,
            cut_type=cut_type,
            na=na,
            force_computation=force_computation,
        )
        to_return = list(
            scipy.stats.chi2_contingency(tab)
        )  # ch2, pval, degree of freedom, expected table
        if rettab:
            to_return += [tab]  # observed table
        return tuple(to_return)
