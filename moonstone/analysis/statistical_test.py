import numpy as np
import pandas as pd
from scipy import stats
from typing import List, Union

from moonstone.utils.pandas.series import SeriesBinning

import logging

logger = logging.getLogger(__name__)
DEFAULT_STATS_TEST = "mann_whitney_u"


def _preprocess_groups_comparison(
    series: pd.Series, group_series: pd.Series, stat_test: str
):
    groups = list(group_series.unique())
    groups.sort()

    list_of_series = []
    new_groups = []
    for i in groups:
        groupi = series[group_series[group_series == i].index].dropna()
        if groupi.size < 1:
            logger.warning(
                f"No observations for group {i} in data. Group dropped."
            )
        elif groupi.size < 20 and stat_test == "mann_whitney_u":
            logger.warning(f"Less than 20 observations for group {i} in data.")
            list_of_series += [groupi]
            new_groups += [i]
        elif groupi.size < 10 and stat_test == "chi2_contingency":
            logger.warning(
                f"Not enough observations (<10) for group {i} in data. Group dropped."
            )
        else:
            list_of_series += [groupi]
            new_groups += [i]

    if len(new_groups) == 0:
        raise RuntimeError("All groups have been dropped: not enough observations by group.")

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
    :param output: {'series', 'dataframe'}
    :param sym: whether generated dataframe (or MultiIndexed series) is symetric or half-full

    In kwargs, you can pass argument for statistical test, like :
    :param equal_var: For ttest_ind, set to True if your samples have the same variance and
    you wish to perform a Student's t-test rather than a Welch's t-test (here default is False)
    :param alternative: {None, 'two-sided', 'less', 'greater'} For Mann - Whitney-U, you can define
    the alternative hypothesis.
    :param bins: For chi2_contingency, you can define bins.
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

    # split dataframe by group + warn and/or drop groups not respecting minimum number of observations
    groups, list_of_series = _preprocess_groups_comparison(
        series, group_series, stat_test
    )

    # if no bins defined, compute automatically bins that will be used to bin every series of group
    if stat_test == "chi2_contingency" and "bins" not in kwargs.keys():
        kwargs["bins"] = _compute_best_bins_values(list_of_series)

    if output == "series":
        dic_df = {}
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                if stat_test == "mann_whitney_u":
                    pval = mann_whitney_u(
                        list_of_series[i], list_of_series[j], **kwargs
                    )[1]
                elif stat_test == "ttest_independence":
                    pval = ttest_independence(
                        list_of_series[i], list_of_series[j], **kwargs
                    )[1]
                elif stat_test == "chi2_contingency":
                    pval = chi2_contingency(
                        list_of_series[i], list_of_series[j], **kwargs
                    )[1]

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

            if stat_test == "mann_whitney_u":
                pval = mann_whitney_u(list_of_series[i], list_of_series[j], **kwargs)[1]
            elif stat_test == "ttest_independence":
                pval = ttest_independence(
                    list_of_series[i], list_of_series[j], **kwargs
                )[1]
            if stat_test == "chi2_contingency":
                pval = chi2_contingency(
                    list_of_series[i], list_of_series[j], **kwargs
                )[1]

            tab[i][j] = pval
            if sym:
                tab[j][i] = pval
    return pd.DataFrame(tab, index=groups, columns=groups)


def mann_whitney_u(
    series1: pd.Series, series2: pd.Series, alternative: str = "two-sided", **kwargs
):
    # alternative = 'two-sided' is the default but using None gives a warning
    return stats.mannwhitneyu(series1, series2, alternative=alternative)


def ttest_independence(
    series1: pd.Series, series2: pd.Series, equal_var: bool = False, **kwargs
):
    # equal_var = False is Welch's t-test
    return stats.ttest_ind(series1, series2, equal_var=equal_var)


def _compute_best_bins_values(list_of_series):
    max_cat = 99
    max = -float("inf")
    min = float("inf")
    for series in list_of_series:
        tmp = int(series.size / 5)  # maybe add other criterium to lower the max_cat
        if tmp < max_cat:
            max_cat = tmp
        tmp = series.min()
        if tmp < min:
            min = tmp
        tmp = series.max()
        if tmp > max:
            max = tmp

    nb_cat = max_cat
    bins = None
    ind = 0
    ind_validated = []
    for series in list_of_series:
        if bins is not None:
            s_binning = SeriesBinning(series)
            s_binning.bins_values = bins
            if (
                s_binning.binned_data[s_binning.binned_data >= 5].size
                == s_binning.binned_data.size
            ):
                ind_validated += [ind]
                ind += 1
                continue
            bins = None
            nb_cat -= 1
        while nb_cat > 1:
            s_binning = SeriesBinning(series)
            bins_values = s_binning.compute_homogeneous_bins(
                min=min, max=max, nb_bins=nb_cat
            )
            s_binning.bins_values = bins_values
            if (
                s_binning.binned_data[s_binning.binned_data >= 5].size
                == s_binning.binned_data.size
            ):
                bins = (
                    s_binning.bins_values
                )  # since bins_values has changed, we need ...
                list_of_series += [
                    list_of_series[i] for i in ind_validated
                ]  # to check new bins with previous series
                ind_validated = [ind]
                break
            nb_cat -= 1
        if nb_cat < 1:
            raise RuntimeError(
                "moonstone wasn't able to compute a contingency table of at least 2 x 2 with the data."
            )
        ind += 1
    return s_binning.bins_values


def chi2_contingency(
    series1: pd.Series,
    series2: pd.Series,
    retbins: bool = False,
    bins: List[Union[int, float]] = None,
    **kwargs,
):
    """
    NB : Cells with 0 raise an error in the scipy.stats.chi2_contingency test. Furthermore, they recommand to use the
    test only if the observed and expected frequencies in each cell are at least 5.

    :param retbins: Whether to return the bins used to make the Chi2 contingency table
    """
    if series1.size < 10 or series2.size < 10:
        logger.warning(
            "Data have less than 10 observations by groups. \
        Another statistical test would be more appropriate to compare the 2 groups."
        )
        return (np.nan, np.nan)

    if bins is None:
        bins = _compute_best_bins_values([series1, series2])

    s1_binning = SeriesBinning(series1)
    s2_binning = SeriesBinning(series2)
    s1_binning.bins_values = bins
    s2_binning.bins_values = bins

    merged_df = pd.concat(
        [s1_binning.binned_data, s2_binning.binned_data],
        axis=1,
        names=["series1", "series2"],
    )
    merged_df = merged_df.fillna(0)
    to_return = list(stats.chi2_contingency(merged_df))
    to_return += s1_binning.bins_values
    return tuple(to_return)
