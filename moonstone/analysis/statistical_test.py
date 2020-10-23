import numpy as np
import pandas as pd
from scipy import stats


def mann_whitney_u_group(series, group_series, alternative='two-sided'):
    groups = list(group_series.unique())
    groups.sort()

    tab = [[np.nan] * len(groups) for _ in range(len(groups))]
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            if i != j:
                pval = stats.mannwhitneyu(series[group_series[group_series == groups[i]].index],
                                          series[group_series[group_series == groups[j]].index],
                                          alternative=alternative)[1]
                tab[i][j] = pval
                tab[j][i] = pval
    mann_whitney_u_df = pd.DataFrame(tab, index=groups, columns=groups)
    return mann_whitney_u_df
