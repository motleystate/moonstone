import numpy as np
import pandas as pd
from scipy import stats
import logging
module_logger = logging.getLogger(__name__)


def count_items(y):
    from collections import Counter

    c = Counter(y)
    c = dict(c)

    sample_count = 0
    for _, value in c.items():
        sample_count += value

    print(f"\nThere appear to be a total of {sample_count} samples.")
    for category, value in c.items():
        print(" %s samples labeled at %s, or %2.1f%s of the total"
              % (value, category, value / sample_count * 100, "%"))


def normalized_stats(x):
    properties = ['Number of Variable', 'Min / Max', 'Mean', 'Variance', 'Skewness', 'Kurtosis / Fisher']
    df_statistics = pd.DataFrame(x).describe().transpose()
    means = np.array(df_statistics['mean'])

    for i in range(len(stats.describe(means))):
        print(f"\t{properties[i]} = {stats.describe(means)[i]}")


class Descriptive(object):
    def __init__(self, df, outdir):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.df = df
        self.outdir = outdir

    def matrix_stats(self, filename):
        output_file = self.outdir+'/'+filename
        df_stats = self.df.describe().transpose()  # Descriptive stats of variables
        df_stats.to_csv(path_or_buf=output_file)  # The above saved to a file (before filtering)
        means = np.array(df_stats['mean'])

        print(f"Saving variables to {output_file} and running descriptive statistics on means:")

        properties = ['Number of Variable', 'Min / Max', 'Mean', 'Variance', 'Skewness', 'Kurtosis / Fisher']
        for i in range(len(stats.describe(means))):
            print(f"\t{properties[i]} = {stats.describe(means)[i]}")

        # print(stats.describe(means))  # Additional stats available through SCIPY

        # Generating a plot could be an option.
        # Probably as part of generalized PDF export.
        # df_means = pd.DataFrame(df_stats['mean'])
        # df_means.boxplot()
        # plt.show()


class FilteringStats(object):
    def __init__(self, df):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.df = df

    def by_mean(self):
        """Compute number of remaining items and reads for different thresholds"""
        items_dict = {}  # key = filtering value, value = # of remaining items
        reads_dict = {}
        counts_df_mean = pd.concat([self.df.sum(axis=1), self.df.mean(axis=1)], axis=1)
        for x in np.arange(.1, 10, .1):  # This fits the current project but need to be generalized.
            counts_df_mean.drop(counts_df_mean[counts_df_mean[1] < x].index, inplace=True)
            items_dict[x] = counts_df_mean.shape[0]
            reads_dict[x] = counts_df_mean[0].sum()
        return items_dict, reads_dict


"""Count data from metagenomic experiments is usually space.
Definitions can vary, but generally >50% zeros who indicate a sparse matrix.
Normalization / Standardization before running algorithms usually results is improved performance and, indeed,
 some techniques require scaled data.
 Methods of scaling count data should take account of sparcity of data.

 This Class will calculate sparsity for a given count matrix"""


class Density(object):
    def __init__(self, count_matrix):
        self.count_matrix = count_matrix

    def is_sparse(self):
        x = self.count_matrix

        # This could be written as a single line. But the calculation is not long and it is clearer
        # if we proceed stepwise.

        # Total number of entries in the matrix
        matrix_entries = x.shape[0]*x.shape[1]

        # The number of NON-zero entries. Numpy is handy here.
        matrix_not_zeros = np.count_nonzero(x)

        percent_zeros = ((matrix_entries-matrix_not_zeros)/matrix_entries)*100

        if percent_zeros >= 50:
            return True
        else:
            return False

    def percent_non_zeros(self):
        # Essentially the same as 'is_sparse' method, but returning the actually percentage of non-zeros
        x = self.count_matrix
        matrix_entries = x.shape[0] * x.shape[1]

        matrix_not_zeros = np.count_nonzero(x)

        percent_zeros = (matrix_entries - matrix_not_zeros) / matrix_entries
        return (1-percent_zeros)*100
