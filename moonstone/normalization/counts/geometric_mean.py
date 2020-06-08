import logging
import math

import numpy as np
import pandas as pd

from moonstone.normalization.base import BaseNormalization

logger = logging.getLogger(__name__)


class GeometricMeanNormalization(BaseNormalization):
    """
    normalization based on the one performed by DeSeq2 : https://bioconductor.org/packages/release/bioc/html/DESeq2.html

    info: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
    """

    def __init__(self, df, log_number=np.e, zero_threshold=80, normalization_level=None, replace_0_to_1=False):
        """
        :param normalization_level: At which level of a multi-index you want the normalization to be perfomed
        """
        super().__init__(df)
        if replace_0_to_1 is True:
            self.df = self.df.replace(0, 1)
        self.log_number = log_number
        self.zero_threshold = zero_threshold
        self.normalization_level = normalization_level
        if normalization_level is not None and isinstance(self.df.index, pd.MultiIndex):
            self.grouped_df = self.df.groupby(level=self.normalization_level).sum()
            logger.info("Normalization on %s level (n=%s)", self.normalization_level, self.grouped_df.shape[0])
        else:
            self.grouped_df = self.df
            logger.info("Normalization on all rows (n=%s)", self.grouped_df.shape[0])

    def non_zero_df(self, df):
        """
        This method removes rows with 0 reads
        """
        threshold = math.ceil(df.shape[1] * self.zero_threshold/100)
        total_nb_rows = df.shape[0]
        non_zero_dataf = df.replace(0, np.nan).dropna(thresh=threshold).astype('float')
        removed_nb_rows = non_zero_dataf.shape[0]
        logger.info("%s/%s rows dropped", total_nb_rows - removed_nb_rows, total_nb_rows)
        if removed_nb_rows / total_nb_rows <= 0.5:
            logger.warning("Zero-filtering has removed %s %% of items!",
                           (total_nb_rows - removed_nb_rows) / total_nb_rows*100)
        self._removed_zero_df = df[~df.index.isin(non_zero_dataf.index)]
        return non_zero_dataf

    def log_df(self, df):
        return df.applymap(lambda x: math.log(x, self.log_number))

    @property
    def removed_zero_df(self):
        """
        gives the dataframe with the rows that were removed for having too many zeros.
        this attribute is computed during the non_zero_df function
        """
        if getattr(self, "_removed_zero_df", None) is None:
            logger.warning("Computing the scaling factors beforehand is required to access this dataframe")
            return None
        else:
            return self._removed_zero_df

    def remove_zero_and_apply_log(self, df):
        return self.log_df(self.non_zero_df(df))

    def calculating_and_substracting_mean_row(self, df):
        """
        Substracting the mean row to original values
        """
        return df.sub(df.mean(axis=1), axis='rows')

    @property
    def scaling_factors(self):
        if getattr(self, "_scaling_factors", None) is None:
            non_zero_log_df = self.remove_zero_and_apply_log(self.grouped_df)
            substracted_mean_df = self.calculating_and_substracting_mean_row(non_zero_log_df)
            while substracted_mean_df.rpow(self.log_number).median().isna().any():
                logging.warning('Zero filtering of %i is too strict to compute scaling factors, trying %i' %
                                (self.zero_threshold, self.zero_threshold - 5))
                self.zero_threshold = self.zero_threshold - 5
                non_zero_log_df = self.remove_zero_and_apply_log(self.grouped_df)
                substracted_mean_df = self.calculating_and_substracting_mean_row(non_zero_log_df)

            scaling_factors = substracted_mean_df.rpow(self.log_number).median()
            setattr(self, "_scaling_factors", scaling_factors)
        return self._scaling_factors

    def normalize(self):
        return self.raw_df.div(self.scaling_factors)
