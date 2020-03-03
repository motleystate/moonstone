import logging

import numpy as np

from moonstone.normalization.base import BaseNormalization

logger = logging.getLogger(__name__)


class TotalCountsNormalization(BaseNormalization):
    """
    normalization based on total counts
    """

    @property
    def scaling_factors(self):
        if getattr(self, "_scaling_factors", None) is None:
            self._scaling_factors = self.raw_df.sum() / self.raw_df.sum().mean()
        return self._scaling_factors

    def normalize(self):
        df = self.raw_df.div(self.scaling_factors)
        df = df.replace(0, np.nan).astype('str').replace(to_replace=r'^0.*$', value='1.0', regex=True)
        return df.astype('float').replace(np.nan, 0)
