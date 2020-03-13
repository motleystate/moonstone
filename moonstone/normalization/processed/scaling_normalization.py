import logging
import numpy as np
from sklearn import preprocessing

from moonstone.normalization.processed.base import BaseScaler
from moonstone.analysis import stats

logger = logging.getLogger(__name__)


class StandardScaler(BaseScaler):
    """
    ML algorithms such as SVM assume that all features are centered around zero and have
    similar variance. Scikit-learn module preprocessing.scale performs this normalization on a single array.
    More info at : https://scikit-learn.org/stable/modules/preprocessing.html
    :return:
    """

    def __init__(self, raw_x):
        super().__init__(raw_x)
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        if not isinstance(self.raw_x, np.ndarray):
            raise ValueError('A NumPy array is required for normalization. Got {}'.format(type(self.raw_x)))

    def scale(self):
        """
        Takes a NumPy array of the independent variables, or features, as 'x' for ML training.
        """

        logger.info("Counts standardized by Standard Scalar. Check Mean ~0.0 and similar Variances:")
        scaled_x = preprocessing.scale(self.raw_x, with_mean=True, with_std=True)
        stats.normalized_stats(scaled_x)
        return scaled_x
