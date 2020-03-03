import logging
import numpy as np
from sklearn import preprocessing

logger = logging.getLogger(__name__)


class ScalingNormalization(object):
    def __init__(self, seq_counts_norm):
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.seq_counts_norm = np.array(seq_counts_norm)

    def standard(self):
        """
        ML algorithms such as SVM assume that all features are centered around zero and have
        similar variance. Scikit-learn module preprocessing.scale performs this normalization on a single array.
        More info at : https://scikit-learn.org/stable/modules/preprocessing.html
        :return:
        """

        ml_counts_norm = preprocessing.scale(self.seq_counts_norm)
        return ml_counts_norm
