import logging
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class MergeDF(object):
    """MergeDF uses the Pandas function 'merge' to combine count data and metadata.
    The purpose of merging in this way is to faithfully combine clinical and sequence data in the same dataframe,
    using the sample names.
    """
    def __init__(self, countfile, metadata, variable):
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.dc = countfile
        self.dm = metadata
        self.variable = variable
        self._merged_df = None

    def merge(self):
        """Performs merging on DadaFrame indexes.
        Note: Pandas data frames requires indexes to be of the same type in order to correctly merge.
        We will check for mismatched index type and try to correct
        """

        logger.info('Merge function called to merge count data and metadata.')
        logger.info(f'Variable {self.variable} from metadata file will be merged with counts.')

        if not isinstance(self.dc.index, type(self.dm.index)):
            logger.warning(f'Index types do not match: {type(self.dc.index)} and {type(self.dm.index)}.')
            self.dc.set_index(np.int64(np.array(self.dc.index)), inplace=True)
            logger.info(f' Indexes reset. Count Index={type(self.dc.index)}, Metadata Index={type(self.dm.index)}')

        df = pd.merge(self.dm[self.variable], self.dc, left_index=True, right_index=True)
        logger.info('Merge function completed. Returning merged data frame.')
        return df

    @property
    def merged_df(self):
        if getattr(self, "_merged_df", None) is None:
            self._merged_df = self.merge()
        return self._merged_df
