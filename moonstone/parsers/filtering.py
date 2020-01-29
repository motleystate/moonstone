"""
@TODO
This needs to be moved and adapted to the filtering module
"""

import logging
module_logger = logging.getLogger(__name__)


class Filtering(object):
    def __init__(self, df, mean_filter):
        # self.logger = logging.getLogger('moonstone.filtering.Filtering')
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.df = df
        self.mean_filter = mean_filter

    def by_mean(self):
        self.logger.info('Started: Filtering by mean number of reads.')
        df_filtered = self.df[self.df.columns[self.df.mean() >= self.mean_filter]]
        num_vars = df_filtered.shape[1]
        self.logger.info('Retained %i taxa with at least %d mean reads.' % (num_vars, self.mean_filter))
        self.logger.info('Completed: Filtering by mean number of reads. Have a nice day!')
        return df_filtered
