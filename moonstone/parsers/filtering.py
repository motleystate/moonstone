import logging
logging.basicConfig(filename='output/info.log', filemode='w', level=logging.DEBUG)


class Filtering(object):
    def __init__(self, df, mean_filter):
        self.df = df
        self.mean_filter = mean_filter

    def by_mean(self):
        logging.info('Started: Filtering by mean number of reads.')
        df_filtered = self.df[self.df.columns[self.df.mean() >= self.mean_filter]]
        num_vars = df_filtered.shape[1]
        logging.info('Retained %i taxa with at least %d mean reads.' % (num_vars, self.mean_filter))
        logging.info('Completed: Filtering by mean number of reads. Have a nice day!')
        return df_filtered
