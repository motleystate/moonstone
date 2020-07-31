import logging
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from moonstone.analysis.stats import (
    FilteringStats
)
from moonstone.core.module_base import ModuleBase

logger = logging.getLogger(__name__)


class Filtering(ModuleBase):
    """
    This class can be used to filter out the data where:
     - Some samples are not included in the metadata.
     - Some rows contain non relevant information

    In additon, it can be used to only get the desired rows to study.
    This class assumes that the samples are placed on the different
    columns and features (or taxa) are specified on the rows. Ex:
    samples  1  2  3  4
    feature1 x  y  t  g
    feature2 f  h  j  k
    feature3 l  m  o  p
    self.steps can be used to save the different argurments used to
    select the desired data.
    """
    def __init__(self, dataframe):
        self.counts_df = dataframe
        self.raw_df = dataframe
        self.steps = []
        self.raw_items_number = self.raw_df.shape[0]
        self.raw_reads_number = self.raw_df.sum().sum()

    def remove_data_without_metadata(self, metadata_df):
        self.counts_df = self.counts_df[list(metadata_df.columns)]
        self.steps.append('remove_data_without_metadata')
        return self.counts_df

    def remove_rows_without_relevant_info(self, name_of_rows, level_name_provided):
        self.counts_df = self.counts_df.drop(index=name_of_rows, level=level_name_provided)
        self.steps.append('remove_rows_without_relevant_info')
        return self.counts_df

    def selecting_rows(self, desired_row_series, level_to_check):
        self.counts_df = self.counts_df[self.counts_df.index.get_level_values(level_to_check).isin(desired_row_series)]
        self.steps.append('selecting_rows')
        return self.counts_df

    def deleting_only_zeros_rows(self, df):
        self.steps.append('deleting_only_zeros_rows')
        return df[df.sum(axis=1) != 0.0]


class FilteringBase(ModuleBase):

    def __init__(self, dataframe):
        self.df = dataframe
        self.steps = []
        self.raw_items_number = self.df.shape[0]
        self.raw_reads_number = self.df.sum().sum()

    @property
    def filtered_df(self):
        if getattr(self, '_filtered_df', None) is None:
            self._filtered_df = self.filter()
        return self._filtered_df

    def filter(self):
        return self.df


class MeanFiltering(FilteringBase):
    """
    This class allows to filter the data according to their mean read count,
    and remove the data that have a poor mean read count.
    You can either give a mean read count threshold or the percentage of data
    that you wish to keep (the threshold will then be computed for you).
    """
    def __init__(self, dataframe, threshold=None, percentage_to_keep=0.9):
        super().__init__(dataframe)
        self.threshold = threshold
        if not threshold:
            self.percentage_to_keep = percentage_to_keep

    def compute_threshold_best_n_percent(self):
        FS_instance = FilteringStats(self.df)
        self.items_dict, self.reads_dict = FilteringStats.by_mean(FS_instance)
        reads_dict_sort = sorted(self.reads_dict.items(), key=lambda x: x[0])  # Sorting by threshold equate to reverse
        #                                                                        sorting by value but a little quicker
        for i in reads_dict_sort:
            if i[1] / self.raw_reads_number > self.percentage_to_keep:  # i[0]:threshold; i[1]:remaining number of reads
                threshold = float('%.1f' % i[0])
                reads = int(i[1])
        logger.info('Computing filtering mean read count threshold...')
        logger.info('Filtering mean read count threshold set to %.2f.' % threshold)
        logger.info('Retaining %.1f%% of the data results in %i retained reads.'
                    % (self.percentage_to_keep * 100, reads))
        self.threshold = threshold
        return threshold

    def filter(self):
        self.steps.append('filtering_by_mean')
        if not self.threshold:
            self.compute_threshold_best_n_percent()
        logger.info('Filtering with threshold set to %.2f...' % self.threshold)
        filtered_df = self.df
        filtered_df.drop(filtered_df[filtered_df.mean(axis=1) < self.threshold].index, inplace=True)
        remaining_items = filtered_df.shape[0]
        remaining_reads = filtered_df.sum().sum()
        logger.info('Started with %i items and a total of %i reads' % (self.raw_items_number, self.raw_reads_number))
        logger.info('Saving new matrix with %i items and %i reads' % (remaining_items, remaining_reads))
        logger.info('Retained %.2f%% of items and %.2f%% of reads' % (remaining_items / self.raw_items_number * 100,
                                                                      remaining_reads / self.raw_reads_number * 100))
        return filtered_df

    def plot_threshold_vs_remaining_data(self, write_html_file, items_name='items'):
        '''  The x and y set below are are are either integers or floats.
        Be aware that some operation will require an np.array  '''
        x_items = list(self.items_dict.keys())  # get x values for plotting
        x_reads = list(self.reads_dict.keys())
        y_items = list(self.items_dict.values())  # get y values for plotting
        y_reads = list(self.reads_dict.values())
        '''  The figure visualizing the data and the filtering.  '''
        fig = make_subplots(specs=[[{"secondary_y": True}]])  # To set a second y-axis.
        fig.add_trace(
            go.Scatter(x=x_items, y=y_items, name="Retained Items: left axis"),
            secondary_y=False,
        )
        fig.add_trace(
            go.Scatter(x=x_reads, y=y_reads, name="Retained Reads: right axis"),
            secondary_y=True,
        )
        fig.add_trace(
            go.Scatter(x=[self.threshold, self.threshold], y=[0, self.raw_reads_number], name="90% Reads Threshold"),
            secondary_y=True,
        )
        # Add figure title
        fig.update_layout(
            title_text="Number of %s and reads in function of the threshold value" % items_name,
            title_x=0.5
        )
        # Set x-axis title
        fig.update_xaxes(title_text="threshold value")
        # Set y-axes titles
        fig.update_yaxes(title_text="number of %s" % items_name, secondary_y=False)
        fig.update_yaxes(title_text="number of reads", secondary_y=True)
        fig.show()

        if write_html_file:
            fig.write_html(write_html_file)

    def visualize(self, write_html_file=False):
        try:
            self.items_dict
        except NameError:
            FS_instance = FilteringStats(self.df)
            self.items_dict, self.reads_dict = FilteringStats.by_mean(FS_instance)
        self.plot_threshold_vs_remaining_data(self, write_html_file)
