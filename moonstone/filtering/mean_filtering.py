import logging
import plotly.graph_objects as go
import plotly.io
from plotly.subplots import make_subplots

from moonstone.analysis.stats import (
    FilteringStats
)
from moonstone.filtering.base import BaseFiltering

logger = logging.getLogger(__name__)


class MeanFiltering(BaseFiltering):
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
        self._items_dict, self._reads_dict = FilteringStats.by_mean(FS_instance)
        reads_dict_sort = sorted(self._reads_dict.items(), key=lambda x: x[0])  # Sorting by threshold equate to reverse
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
        self.remaining_items_number = filtered_df.shape[0]
        self.remaining_reads_number = filtered_df.sum().sum()
        logger.info('Started with %i items and a total of %i reads' % (self.raw_items_number, self.raw_reads_number))
        logger.info('Saving new matrix with %i items and %i reads' % (self.remaining_items_number,
                                                                      self.remaining_reads_number))
        logger.info('Retained %.2f%% of items and %.2f%% of reads' % (
            self.remaining_items_number / self.raw_items_number * 100,
            self.remaining_reads_number / self.raw_reads_number * 100)
            )
        return filtered_df

    def plot_threshold_vs_remaining_data(self, write_html_file, items_name='items'):
        '''  The x and y set below are are are either integers or floats.
        Be aware that some operation will require an np.array  '''
        x_items = list(self._items_dict.keys())  # get x values for plotting
        x_reads = list(self._reads_dict.keys())
        y_items = list(self._items_dict.values())  # get y values for plotting
        y_reads = list(self._reads_dict.values())
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
            plotly.io.write_html(fig, write_html_file)

    def visualize(self, write_html_file=False):
        print(write_html_file)
        try:
            self._items_dict
        except AttributeError:
            FS_instance = FilteringStats(self.df)
            self._items_dict, self._reads_dict = FilteringStats.by_mean(FS_instance)
        self.plot_threshold_vs_remaining_data(write_html_file)

    def generate_report_data(self):
        data_dic = {
            "threshold": self.threshold,
            "n_items_removed": self.raw_items_number-self.remaining_items_number,
            "n_reads_removed": self.raw_reads_number-self.remaining_reads_number
            }
        try:
            self.percentage_to_keep
            data_dic['percentage_to_keep'] = self.percentage_to_keep
        except AttributeError:
            pass
        return {
            "title": "Filtering by mean",
            "data": data_dic
        }
