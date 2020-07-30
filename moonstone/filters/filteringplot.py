import plotly.graph_objects as go
from plotly.subplots import make_subplots


class FilteringPlot:

    def __init__(self, dataframe):
        self.counts_df = dataframe
        self.raw_df = dataframe
        self.raw_items_number = self.raw_df.shape[0]
        self.raw_reads_number = self.raw_df.sum().sum()

    def threshold_vs_remaining_data(self, items_dict, reads_dict, threshold,
                                    items_name='items'):
        '''  The x and y set below are are are either integers or floats.
        Be aware that some opperation will require an np.array  '''
        x_items = list(items_dict.keys())  # get x values for plotting
        x_reads = list(reads_dict.keys())
        y_items = list(items_dict.values())  # get y values for plotting
        y_reads = list(reads_dict.values())
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
            go.Scatter(x=[threshold, threshold], y=[0, self.raw_reads_number], name="90% Reads Threshold"),
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
