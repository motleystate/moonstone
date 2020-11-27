import pandas as pd
import numpy as np
import plotly.graph_objects as go


def plotting_info(input_dict) -> pd.DataFrame:
    """Takes info from the file_info_dict dictionary, and constructs a Pandas data frame of filenames
    and their associated number of reads."""
    files = []
    reads = []
    for key in input_dict:
        files.append(key)
        reads.append(input_dict[key][2])

    df_info = pd.DataFrame(index=files, data=reads, columns=['reads'])
    return df_info


def get_reads(df_info) -> [list, list]:
    """Takes the dataframe from plotting_info and returns the filenames and number of reads for each entry."""
    files = list(df_info.index)
    reads = np.array(df_info['reads'])
    return files, reads


def make_annotations(df_info) -> str:
    """Generates a string to be used for annotation the graph. Data is generated from the `describe` function
    in Pandas, and then parsed from the new dataframe.
    The annotation string is interpreted as HTML in the Plotly output."""
    total_reads = df_info.sum()
    ds = df_info.describe(percentiles=[.05, .1, .25, .5, .75])

    sample_number = ds.loc['count'][0]
    mean_reads = ds.loc['mean'][0]
    std = ds.loc['std'][0]
    min_reads = ds.loc['min'][0]
    max_reads = ds.loc['max'][0]
    ninety_five = ds.loc['5%']
    ninety = ds.loc['10%']
    seventy_five = ds.loc['25%']
    fifty = ds.loc['50%']
    twenty_five = ds.loc['75%']

    annotation = 'Number of Samples = %i<br>Total Reads = %i<br>' \
                 'Mean Reads = %i<br>Standard Deviation = %i<br>Min Reads = %i<br>' \
                 'Max Reads = %i<br><br>95%% of samples have at least %i reads<br>90%% of samples have at least %i' \
                 'reads<br>75%% of samples have at least %i reads<br>50%% of samples have at least %i reads<br>25%%' \
                 'of samples have at least %i reads<br>' \
                 % (sample_number, total_reads, mean_reads, std, min_reads, max_reads, ninety_five,
                    ninety, seventy_five, fifty, twenty_five)
    return annotation


class PlotReads:
    def __init__(self, file_info_dict):
        self.file_info_dict = file_info_dict

    def plot_figure(self):
        df_info = plotting_info(self.file_info_dict)
        files, reads = get_reads(df_info)
        annotation = make_annotations(df_info)

        trace1 = go.Box(
            name='All Samples',
            x=reads,
            quartilemethod='inclusive',
            boxmean='sd',
            boxpoints='all',
            notched=False,
            pointpos=0,
            marker_color="rgb(0, 0, 0)",
            fillcolor='rgba(18,23,59, 0.5)',
            marker=dict(
                size=7,
                color='rgb(0, 0, 0)'
            ),
            width=0
        )
        data = [trace1]

        layout = go.Layout(
            title=dict(
                text='Reads per Sample',
                font_size=24,
                xanchor='center',
                x=.5,
                yanchor='bottom',
                y=.85
            ),
            annotations=[dict(
                    text=annotation,
                    showarrow=False,
                    align='left'
                )],
            xaxis=dict(
                type='log',
                rangemode='normal',
                title='Number of Reads'
            )
        )

        fig = go.Figure(data=data, layout=layout)
        fig.show()
