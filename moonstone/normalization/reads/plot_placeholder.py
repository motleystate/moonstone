import pandas as pd
import numpy as np
import plotly.graph_objects as go


def plotting_info(input_dict):
    files = []
    reads = []
    for key in input_dict:
        files.append(key)
        reads.append(input_dict[key][2])

    df = pd.DataFrame(index=files, data=reads, columns=['reads'])
    reads = np.array(df['reads'])
    
    ds = df.describe(percentiles=[.05, .1, .25, .5, .75])
    sample_number = ds.loc['count'][0]
    mean_reads = ds.loc['mean'][0]
    std = ds.loc['std'][0]
    min_reads = ds.loc['min'][0]
    max_reads = ds.loc['max'][0]
    annotation = 'Number of Samples = %i<br>Mean Reads = %i<br>Standard Deviation = %i<br>Min Reads = %i<br>' \
                 'Max Reads = %i' % (sample_number, mean_reads, std, min_reads, max_reads)
    return reads, annotation


class PlotReads:
    def __init__(self, file_info_dict):
        self.file_info_dict = file_info_dict

    def plot_figure(self):
        reads, annotation = plotting_info(self.file_info_dict)

        trace1 = go.Box(
            name='All Samples',
            x=reads,
            quartilemethod='inclusive',
            boxmean='sd',
            boxpoints='all',
            notched=False,
            pointpos=0,
            marker_color='rgb(0, 0, 0)',
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
                text=annotation
            )]
        )

        fig = go.Figure(data=data, layout=layout)
        fig.show()





