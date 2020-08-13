import pandas as pd
from typing import Optional

from moonstone.plot.plot_template import (
    HistogramTemplate
)

"""
plot that could be used on any dataframe no matter the module
Or should we put that in the analysis module?
"""


# What people might want to visualize?

class PlotStatsData():

    def __init__(self, dataframe: pd.DataFrame):
        self.df = dataframe

    def plot_mean(self, output_file: Optional[str] = False):
        df_mean = self.df.mean(axis=1)
        binsvalues = HistogramTemplate.compute_asymetric_bins(
            df_mean.max(),
            negative=False
            )
        binned_df = HistogramTemplate.in_bins(
            self.df_mean,
            binsvalues
            )
        HistogramTemplate.one_hist(binned_df, log=True, output_file=output_file)


class PlotStatsMetadata():

    def __init__(self, metadata_dataframe: pd.DataFrame):
        self.metadata_df = metadata_dataframe

    # plot stats about metadata : proportions M/F etc. ?
