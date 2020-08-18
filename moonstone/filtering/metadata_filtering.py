from typing import List

import pandas as pd

from moonstone.core.module_base import BaseModule
from moonstone.filtering.base import BothAxisFiltering


class Filtering(BaseModule):
    """
    @TODO to be removed
    This is a deprecated class that remains as long as all part of moonstone using it is updated:
    * main.py
    * differentialanalysis.py

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
