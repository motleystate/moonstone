class FilteringOptions:
    """
    This class can be used to filter out the data where:
     - Some samples are not included in the metadata.
     - Some rows contain non relevant information

    In additon, it can be used to only get the desired rows to study.
    """
    def __init__(self, dataframe):
        self.counts_df = dataframe
        self.raw_dataframe = dataframe
        self.steps = []

    def remove_data_without_metadata(self, metadata_df):
        self.counts_df = self.counts_df[list(metadata_df.columns)]
        self.steps.append('remove_data_without_metadata')
        return self.counts_df

    def remove_rows_wihtout_relevant_info(self, name_of_rows, level_name_provided):
        self.counts_df = self.counts_df.drop(index=name_of_rows, level=level_name_provided)
        self.steps.append('remove_rows_wihtout_relevant_info')
        return self.counts_df

    def selecting_certain_rows(self, desired_row_series, level_to_check):
        self.counts_df = self.counts_df[self.counts_df.index.get_level_values(level_to_check).isin(desired_row_series)]
        self.steps.append('selecting_certain_rows')
        return self.counts_df
