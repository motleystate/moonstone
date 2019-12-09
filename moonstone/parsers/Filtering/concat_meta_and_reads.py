import pandas as pd


class ConcatMetaAndReads:

    def __init__(self, metadata_df, reads_df):
        self.metadata = metadata_df
        self.reads = reads_df
        self.metadata.columns = self.metadata.columns.astype(str)
        self.reads.columns = self.reads.columns.astype(str)

    def check_column_names(self):
        return self.metadata.columns.equals(self.reads.columns)

    @property
    def full_dataframe(self):
        if getattr(self, "_full_dataframe", None) is None:
            if self.check_column_names() is True:
                meta_and_reads_df = pd.concat([self.metadata, self.reads], axis=0, sort=False)
                setattr(self, "_full_dataframe", meta_and_reads_df)
                return self._full_dataframe
            else:
                return 'Cannot concat metadata dataframe and reads dataframe because\
                        column names or number of samples do not match.'

    @property
    def full_dataframe_with_features_in_columns(self):
        if getattr(self, "_full_dataframe_with_features_in_columns", None) is None:
            full_dataframe_with_features_in_columns = self.full_dataframe.transpose()
            setattr(self, "_full_dataframe_with_features_in_columns", full_dataframe_with_features_in_columns)
        return self._full_dataframe_with_features_in_columns
