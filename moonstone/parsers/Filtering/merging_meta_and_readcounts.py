import pandas as pd


class MergingMetaAndReadCounts:
    """
    Both, on the metadata_df and the reads_df assumes that the different samples
    are displayed horizontally and the features(or taxa) vertically:

    samples  1  2  3
    feature1 x  t  c
    feature2 y  g  a
    feature3 v  b  n
    """

    def __init__(self, metadata_df, reads_df):
        self.metadata_df = metadata_df
        self.read_count_df = reads_df
        self.metadata_df.columns = self.metadata.columns.astype(str)
        self.read_count_df.columns = self.reads.columns.astype(str)

    def check_column_names(self):
        return self.metadata.columns.equals(self.reads.columns)

    @property
    def full_dataframe(self):
        try: 
            if getattr(self, "_full_dataframe", None) is None:
                if self.check_column_names() is True:
                    meta_and_reads_df = pd.concat([self.metadata_df, self.read_count_df], axis=0, sort=False)
                    setattr(self, "_full_dataframe", meta_and_reads_df)
            return self._full_dataframe
        except:
                print('Cannot concat metadata dataframe and reads dataframe because\
                        column names or number of samples do not match.')

    @property
    def full_dataframe_with_features_in_columns(self):
        if getattr(self, "_full_dataframe_with_features_in_columns", None) is None:
            full_dataframe_with_features_in_columns = self.full_dataframe.transpose()
            setattr(self, "_full_dataframe_with_features_in_columns", full_dataframe_with_features_in_columns)
        return self._full_dataframe_with_features_in_columns
