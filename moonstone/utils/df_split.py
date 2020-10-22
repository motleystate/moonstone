import pandas as pd
from typing import Union

from moonstone.core.module_base import BaseDF


class DivideByGroup(BaseDF):

    def __init__(self, dataframe: Union[pd.Series, pd.DataFrame],
                 metadata_dataframe: pd.DataFrame):
        super().__init__(dataframe)
        self.metadata_df = metadata_dataframe

    def split_df(self, group_col: str, division_seq: str = None):
        """
        :param division_seq: example '1_2_3_4' will return 4 dataframes, whereas '1-2_3-4' will return 2 dataframes,
        one with group_col = 1 and 2 and the other with group_col = 3 or 4
        """
        if division_seq is None:
            division_seq = list(self.metadata_df[group_col].unique())
            division_seq = [[el] for el in division_seq]
        else:
            division_seq = division_seq.split("_")
            division_seq = [x.split('-') for x in division_seq]

        # Transform multiindex into simple index (or else merge doesn't work)
        if isinstance(self.df.index, pd.MultiIndex):
            typ = 1
            index_names = self.df.index.names
            df = self.df.reset_index()
            df_index = df[index_names]
            df = df.transpose()
        else:
            typ = 0
            df = self.df.transpose()

        merged_df = df.merge(self.metadata_df[group_col], how='inner', left_index=True, right_index=True)
        list_df = []
        for i in division_seq:
            tmp = merged_df[merged_df[group_col].isin(i)]
            if typ == 1:
                tmp = df_index.merge(tmp.transpose(), left_index=True, right_index=True)    # group_col dropped
                tmp = tmp.set_index(index_names)
            else:
                tmp = tmp.drop([group_col], axis=1).transpose()
            list_df += [tmp]
        return tuple(list_df)
