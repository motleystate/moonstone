"""Here is where the normalized count data will be opened.
 The file needs to be opened formatted and some basic statistics calculated"""

import pandas as pd


class Inputs(object):
    def __init__(self, countfile):
        self.countfile = countfile

    def opencounts(self):
        # Normalized Counts from SHAMAN
        # This should have been processed with Deseq2
        df = pd.read_csv(self.countfile, sep=',')
        # df = df.rename(columns={"Unnamed: 0": "sample"})  # Samples are numbered but the column name is empty.
        df = df.transpose()  # OTUs are now columns, and samples are rows.

        # OTU report can generate unidentified taxa, which are labeled 'NA' by SHAMAN
        # These need to be identified and removed.
        # The following statements looks for 'NaN' in row 1.
        # If 'NaN' is found, the whole column is removed.
        if df.iloc[0].isnull().any():  # NEED TO FILTER OUT NAN BEFORE PROCEEDING!
            print(f'Removing {df.iloc[0].isnull().any().sum()} unnamed/unidentified variable!')
            # print(df.columns[df.iloc[0].isnull()])
            df.drop(df.columns[df.iloc[0].isnull()], axis=1, inplace=True)

        df.columns = df.iloc[0]
        # After transposition, the far-left index numbers become the column names.
        # This corrected by renaming columns by the 1st item which is the OTU name.

        df = df.drop(df.index[0], axis=0)
        # formerly written as 'df = df.drop('Samples')'
        # Removes the first row (copied to be the column names).
        # Resulting matrix has OTUs (taxa) as columns and samples as rows.

        df.rename_axis("", axis="columns", inplace=True)  # Should remove extraneous column axis name
        df.rename_axis("sample", inplace=True)  # Gives the row axis the proper 'sample' name
        df = df.applymap(float)  # Converts all values to floats in preparation of analysis
        return df

