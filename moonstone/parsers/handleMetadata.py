"""
Here is where the clinical data for the samples will be opened.
This file is likely to be highly variable depending on the project
    and will likely need several versions of the import code.
At the least, samples names/references should match with the 'counts' file
"""

import pandas as pd


class Inputs(object):
    def __init__(self, metafile):
        self.metafile = metafile

    def openmeta(self):
        df = pd.read_csv(self.metafile, sep=',', index_col='sample')
        # df.drop(['SUBJID'], axis=1, inplace=True)
        return df
