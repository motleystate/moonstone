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

    def openmeta(self, outdir_path):
        df = pd.read_csv(self.metafile, sep=',', index_col=0)
        df.rename_axis("sample", inplace=True)
        df.to_csv(path_or_buf=outdir_path + '/' + 'importedMetaData.csv')
        return df
