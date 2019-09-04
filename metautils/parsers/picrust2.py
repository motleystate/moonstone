import pandas as pd


class Picrust2PathwaysParser(object):
    """
    Format is the following:

    pathways    sample_1    sample_2    sample_3
    pathway_1   14.3    123.4   12
    pathway_2   94.1    1231.1  124.2
    pathway_4   09.4    15.5    12.2
    """

    def __init__(self, file_path):
        self.file_path = file_path

    @property
    def dataframe(self):
        if not getattr(self, '_dataframe', None):
            self.to_dataframe()
        return self._dataframe

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path, sep='\t', index_col=0)
        
