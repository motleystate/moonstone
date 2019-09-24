import pandas as pd

from .base import BaseParser


class Picrust2PathwaysParser(BaseParser):
    """
    Format is the following:

    pathways    sample_1    sample_2    sample_3
    pathway_1   14.3    123.4   12
    pathway_2   94.1    1231.1  124.2
    pathway_4   09.4    15.5    12.2
    """

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path, sep='\t', index_col=0)
