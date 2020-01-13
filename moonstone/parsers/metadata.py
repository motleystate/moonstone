from moonstone.analysis.columns_statistics import DataframeStatistics
from moonstone.transform.cleaning import DataFrameCleaner
from .base import BaseParser


class MetadataParser(BaseParser):

    def __init__(self, *args, cleaning_operations=None, **kwargs):
        """
        Cleaning operations are based on DataFrameCleaner object that allow to perform transformation
        operations on different columns.
        Format is the following:
        {'col_name': [
            ('operation1', 'operation1_options'),
            ('operation2', 'operation2_options')
        ]}
        """
        self.cleaning_operations = cleaning_operations
        if self.cleaning_operations is None:
            self.cleaning_operations = {}
        super().__init__(*args, **kwargs)

    def to_dataframe(self):
        dataframe = super().to_dataframe()
        df_cleaner = DataFrameCleaner(dataframe)
        for col_name, transformations in self.cleaning_operations.items():
            for transformation in transformations:
                transf_name = transformation[0]
                transf_options = transformation[1]
                df_cleaner.run_transform(col_name, transf_name, transf_options)
        return df_cleaner.df

    def get_stats(self):
        """
        :return: list of dict containing statistics about each column
        :rtype: list(dict)
        """
        return DataframeStatistics(self.dataframe).get_stats()
