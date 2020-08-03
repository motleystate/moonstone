import pandas as pd


class BaseParser:

    def __init__(self, file_path: str, sep: str = '\t', no_header: bool = False, parsing_options: dict = None):
        """
        :param file_path: path of the input file to be parsed
        :param sep: delimiter to use (same behaviour as ``read_csv`` from pandas_)
        :param no_header: set to True if table has no header
        :param parsing_options: Extra parsing options for ``read_csv`` method, (see pandas_ documentation)

        .. _pandas: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html
        """
        self.file_path = file_path
        self.sep = sep
        self.header = 'infer'
        if no_header is True:
            self.header = None
        self.parsing_options = parsing_options
        if self.parsing_options is None:
            self.parsing_options = {}

    @property
    def dataframe(self) -> pd.DataFrame:
        """
        retrieve the pandas dataframe constructed from the input file
        """
        if getattr(self, '_dataframe', None) is None:
            self._dataframe = self.to_dataframe()
        return self._dataframe

    def to_dataframe(self) -> pd.DataFrame:
        """
        method that handles the loading and parsing of your file into a pandas dataframe
        """
        return pd.read_csv(self.file_path, sep=self.sep, header=self.header, **self.parsing_options)
