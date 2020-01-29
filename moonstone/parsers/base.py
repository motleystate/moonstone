import pandas as pd


class BaseParser(object):

    def __init__(self, file_path, sep='\t', no_header=False, parsing_options=None):
        """
        :param parsing_options: Extra parsing options for `read_csv` method, see pandas documentation
        :type parsing_options: DICT
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
    def dataframe(self):
        if getattr(self, '_dataframe', None) is None:
            self._dataframe = self.to_dataframe()
        return self._dataframe

    def to_dataframe(self):
        return pd.read_csv(self.file_path, sep=self.sep, header=self.header, **self.parsing_options)
