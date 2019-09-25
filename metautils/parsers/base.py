import pandas as pd


class BaseParser(object):

    def __init__(self, file_path):
        self.file_path = file_path

    @property
    def dataframe(self):
        if getattr(self, '_dataframe', None) is None:
            self.to_dataframe()
        return self._dataframe

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path)
