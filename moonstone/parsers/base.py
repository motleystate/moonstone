from typing import Union

import pandas as pd


class BaseParser:
    def __init__(
        self,
        file_path: str,
        sep: str = "\t",
        no_header: bool = False,
        parsing_options: dict = None,
    ):
        """
        :param file_path: path of the input file to be parsed
        :param sep: delimiter to use (same behaviour as ``read_csv`` from pandas_)
        :param no_header: set to True if table has no header
        :param parsing_options: Extra parsing options for ``read_csv`` method, (see pandas_ documentation)

        .. _pandas: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_csv.html
        """
        self.file_path = file_path
        self.sep = sep
        self.header: Union[str, None] = "infer"
        if no_header is True:
            self.header = None
        self.parsing_options = parsing_options
        if self.parsing_options is None:
            self.parsing_options = {}

    @property
    def dataframe(self) -> pd.DataFrame:
        """Retrieve the pandas dataframe constructed from the input file."""
        if getattr(self, "_dataframe", None) is None:
            self._dataframe = self._load_data()
        return self._dataframe

    @dataframe.setter
    def dataframe(self, df: pd.DataFrame):
        """Manually set dataframe."""
        self._dataframe = df

    def _load_data(self) -> pd.DataFrame:
        """
        method that handles the loading and parsing of your file into a pandas dataframe.
        """
        return pd.read_csv(
            self.file_path, sep=self.sep, header=self.header, **self.parsing_options
        )

    @property
    def plotter(self):
        """Access to instance dedicated to visualization for this type of data."""
        if getattr(self, "_plotter", None) is None:
            self._plotter = self._instantiate_plot()
        return self._plotter

    def _instantiate_plot(self) -> pd.DataFrame:
        """
        method that handles the loading and parsing of your file into a pandas dataframe
        """
        if getattr(self, "PLOT_CLASS", None) is None:
            raise NotImplementedError("no PLOT_CLASS in %s", __class__.__name__)
        return self.PLOT_CLASS(self.dataframe)
