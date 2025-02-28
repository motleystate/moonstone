import logging
import itertools
import pandas as pd
import numpy as np
import scipy
from typing import Union

logger = logging.getLogger(__name__)


class StructureRemodelling(object):

    _SYM_BOOL_DIC = {True: "symmetric", False: "asymmetric"}

    def __init__(self, data_structure: Union[pd.Series, pd.DataFrame], sym: bool = None):
        """
        Args:
            data_structure: pandas dataframe or series to remodel
            sym: is the data_structure symmetrical. {None (default), True, False}.
                 If set to None, then the symmetry will be checked beforehand
        """
        sym = self._valid_sym_param_init(sym=sym)

        if isinstance(data_structure, pd.Series):
            if sym is None:
                sym = self.check_symmetric_series(data_structure)
            if sym is True:
                print("object stored as a symmetric series")
                self._symmetric_series = data_structure
            else:
                print("object stored as an asymmetric series")
                self._asymmetric_series = data_structure

        elif isinstance(data_structure, pd.DataFrame):
            if sym is None:
                sym = self.check_symmetric_dataframe(data_structure)
            if sym is True:
                print("object stored as a symmetric dataframe")
                self._symmetric_dataframe = data_structure
            else:
                print("object stored as an asymmetric dataframe")
                self._asymmetric_dataframe = data_structure
    
    def _valid_sym_param_init(self, sym: bool):
        if isinstance(sym, bool) or sym is None:
            return sym
        raise ValueError("sym='{}' not valid. Valid choices: [True, False, None]. \
If you're not sure of the symmetry of your data structure, use None.".format(sym))
        return None

    def _valid_sym_param_get_from_arguments(self, sym: bool):
        if isinstance(sym, bool):
            return sym
        raise ValueError("sym='{}' not valid. Valid choices: [True, False].".format(sym))
    
    def _valid_structure_param(self, structure: str):
        choices = ["series", "dataframe"]
        if structure in choices:
            return structure
        raise ValueError("structure='{}' not valid. Valid choices: {}".format(structure, choices))

    def get_from_arguments(self, structure: str, sym: bool):
        sym = self._valid_sym_param_get_from_arguments(sym=sym)
        structure = self._valid_structure_param(structure=structure)
        return getattr(self, self._SYM_BOOL_DIC[sym]+"_"+structure) 

    def check_symmetric_series(self, ser):
        return all(ser.get((a, b), None) == ser.get((b, a), None) for a, b in ser.index)

    def check_symmetric_dataframe(self, df):
        return scipy.linalg.issymmetric(df.to_numpy())

    def _make_series_symmetric(self, ser):
        ser = pd.concat([ser, ser.reorder_levels([1, 0])])
        ser = ser[~ser.index.duplicated(keep='first')]  # if onself there (ex: ("A", "A"))
        labels = set(ser.index.get_level_values(0))
        all_pairs = list(itertools.product(labels, repeat=2))
        return ser.reindex(all_pairs, fill_value=np.nan).sort_index()

    def _make_dataframe_symmetric(self, df):
        # If the first column and last line aren't there. We need to add it before making it symmetric.
        # Example:
        #   B   C   D           A   B   C   D           A   B   C   D
        # A 1   2   3       A   NaN 1   2   3       A   NaN 1   2   3
        # B NaN 4   5   =>  B   NaN NaN 4   5   =>  B   1   NaN 4   5
        # C NaN NaN 6       C   NaN NaN NaN 6       C   2   4   NaN 6
        #                   D   NaN NaN NaN NaN     D   3   5   6   NaN
        df = df.copy()
        dc = list(set(df.columns).difference(df.index))
        if len(dc) == 1:
            df.loc[dc[0]] = np.nan
        dr = list(set(df.index).difference(df.columns))
        if len(dr) == 1:
            df[dr[0]] = np.nan
            df = df[df.index]   # reordering not essential but so it's looks more like
            # a symetric dataframe instead of a messy one
        df = df.fillna(df.transpose())
        return df

    def _make_dataframe_asymmetric(self, df):
        # https://stackoverflow.com/questions/74767255/how-to-transform-a-multi-dimensional-symmetric-dataframe-to-a-2-column-dataframe
        return df.where(np.triu(np.ones(df.shape, dtype=bool))).dropna(how="all", axis=0).dropna(how="all", axis=1)

    def _make_series_asymmetric(self, ser):
        # the goal of an asymmetric data structure is to reduce the size so we remove all the redundant data
        # as well as all the NaN.
        return ser[~ser.index.map(frozenset).duplicated()].dropna()

    def _transform_series_to_dataframe(self, ser):
        df = ser.unstack(level=1)
        df.index.name = None
        df.columns.name = None
        return df

    def _transform_dataframe_to_series(self, df, dropna: bool = False, onself: bool = True):
        """
        Args:
            dropna: dropna argument of panda's stack function.
            onself: whether to keep values on self (True) or not (False).
                    Example of value on self: ("A", "A")
        """
        ser = df.stack(dropna=dropna)
        if onself is True:
            return ser
        else:  # onself is False
            return ser[ser.index.get_level_values(0) != ser.index.get_level_values(1)]

    @property
    def symmetric_series(self):
        if getattr(self, '_symmetric_series', None) is None:
            if getattr(self, '_asymmetric_series', None) is not None:
                self._symmetric_series = self._make_series_symmetric(self.asymmetric_series)
            else:
                self._symmetric_series = self._transform_dataframe_to_series(self.symmetric_dataframe, dropna=False)
        return self._symmetric_series

    @property
    def symmetric_dataframe(self):
        if getattr(self, '_symmetric_dataframe', None) is None:
            if getattr(self, '_asymmetric_dataframe', None) is not None:
                self._symmetric_dataframe = self._make_dataframe_symmetric(self.asymmetric_dataframe)
            else:
                self._symmetric_dataframe = self._transform_series_to_dataframe(self.symmetric_series)
        return self._symmetric_dataframe

    @property
    def asymmetric_series(self):
        if getattr(self, '_asymmetric_series', None) is None:
            if getattr(self, '_symmetric_series', None) is not None:
                self._asymmetric_series = self._make_series_asymmetric(self.symmetric_series)
            else:
                self._asymmetric_series = self._transform_dataframe_to_series(self.asymmetric_dataframe, dropna=True)
        return self._asymmetric_series

    @property
    def asymmetric_dataframe(self):
        if getattr(self, '_asymmetric_dataframe', None) is None:
            if getattr(self, '_symmetric_dataframe', None) is not None:
                self._asymmetric_dataframe = self._make_dataframe_asymmetric(self.symmetric_dataframe)
            else:
                self._asymmetric_dataframe = self._transform_series_to_dataframe(self.asymmetric_series)
        return self._asymmetric_dataframe
