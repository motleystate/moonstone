import math
import numpy as np
import pandas as pd
from typing import Union

from moonstone.utils.plot import check_list_type
from moonstone.utils.convert import pandas_to_python_type


class SeriesStatsBuilder:
    def __init__(self, series):
        self.series = series

    def _build_base_stats(self):
        repartition = self.series.value_counts()
        return {
            "col_name": self.series.name,
            "col_type": str(self.series.dtype),
            "python_col_type": pandas_to_python_type(self.series.dtype),
            "n_values": self.series.size,
            "n_uniq_values": len(self.series.value_counts()),
            "values_repartition": {i: repartition[i] for i in repartition.index},
        }

    def _build_object_stats(self):
        return self._build_base_stats()

    def _build_number_stats(self):
        stats_dict = self._build_base_stats()
        stats_dict["mean"] = round(self.series.mean(), 2)
        return stats_dict

    def _build_int64_stats(self):
        return self._build_number_stats()

    def _build_float64_stats(self):
        return self._build_number_stats()

    def build_stats(self):
        return getattr(
            self, f"_build_{str(self.series.dtype)}_stats", self._build_base_stats
        )()


class SeriesBinning:
    def __init__(self, series: pd.Series):
        self.data = series

    def compute_heterogeneous_bins(self):
        """Logish bins"""
        max = self.data.max()
        magnitude = int(math.log10(max))
        bval = [-0.1, 1]  # to have the 0 value
        i = 0
        while i < magnitude:
            bval += list(np.arange(2 * 10 ** i, 10 ** (i + 1) + 1, 10 ** i))
            i += 1
        # i=magnitude
        bval += list(np.arange(2 * 10 ** i, max + 10 ** i, 10 ** i))  # up to maximum
        return bval

    def compute_homogeneous_bins(
        self,
        min: Union[int, float] = 0,
        max: Union[int, float] = None,
        nb_bins: int = None,
    ):
        """
        :param min: lower edge of bins.
        :param max: higher edge of bins (default is the dataframe's maximum).
        :param nb_bins: number of bins.
        """
        if max is None:
            max = self.data.max()
        bval = [min - 0.001]  # to have the minimum value
        if nb_bins is None:
            magnitude = int(math.log10(max))
            step = 10 ** magnitude
        else:
            step = (max - min) / nb_bins
        bval += list(np.arange(min + step, max + step, step))
        return bval

    @property
    def bins_values(self):
        """
        retrieves bins_values, and compute it if no values given
        """
        if getattr(self, "_bins_values", None) is None:
            if getattr(self, "heterogeneous", None):
                self.bins_values = self.compute_heterogeneous_bins()
            else:
                self.bins_values = self.compute_homogeneous_bins()
        return self._bins_values

    @bins_values.setter
    def bins_values(self, bins_values):
        if type(bins_values) == list and check_list_type(
            bins_values, (int, float, np.integer)
        ):
            self._bins_values = bins_values
        else:
            raise ValueError(
                "Error : expecting list of numerical values (int, float) in bins_values."
            )

    def compute_binned_data(self, normalize: bool = False, heterogeneous: bool = False):
        """
        :param heterogeneous: set to True, if you wish for heterogenous bins
        """
        self.heterogeneous = heterogeneous

        if getattr(self, "nb_bins", None) is not None:
            binned_df, bins_values = pd.cut(self.data, bins=self.nb_bins, retbins=True)
            self.bins_values = list(bins_values)
        else:
            binned_df = pd.cut(
                self.data, bins=self.bins_values
            )  # put every items in the appropriate bin
        data = pd.value_counts(binned_df, normalize=normalize)
        data = data.reindex(binned_df.cat.categories)
        new_xnames = list(data.index.astype(str))
        new_xnames[0] = new_xnames[0].replace("(-0.001", "[0.0")
        new_xnames = [new_xnames[i].replace("(", "]") for i in range(len(new_xnames))]
        data.index = new_xnames
        return data

    @property
    def binned_data(self, normalize: bool = False, heterogeneous: bool = False):
        """
        :param heterogeneous: set to True, if you wish for heterogenous bins
        """
        if getattr(self, "_binned_data", None) is None:
            self._binned_data = self.compute_binned_data(normalize, heterogeneous)
        return self._binned_data
