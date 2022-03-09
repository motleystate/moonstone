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
    def __init__(
        self,
        series: pd.Series,
        nbins: int = None,
        cut_type: str = "equal-width",
        bins_values: list = None,
        farleft: bool = True,
        heterogenous: bool = False,
    ):
        """
        Args:
            series: Series to bin.
            bins_values: List of bins' edges.
            nbins: Number of bins. (skipped by bins_values)
            cut_type: {"equal-width" (default), "equal-size"} how the bins are defined. (used with nbins and skipped by
              bins_values)
            farleft: To use with bins_values. If set to True (default), the far left edges is included in the first bin.

        If nbins and bins_values left undefined, bins_values will be automatically defined starting at 0 and
        ending with the bin containing the maximum value in the data.
        """
        self.data = series
        self.nbins = nbins
        self.cut_type = cut_type
        if bins_values is not None:
            self.bins_values = bins_values  # setter doesn't accept None
        self.farleft = farleft
        self.heterogeneous = heterogenous

    def compute_heterogeneous_bins(
        self, min: Union[int, float] = 0, max: Union[int, float] = None
    ):
        """
        Logish bins

        :param min: lower edge of bins.
        :param max: higher edge of bins (default is the dataframe's maximum).
        """
        if max is None:
            max = self.data.max()
        magnitude = int(math.log10(max))
        if min == 0:
            bval = [0]
            min = 1
        else:
            bval = []
        i = int(math.log10(min))
        bval += list(
            np.arange(int(str(min)[0]) * 10**i, 10 ** (i + 1) + 1, 10**i)
        )  # starting at min
        i += 1
        while i < magnitude:
            bval += list(np.arange(2 * 10**i, 10 ** (i + 1) + 1, 10**i))
            i += 1
        # i=magnitude
        bval += list(np.arange(2 * 10**i, max + 10**i, 10**i))  # up to maximum
        return bval

    def compute_homogeneous_bins(
        self,
        min: Union[int, float] = 0,
        max: Union[int, float] = None,
        nbins: int = None,  # not using self.nbins to have the ability to call on its own
    ):
        """
        :param min: lower edge of bins.
        :param max: higher edge of bins (default is the dataframe's maximum).
        :param nbins: number of bins.
        """
        if max is None:
            max = self.data.max()
        bval = [min]
        # bval = [min - 0.001]  # to have the minimum value
        # self._farleftedges = "["+str(min)
        if nbins is None:
            magnitude = math.log10(max)
            if magnitude.is_integer():
                magnitude = int(magnitude) - 1  # so that 0 to 10, gives [0, 1, ..., 10]
                # instead of [0, 10]
            else:
                magnitude = int(magnitude)
            step = 10**magnitude
        else:
            step = (max - min) / nbins
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
        if self.farleft:
            self._farleftedges = "[" + str(self._bins_values[0])
            return [self._bins_values[0] - 0.001] + self._bins_values[1:]
        self._farleftedges = "]" + str(self._bins_values[0])
        return self._bins_values

    @bins_values.setter
    def bins_values(self, bins_values):
        if isinstance(bins_values, list) and check_list_type(
            bins_values, (int, float, np.integer)
        ):
            self._bins_values = bins_values
        else:
            raise ValueError(
                "Error : expecting list of numerical values (int, float) in bins_values."
            )

    def compute_binned_data(self):
        """
        If no bins_values given during instantiation, will call
        :param heterogeneous: set to True, if you wish for heterogenous bins
        """
        if getattr(self, "_bins_values", None) is None and self.nbins is not None:
            if self.cut_type == "equal-size":
                binned_df, bins_values = pd.qcut(self.data, q=self.nbins, retbins=True)
            else:  # default = "equal-width"
                binned_df, bins_values = pd.cut(
                    self.data, bins=self.nbins, retbins=True
                )
            self.bins_values = list(bins_values)
        else:
            binned_df = pd.cut(
                self.data, bins=self.bins_values
            )  # put every items in the appropriate bin
        data = pd.value_counts(binned_df)
        data = data.reindex(binned_df.cat.categories)
        new_xnames = list(data.index.astype(str))
        tmp = new_xnames[0].split(",")
        if getattr(self, "_farleftedges", None) is not None:
            tmp[0] = self._farleftedges
        else:
            tmp[0] = "[" + str(
                float(self.data.min())
            )  # float because all numbers in the other intervals written as float
        new_xnames[0] = ",".join(tmp)
        new_xnames = [new_xnames[i].replace("(", "]") for i in range(len(new_xnames))]
        data.index = new_xnames
        return data

    @property
    def binned_data(self):
        """
        :param heterogeneous: set to True, if you wish for heterogenous bins
        """
        if getattr(self, "_binned_data", None) is None:
            self._binned_data = self.compute_binned_data()
        return self._binned_data
