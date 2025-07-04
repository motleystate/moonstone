from unittest import TestCase

import pandas as pd

from moonstone.utils.pandas.series import SeriesStatsBuilder, SeriesBinning


class TestSeriesStatsBuilder(TestCase):
    def test_build_stats_str(self):
        series = pd.Series(["M", "F", "M"], name="test_str")
        expected_dict = {
            "col_name": "test_str",
            "col_type": "object",
            "python_col_type": "str",
            "n_values": 3,
            "n_uniq_values": 2,
            "values_repartition": {"M": 2, "F": 1},
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)

    def test_build_stats_int(self):
        series = pd.Series([12, 41, 15], name="test_int")
        expected_dict = {
            "col_name": "test_int",
            "col_type": "int64",
            "python_col_type": "int",
            "n_values": 3,
            "n_uniq_values": 3,
            "mean": 22.67,
            "values_repartition": {
                12: 1,
                41: 1,
                15: 1,
            },
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)

    def test_build_stats_float(self):
        series = pd.Series([13.3, 15.3, 19.1], name="test_float")
        expected_dict = {
            "col_name": "test_float",
            "col_type": "float64",
            "python_col_type": "float",
            "n_values": 3,
            "n_uniq_values": 3,
            "mean": 15.9,
            "values_repartition": {
                13.3: 1,
                15.3: 1,
                19.1: 1,
            },
        }
        stats_builder = SeriesStatsBuilder(series)
        tested_dict = stats_builder.build_stats()
        self.assertDictEqual(tested_dict, expected_dict)


class TestSeriesBinning(TestCase):
    def test_compute_homogeneous_bins_under_10(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_homogeneous_bins(0, 4)
        expected_object = [0, 1, 2, 3, 4]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_homogeneous_bins_0_to_10(self):
        # testing the magnitude.is_integer() condition
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_homogeneous_bins(0, 10)
        expected_object = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_homogeneous_bins_over_10(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_homogeneous_bins(0, 14)
        expected_object = [0, 10, 20]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_homogeneous_bins_under_1(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_homogeneous_bins(0, 1)
        expected_object = [
            0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6, 0.7000000000000001, 0.8, 0.9, 1.0,
        ]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_heterogeneous_bins_0_to_54(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_heterogeneous_bins(0, 54)
        expected_object = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_heterogeneous_bins_1_to_54(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_heterogeneous_bins(1, 54)
        expected_object = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60]
        self.assertListEqual(observed_object, expected_object)

    def test_compute_heterogeneous_bins_87_to_1354(self):
        tested_object_instance = SeriesBinning(
            pd.Series(["not empty to avoid warning"])
        )
        observed_object = tested_object_instance.compute_heterogeneous_bins(87, 1354)
        expected_object = [80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000]
        self.assertListEqual(observed_object, expected_object)

    def test_bins_values_homogeneous(self):
        tested_object = pd.Series(
            {
                "gene_1": 10.5,
                "gene_2": 5.9,
                "gene_3": 9,
            }
        )
        tested_object.name = "mean read count"
        tested_object_instance = SeriesBinning(tested_object)
        tested_object_instance.bins_values  # call compute_homogeneous_bins()
        expected_object = [0, 10, 20]
        self.assertListEqual(tested_object_instance.bins_values, expected_object)

    def test_bins_values_heterogeneous(self):
        tested_object = pd.Series(
            {
                "gene_1": 10.5,
                "gene_2": 5.9,
                "gene_3": 9,
            }
        )
        tested_object.name = "mean read count"
        tested_object_instance = SeriesBinning(tested_object)
        tested_object_instance.heterogeneous = True
        tested_object_instance.bins_values  # call compute_heterogeneous_bins()
        expected_object = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20.0]
        self.assertListEqual(tested_object_instance.bins_values, expected_object)

    def test_compute_binned_data_giving_bins_values(self):
        series = pd.Series(
            [5, 8, 10], index=['i1', 'i2', 'i3']
        )
        expected_object = pd.Series(
            [1, 2], index=[']0, 5]', ']5, 10]']
        )
        expected_object.name = "count"
        tested_object_instance = SeriesBinning(series)
        tested_object_instance.bins_values = [0, 5, 10]
        tested_object_instance.farleft = False
        tested_object = tested_object_instance.compute_binned_data()
        pd.testing.assert_series_equal(tested_object, expected_object)
        pd.testing.assert_series_equal(
            tested_object_instance.binned_data, expected_object   # take the opportunity to check binned_data
        )

    def test_compute_binned_data_nbins_equalwidth(self):
        series = pd.Series(
            [0, 8, 10, 13, 14, 18, 19, 20],
            index=["i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"]
            )
        expected_object = pd.Series([3, 5], index=["[0.0, 10.0]", "]10.0, 20.0]"])
        expected_object.name = "count"
        tested_object_instance = SeriesBinning(series, nbins=2)
        tested_object = tested_object_instance.compute_binned_data()
        pd.testing.assert_series_equal(tested_object, expected_object)

    def test_compute_binned_data_nbins_equalsize(self):
        series = pd.Series(
            [0, 8, 10, 13, 14, 18, 19, 20],
            index=["i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8"]
            )
        expected_object = pd.Series([4, 4], index=["[0.0, 13.5]", "]13.5, 20.0]"])
        expected_object.name = "count"
        tested_object_instance = SeriesBinning(series, nbins=2, cut_type="equal-size")
        tested_object = tested_object_instance.compute_binned_data()
        pd.testing.assert_series_equal(tested_object, expected_object)
