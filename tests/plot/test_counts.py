from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.plot.counts import PlotCountsStats, PlotTaxonomyCounts


class TestPlotCountsStats(TestCase):
    def test_plot_mean_distribution(self):
        """
        Making sure that at least everything runs without errors.
        """
        tested_object = pd.DataFrame.from_dict(
            {
                "specie_1": [10, 0, 0, 2],
                "specie_2": [25, 6, 3, 9],
                "specie_3": [9, 7, 8, 3],
                "specie_4": [3, 2, 1, 0],
                "specie_5": [8, 3, 0, 1],
            },
            orient="index",
            columns=["1", "2", "3", "4"],
        )
        tested_object.columns.name = "sample"
        tested_object_instance = PlotCountsStats(tested_object, items_name="species")
        tested_object_instance.plot_mean_distribution(show=False)


class TestPlotTaxonomyCounts(TestCase):
    def setUp(self):
        self.tested_object = pd.DataFrame(
            [
                [
                    "Bacteria",
                    "Actinobacteria",
                    "Actinobacteria",
                    "Actinomycetales",
                    "Actinomycetaceae",
                    "Actinobaculum",
                    "Actinobaculum_massiliense",
                    1.0,
                    2.0,
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Lactobacillaceae",
                    "Lactobacillus",
                    "Lactobacillus (genus)",
                    0,
                    16.0,
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus (genus)",
                    1.3,
                    0.4,
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus_thermophilus",
                    1.7,
                    0.7,
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus_salivarius",
                    3.3,
                    1.2,
                ],
            ],
            columns=[
                "kingdom",
                "phylum",
                "class",
                "order",
                "family",
                "genus",
                "species",
                "SAMPLE_1",
                "SAMPLE_2",
            ],
        )
        self.tested_object = self.tested_object.set_index(
            ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        )
        self.tested_object.columns.name = "sample"
        self.tested_instance = PlotTaxonomyCounts(self.tested_object)

    def test_compute_relative_abundances_taxa_dataframe_for_n_taxa_simple_case(self):
        taxa_level = "genus"
        expected_df = pd.DataFrame(
            [
                [86.30137, 11.330049],
                [0.0, 78.817734],
            ],
            index=["Streptococcus", "Lactobacillus"],
            columns=["SAMPLE_1", "SAMPLE_2"],
        )
        expected_df.columns.name = "sample"
        expected_df.index.name = taxa_level
        pd.testing.assert_frame_equal(
            self.tested_instance._compute_relative_abundances_taxa_dataframe(
                taxa_level, 2,
            )[0],
            expected_df,
        )

    def test_compute_relative_abundances_taxa_dataframe_higher_classification(self):
        # testing higher_classification argument that allows to not take into account reads
        # defined only at an higher taxa level
        taxa_level = "species"
        expected_df = pd.DataFrame(
            [
                [45.205479, 5.911330],
                [23.287671, 3.448276],
            ],
            index=["Streptococcus_salivarius", "Streptococcus_thermophilus"],
            columns=["SAMPLE_1", "SAMPLE_2"],
        )
        expected_df.columns.name = "sample"
        expected_df.index.name = "species"
        pd.testing.assert_frame_equal(
            self.tested_instance._compute_relative_abundances_taxa_dataframe(
                taxa_level, 2, higher_classification=False
            )[0],
            expected_df,
        )

    def test_reorder_samples_clustering(self):
        origin_df = pd.DataFrame(
            [
                [1.0, 2.0, 15.0, 12.0, 16.0],
                [0.0, 16.0, 5.0, 2.0, 1.0],
                [99.0, 82.0, 80.0, 86.0, 83.0],
            ],
            index=["Actinobaculum", "Lactobacillus", "Others"],
            columns=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4", "SAMPLE_5"],
        )
        expected_df = pd.DataFrame(
            [
                [2.0, 1.0, 15.0, 12.0, 16.0],
                [16.0, 0.0, 5.0, 2.0, 1.0],
                [82.0, 99.0, 80.0, 86.0, 83.0],
            ],
            index=["Actinobaculum", "Lactobacillus", "Others"],
            columns=["SAMPLE_2", "SAMPLE_1", "SAMPLE_3", "SAMPLE_4", "SAMPLE_5"],
        )
        pd.testing.assert_frame_equal(
            self.tested_instance._cluster_samples(origin_df),
            expected_df,
        )

    def test_get_percentage_presence_genus_level_simple_case(self):
        taxa_level = "genus"
        expected_series = pd.Series(  # noqa
            [50.0, 100.0],
            index=("Lactobacillus (mean=8.00)", "Streptococcus (mean=4.30)"),
        )
        """
        pd.testing.assert_series_equal(
            self.tested_instance._get_percentage_presence(
                self.tested_object, taxa_level, 2
            ),
            expected_series,
        )
        """
        fig = self.tested_instance._plot_most_prevalent_bargraph(  # noqa
                taxa_level, 2, show=False
        )
        # print(self.tested_instance.prevalence_series)
        # print(fig)

    def test_divide_samples_into_subgroups(self):
        origin_df = pd.DataFrame(
            [
                [1.0, 2.0, 15.0, 12.0, 16.0],
                [0.0, 16.0, 5.0, 2.0, 1.0],
                [99.0, 82.0, 80.0, 86.0, 83.0],
            ],
            index=["Actinobaculum", "Lactobacillus", "Others"],
            columns=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4", "SAMPLE_5"],
        )
        sep_series = pd.Series(
            ["A", "B", "B", "A", "B"],
            index=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4", "SAMPLE_5"]
        )
        expected_df = pd.DataFrame(
            [
                [1.0, 12.0, 2.0, 15.0, 16.0],
                [0.0, 2.0, 16.0, 5.0, 1.0],
                [99.0, 86.0, 82.0, 80.0, 83.0],
            ],
            index=["Actinobaculum", "Lactobacillus", "Others"],
            columns=["SAMPLE_1", "SAMPLE_4", "SAMPLE_2", "SAMPLE_3", "SAMPLE_5"],
        )
        expected_x_coor = [(-0.5, 0.5, 1.5), (1.5, 3.0, 4.5)]
        expected_subgps = np.array(['A', 'B'], dtype=object)

        reordered_df, x_coor, subgps = self.tested_instance._divide_samples_into_subgroups_and_reorder(
            origin_df, sep_series
        )
        pd.testing.assert_frame_equal(
            reordered_df,
            expected_df,
        )
        np.testing.assert_array_equal(subgps, expected_subgps)
        self.assertListEqual(x_coor, expected_x_coor)
