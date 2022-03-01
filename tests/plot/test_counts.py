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
                    0.0,
                    4.2,
                    0.0,
                    2.0
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Lactobacillaceae",
                    "Lactobacillus",
                    "Lactobacillus (genus)",
                    0.0,
                    16.0,
                    8.0,
                    9.0
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
                    0.0,
                    3.0
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus_thermophilus",
                    0.8,
                    0.2,
                    0.7,
                    0.1
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
                    0.0,
                    2.3
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
                "SAMPLE_3",
                "SAMPLE_4"
            ],
        )
        self.tested_object = self.tested_object.set_index(
            ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        )
        self.tested_object.columns.name = "sample"
        self.tested_instance = PlotTaxonomyCounts(self.tested_object)
        # -------------
        self.mean_relab_ser_species = pd.Series({
            'Actinobaculum_massiliense': 7.821508,
            'Lactobacillus (genus)': 54.889836,
            'Streptococcus (genus)': 11.046235,
            'Streptococcus_salivarius': 20.147512,
            'Streptococcus_thermophilus': 6.094910
            })
        self.mean_relab_ser_species.index.names = ["species"]
        self.prev_ser_species = pd.Series({
            'Actinobaculum_massiliense': 50.0,
            'Lactobacillus (genus)': 75.0,
            'Streptococcus (genus)': 75.0,
            'Streptococcus_salivarius': 75.0,
            'Streptococcus_thermophilus': 100.0
            })
        self.prev_ser_species.index.names = ["species"]

    def test_relative_abundance_dataframe(self):
        expected_df = pd.DataFrame(
            [
                [
                    "Bacteria",
                    "Actinobacteria",
                    "Actinobacteria",
                    "Actinomycetales",
                    "Actinomycetaceae",
                    "Actinobaculum",
                    "Actinobaculum_massiliense",
                    0.0,
                    19.090909,
                    0.0,
                    12.195122
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Lactobacillaceae",
                    "Lactobacillus",
                    "Lactobacillus (genus)",
                    0.0,
                    72.727273,
                    91.954023,
                    54.878049
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus (genus)",
                    24.074074,
                    1.818182,
                    0.0,
                    18.292683
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus_thermophilus",
                    14.814815,
                    0.909091,
                    8.045977,
                    0.609756
                ],
                [
                    "Bacteria",
                    "Firmicutes",
                    "Bacilli",
                    "Lactobacillales",
                    "Streptococcaceae",
                    "Streptococcus",
                    "Streptococcus_salivarius",
                    61.111111,
                    5.454545,
                    0.0,
                    14.024390
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
                "SAMPLE_3",
                "SAMPLE_4"
            ],
        )
        expected_df = expected_df.set_index(
            ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        )
        expected_df.columns.name = "sample"
        pd.testing.assert_frame_equal(
            self.tested_instance.relative_abundance_dataframe,
            expected_df,
        )

    def test_prevalence_series(self):
        expected_ser = pd.Series({
            ('Bacteria', 'Actinobacteria', 'Actinobacteria', 'Actinomycetales',
             'Actinomycetaceae', 'Actinobaculum', 'Actinobaculum_massiliense'): 50.0,
            ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae',
             'Lactobacillus', 'Lactobacillus (genus)'): 75.0,
            ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae',
             'Streptococcus', 'Streptococcus (genus)'): 75.0,
            ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae',
             'Streptococcus', 'Streptococcus_thermophilus'): 100.0,
            ('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae',
             'Streptococcus', 'Streptococcus_salivarius'): 75.0
            })
        expected_ser.index.names = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
        pd.testing.assert_series_equal(
            self.tested_instance.prevalence_series,
            expected_ser,
        )

    def test_generate_list_species_to_plot_abundance(self):
        expected_most_abundant_species = pd.Index(['Lactobacillus (genus)', 'Streptococcus_salivarius'])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,
            determining_threshold=None,        # average_relative_abundance_threshold
            higher_classification=True,
            threshold_on_other_variable=None,  # prevalence_threshold
            ascending=False
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_generate_list_species_to_plot_abundance_no_higher_classification(self):
        expected_most_abundant_species = pd.Index(['Streptococcus_salivarius', 'Actinobaculum_massiliense'])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,
            determining_threshold=None,        # average_relative_abundance_threshold
            higher_classification=False,
            threshold_on_other_variable=None,  # prevalence_threshold
            ascending=False
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_generate_list_species_to_plot_abundance_determining_threshold(self):
        # species in top are all the species with an average relative abundance > 7.0%
        expected_most_abundant_species = pd.Index([
            'Lactobacillus (genus)', 'Streptococcus_salivarius',
            'Streptococcus (genus)', 'Actinobaculum_massiliense'
            ])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,                    # skipped
            determining_threshold=7.0,        # average_relative_abundance_threshold
            higher_classification=True,
            threshold_on_other_variable=None,  # prevalence_threshold
            ascending=False
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_generate_list_species_to_plot_abundance_threshold_on_other_variable(self):
        # species in top 2 need to be in at least 80% of samples
        expected_most_abundant_species = pd.Index(['Streptococcus_thermophilus'])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,
            determining_threshold=None,        # average_relative_abundance_threshold
            higher_classification=True,
            threshold_on_other_variable=80,  # prevalence_threshold
            ascending=False
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_generate_list_species_to_plot_abundance_ascending(self):
        expected_most_abundant_species = pd.Index(['Streptococcus_salivarius', 'Lactobacillus (genus)'])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,
            determining_threshold=None,        # average_relative_abundance_threshold
            higher_classification=True,
            threshold_on_other_variable=None,  # prevalence_threshold
            ascending=True
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_generate_list_species_to_plot_abundance_both_threshold(self):
        # all species with average relative abundance > 7.0% and that are present in at least 60% of the samples
        expected_most_abundant_species = pd.Index([
            'Lactobacillus (genus)', 'Streptococcus_salivarius', 'Streptococcus (genus)'
            ])
        expected_most_abundant_species.name = "species"
        tested_most_abundant_species = self.tested_instance._generate_list_species_to_plot(
            determining_ser_taxa=self.mean_relab_ser_species,
            other_variable_ser_taxa=self.prev_ser_species,
            taxa_number=2,                   # skipped
            determining_threshold=7.0,       # average_relative_abundance_threshold
            higher_classification=True,
            threshold_on_other_variable=60,  # prevalence_threshold
            ascending=False
        )
        pd.testing.assert_index_equal(
            expected_most_abundant_species,
            tested_most_abundant_species
        )

    def test_compute_relative_abundances_taxa_dataframe_species(self):
        # both threshold + no higher_classification
        taxa_level = "species"
        expected_df = pd.DataFrame(
            [
                [61.111111, 5.454545, 0.0, 14.02439],
            ],
            index=["Streptococcus_salivarius"],
            columns=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4"],
        )
        expected_df.columns.name = "sample"
        expected_df.index.name = taxa_level

        expected_ser = pd.Series(
            {'Streptococcus_salivarius': 20.14751170238975}
        )
        expected_ser.index.name = taxa_level

        tested_relab_df_taxa_top_sp, tested_taxa_number, tested_mean_relab_ser_taxa_top_sp = \
            self.tested_instance._compute_relative_abundances_taxa_dataframe(
                taxa_level=taxa_level,
                taxa_number=2,                   # skipped
                average_relative_abundance_threshold=7.0,
                higher_classification=False,
                prevalence_threshold=60,
                ascending=False
                )

        pd.testing.assert_frame_equal(tested_relab_df_taxa_top_sp, expected_df)
        pd.testing.assert_series_equal(tested_mean_relab_ser_taxa_top_sp, expected_ser)
        self.assertEqual(tested_taxa_number, 1)

    def test_compute_relative_abundances_taxa_dataframe_genus(self):
        taxa_level = "genus"
        expected_df = pd.DataFrame(
            [
                [0.0, 72.727273, 91.954023, 54.878049],
                [100.0, 8.181818, 8.045977, 32.926829]
            ],
            index=["Lactobacillus", "Streptococcus"],
            columns=["SAMPLE_1", "SAMPLE_2", "SAMPLE_3", "SAMPLE_4"],
        )
        expected_df.columns.name = "sample"
        expected_df.index.name = taxa_level

        expected_ser = pd.Series(
            {'Lactobacillus': 54.889836124066576, 'Streptococcus': 37.28865611540128}
        )
        expected_ser.index.name = taxa_level

        tested_relab_df_taxa_top_genus, tested_taxa_number, tested_mean_relab_ser_taxa_top_genus = \
            self.tested_instance._compute_relative_abundances_taxa_dataframe(
                taxa_level=taxa_level,
                taxa_number=2,                   # skipped
                average_relative_abundance_threshold=7.0,
                higher_classification=False,
                prevalence_threshold=60,
                ascending=False
                )

        pd.testing.assert_frame_equal(tested_relab_df_taxa_top_genus, expected_df)
        pd.testing.assert_series_equal(tested_mean_relab_ser_taxa_top_genus, expected_ser)
        self.assertEqual(tested_taxa_number, 2)

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

    def test_plot_most_abundant_taxa_bargraph_descending(self):
        fig = self.tested_instance._plot_most_abundant_taxa_bargraph(
            taxa_number=3,
            prevalence_threshold=None,
            ascending=False,
            show=False
        )
        expected_x = (11.04623470477129, 20.14751170238975, 54.889836124066576)
        expected_y = ("<i>Streptococcus</i> (genus)", "<i>Streptococcus salivarius</i>", 
                      "<i>Lactobacillus</i> (genus)")

        self.assertTupleEqual(fig['data'][0]['x'], expected_x)
        self.assertTupleEqual(fig['data'][0]['y'], expected_y)
        self.assertEqual(fig['layout']['title']['text'], '3 most abundant species')
        self.assertEqual(fig['layout']['yaxis']['title']['text'], 'Species')

    def test_plot_most_abundant_taxa_bargraph_ascending(self):
        fig = self.tested_instance._plot_most_abundant_taxa_bargraph(
            taxa_number=3,
            ascending=True,
            show=False
        )
        expected_x = (54.889836124066576, 20.14751170238975, 11.04623470477129)
        expected_y = ("<i>Lactobacillus</i> (genus)", "<i>Streptococcus salivarius</i>", 
                      "<i>Streptococcus</i> (genus)")

        self.assertTupleEqual(fig['data'][0]['x'], expected_x)
        self.assertTupleEqual(fig['data'][0]['y'], expected_y)

    def test_plot_most_abundant_taxa_bargraph_prevalence_threshold(self):
        fig = self.tested_instance._plot_most_abundant_taxa_bargraph(
            taxa_number=3,
            prevalence_threshold=80,
            show=False
        )
        expected_x = tuple([6.0949097082402375])
        expected_y = tuple(["<i>Streptococcus thermophilus</i>"])

        self.assertTupleEqual(fig['data'][0]['x'], expected_x)
        self.assertTupleEqual(fig['data'][0]['y'], expected_y)
        self.assertEqual(
            fig['layout']['title']['text'],
            '1 most abundant species (present in at least 80% of samples)'
            )

    def test_plot_most_abundant_taxa_boxplot(self):
        fig = self.tested_instance._plot_most_what_taxa_boxplot_or_violin(
            "abundant", "boxplot",
            taxa_number=2,
            show=False
        )

        expected_x_Ss = [61.11111111, 5.45454545, 0.0, 14.02439024]
        expected_y_Ss = [
            '<i>Streptococcus salivarius</i>', '<i>Streptococcus salivarius</i>',
            '<i>Streptococcus salivarius</i>', '<i>Streptococcus salivarius</i>'
            ]
        expected_x_L = [0.0, 72.72727273, 91.95402299, 54.87804878]
        expected_y_L = [
            '<i>Lactobacillus</i> (genus)', '<i>Lactobacillus</i> (genus)',
            '<i>Lactobacillus</i> (genus)', '<i>Lactobacillus</i> (genus)'
            ]

        np.testing.assert_allclose(fig['data'][0]['x'], expected_x_Ss)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y_Ss)
        np.testing.assert_allclose(fig['data'][1]['x'], expected_x_L)
        np.testing.assert_array_equal(fig['data'][1]['y'], expected_y_L)
        self.assertEqual(
            fig['layout']['title']['text'],
            'Relative abundance of the 2 most abundant microbial genomes among individuals \
of the cohort'
            )

    def test_plot_most_abundant_taxa_violin_and_prevalence_threshold(self):
        fig = self.tested_instance._plot_most_what_taxa_boxplot_or_violin(
            "abundant", "violin",
            taxa_number=2,
            threshold_on_other_variable=80,
            show=False
        )

        expected_x_St = [14.81481481, 0.90909091, 8.04597701, 0.6097561]
        expected_y_St = [
            '<i>Streptococcus thermophilus</i>', '<i>Streptococcus thermophilus</i>',
            '<i>Streptococcus thermophilus</i>', '<i>Streptococcus thermophilus</i>'
            ]

        np.testing.assert_allclose(fig['data'][0]['x'], expected_x_St)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y_St)
        self.assertEqual(
            fig['layout']['title']['text'],
            'Relative abundance of the 1 most abundant microbial genomes among individuals \
of the cohort (present in at least 80% of samples)'
            )

    def test_plot_most_prevalent_taxa_boxplot_mean_threshold_mean_info(self):
        fig = self.tested_instance._plot_most_what_taxa_boxplot_or_violin(
            "prevalent", "boxplot",
            taxa_number=2,
            threshold_on_other_variable=2.0,
            mean_info=True,
            show=False
        )

        expected_x_L = [0.0, 72.72727273, 91.95402299, 54.87804878]
        expected_y_L = [
            '<i>Lactobacillus</i> (genus) (mean=8.25)',
            '<i>Lactobacillus</i> (genus) (mean=8.25)',
            '<i>Lactobacillus</i> (genus) (mean=8.25)',
            '<i>Lactobacillus</i> (genus) (mean=8.25)'
            ]

        np.testing.assert_allclose(fig['data'][0]['x'], expected_x_L)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y_L)
        self.assertEqual(
            fig['layout']['title']['text'],
            'Relative abundance of the 1 most prevalent microbial genomes among individuals \
of the cohort (with mean among samples > 2.0)'
            )

    def test_plot_most_prevalent_taxa_bargraph_mean_threshold_mean_info(self):
        fig = self.tested_instance._plot_most_prevalent_taxa_bargraph(
            taxa_number=2,
            mean_threshold=2.0,
            mean_info=True,
            show=False
        )
        
        expected_x = [75.0]
        expected_y = [
            '<i>Lactobacillus</i> (genus) (mean=8.25)'
            ]

        np.testing.assert_allclose(fig['data'][0]['x'], expected_x)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y)
        self.assertEqual(
            fig['layout']['title']['text'],
            '1 most prevalent species (with mean among samples > 2.0)'
            )

    def test_plot_most_prevalent_taxa_modebargraph_plotting_options(self):
        fig = self.tested_instance.plot_most_prevalent_taxa(
            taxa_number=2,
            show=False,
            mode="bar",                 # self._valid_mode_param accepts it as "bargraph"
            plotting_options={
                'xaxes': {'type': 'log'}
            }
        )
        
        expected_x = [75.0, 100.0]
        expected_y = [
            '<i>Streptococcus salivarius</i>', '<i>Streptococcus thermophilus</i>'
            ]

        np.testing.assert_allclose(fig['data'][0]['x'], expected_x)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y)
        self.assertEqual(
            fig['layout']['title']['text'],
            '2 most prevalent species'
            )
        self.assertEqual(
            fig['layout']['xaxis']['type'],
            'log'
            )

    def test_plot_most_prevalent_taxa_modeboxplot_warning(self):
        with self.assertLogs('moonstone.plot.counts', level='WARNING') as log:
            fig = self.tested_instance.plot_most_prevalent_taxa(
                taxa_number=2,
                mean_threshold=50.0,
                show=False,
                mode="box",                 # self._valid_mode_param accepts it as "boxplot"
            )
            self.assertEqual(
                fig['layout']['title']['text'],
                'Relative abundance of the 0 most prevalent microbial genomes among individuals \
of the cohort (with mean among samples > 50.0)'
                )
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:moonstone.plot.counts:No species abide by the threshold(s) given. \
You may want to try to lower your threshold(s).", log.output)

    def test_plot_most_abundant_taxa_invalidmode(self):
        # If the mode is invalid, the graph will be plotted as a bargraph, which is the default mode
        expected_x = [20.14751170238975, 7.8215077605321515]
        expected_y = [
            '<i>Streptococcus salivarius</i>', '<i>Actinobaculum massiliense</i>'
            ]

        with self.assertLogs('moonstone.plot.counts', level='WARNING') as log:
            fig = self.tested_instance.plot_most_abundant_taxa(
                taxa_number=2,
                higher_classification=False,
                ascending=True,
                show=False,
                mode="invalidmode",
            )
            np.testing.assert_allclose(fig['data'][0]['x'], expected_x)
            np.testing.assert_array_equal(fig['data'][0]['y'], expected_y)
            self.assertEqual(
                fig['layout']['title']['text'],
                '2 most abundant species'
                )
            self.assertEqual(len(log.output), 1)
            self.assertIn(
                "WARNING:moonstone.plot.counts:mode='invalidmode' not valid, set to default (bargraph).", 
                log.output
                )

    def test_plot_most_abundant_taxa_modeviolin(self):
        fig = self.tested_instance.plot_most_abundant_taxa(
            taxa_level="genus",
            taxa_number=2,
            show=False,
            ascending=True,
            mode="violingraph",           # self._valid_mode_param accepts it as "violin"
            plotting_options={
                'xaxes':{'type': 'linear'}
            }
        )

        expected_x_L = [0.0 , 72.72727273, 91.95402299, 54.87804878]
        expected_y_L = ['<i>Lactobacillus</i>', '<i>Lactobacillus</i>', '<i>Lactobacillus</i>',
                        '<i>Lactobacillus</i>']
        expected_x_S = [100.0 , 8.18181818, 8.04597701, 32.92682927]
        expected_y_S = ['<i>Streptococcus</i>', '<i>Streptococcus</i>', '<i>Streptococcus</i>',
                        '<i>Streptococcus</i>']
        np.testing.assert_allclose(fig['data'][0]['x'], expected_x_L)
        np.testing.assert_array_equal(fig['data'][0]['y'], expected_y_L)
        np.testing.assert_allclose(fig['data'][1]['x'], expected_x_S)
        np.testing.assert_array_equal(fig['data'][1]['y'], expected_y_S)
        self.assertEqual(
            fig['layout']['title']['text'],
            'Relative abundance of the 2 most abundant microbial genomes among individuals \
of the cohort'
            )
        self.assertEqual(
            fig['layout']['xaxis']['type'],
            'linear'
            )