from unittest import TestCase

import pandas as pd

from moonstone.plot.global_plot import (
    PlotStatsData,
    PlotStatsMetadata
)


class TestPlotStatsData(TestCase):

    def test_plot_mean(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [25, 6, 3, 9],
                'specie_3': [9, 7, 8, 3],
                'specie_4': [3, 2, 1, 0],
                'specie_5': [8, 3, 0, 1],
            },
            orient='index', columns=['1', '2', '3', '4'])
        tested_object.columns.name = 'sample'
        tested_object_instance = PlotStatsData(tested_object, items_name='species')
        tested_object_instance.plot_mean(show=False)
        # expected_object = ?????
        # self.assertEqual(tested_object_instance.plot_mean(),
        #                 expected_object)


class TestPlotStatsMetadata(TestCase):

    def test_plot_sex(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'A': [25, 'F', 'y'],
                'B': [32, 'M', 'n'],
                'C': [48, 'M', 'y'],
                'D': [68, 'F', 'n'],
                'E': [18, 'M', 'n']
            },
            orient='index', columns=['age', 'sex', 'smoker'])
        tested_object_instance = PlotStatsMetadata(tested_object)
        tested_object_instance.plot_sex(show=False)
        # expected_object = ?????

    def test_plot_age(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'A': [25, 'F', 'y'],
                'B': [32, 'M', 'n'],
                'C': [48, 'M', 'y'],
                'D': [68, 'F', 'n'],
                'E': [18, 'M', 'n']
            },
            orient='index', columns=['age', 'sex', 'smoker'])
        tested_object_instance = PlotStatsMetadata(tested_object)
        tested_object_instance.plot_age(show=False,
                                        plotting_options={'log': 'zzeezez', 'colorbar': 63, 'tickangle': 50})
        # expected_object = ?????
