from unittest import TestCase

import pandas as pd

from moonstone.plot.counts import PlotCountsStats


"""
Not real tests but make sure that at least everything runs without errors
"""


class TestPlotCountsStats(TestCase):

    def test_plot_mean_distribution(self):
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
        tested_object_instance = PlotCountsStats(tested_object, items_name='species')
        tested_object_instance.plot_mean_distribution(show=False)
