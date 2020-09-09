from unittest import TestCase

import pandas as pd

from moonstone.plot.metadata import PlotMetadataStats

"""
Not real tests but make sure that at least everything runs without errors
"""


class TestPlotMetadataStats(TestCase):

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
        tested_object_instance = PlotMetadataStats(tested_object)
        tested_object_instance.plot_sex(show=False)

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
        tested_object_instance = PlotMetadataStats(tested_object)
        tested_object_instance.plot_age(show=False)
