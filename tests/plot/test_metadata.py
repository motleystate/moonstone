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
        tested_object_instance.plot_age(show=False,
                                        plotting_options={'layout': {'yaxis_type': 'log',
                                                                     'shapes': [
                                                                         {"x0": 0, "y0": 1, "x1": 1, "y1": 1,
                                                                          "xref": "paper",
                                                                          "line": {"dash": "dash",
                                                                                   "color": "red", "width": 2}},
                                                                         {"x0": 0.5, "y0": 0, "x1": 0.5, "y1": 1,
                                                                          "xref": "paper", "yref": "paper",
                                                                          "line": {"color": "blue", "width": 2}}
                                                                          ]},
                                                          'traces': {'marker_color': 'red'},
                                                          'xaxes': {'tickangle': 50}
                                                          }
                                        )
