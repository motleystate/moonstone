from unittest import TestCase

import pandas as pd

from moonstone.plot.graphs.histogram import Histogram


class TestHistogram(TestCase):

    def setUp(self):
        self.series = pd.Series({
            'sample1': 47,
            'sample2': 33,
            'sample3': 17,
            'sample4': 23,
            'sample5': 41,
            'sample6': 25
        })

    def test_with_bins(self):
        # Not real test but make sure that visualize() runs without errors
        plot = Histogram(self.series)
        plot.plot_one_graph(bins_size=10, show=False)
