from unittest import TestCase

import pandas as pd

from moonstone.plot.graphs.bargraph import BarGraph


class TestBarGraph(TestCase):

    def test_get_chart(self):
        series = pd.Series(
            [1, 2, 3], index=['i1', 'i2', 'i3']
        )
        expected_x = ('i1', 'i2', 'i3')
        expected_y = (1, 2, 3)
        plot = BarGraph(series)
        tested_graph = plot._get_chart()
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)

    def test_get_chart_horizontal(self):
        series = pd.Series(
            [1, 2, 3], index=['i1', 'i2', 'i3']
        )
        expected_x = (1, 2, 3)
        expected_y = ('i1', 'i2', 'i3')
        plot = BarGraph(series)
        tested_graph = plot._get_chart(orientation='h')
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)
