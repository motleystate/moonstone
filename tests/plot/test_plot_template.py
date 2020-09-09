from unittest import TestCase

import pandas as pd

from moonstone.plot.plot_template import (
    CategoryBarGraph, DistributionBarGraph
)


class TestCategoryBarGraph(TestCase):

    def test_get_chart(self):
        series = pd.Series(
            [1, 2, 3], index=['i1', 'i2', 'i3']
        )
        expected_x = ('i1', 'i2', 'i3')
        expected_y = (1, 2, 3)
        plot = CategoryBarGraph(series)
        tested_graph = plot._get_chart()
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)

    def test_get_chart_horizontal(self):
        series = pd.Series(
            [1, 2, 3], index=['i1', 'i2', 'i3']
        )
        expected_x = (1, 2, 3)
        expected_y = ('i1', 'i2', 'i3')
        plot = CategoryBarGraph(series)
        tested_graph = plot._get_chart(orientation='h')
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)


class TestDistributionBarGraph(TestCase):

    def test_compute_heterogeneous_bins(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = DistributionBarGraph(tested_object)
        tested_object_instance.bins_values     # call compute_heterogeneous_bins()
        expected_object = [-0.1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20.0]
        self.assertListEqual(tested_object_instance.bins_values, expected_object)

    def test_get_chart(self):
        series = pd.Series(
            [5, 8, 10], index=['i1', 'i2', 'i3']
        )
        expected_x = (']0, 5]', ']5, 10]')
        expected_y = (1, 2)
        plot = DistributionBarGraph(series)
        plot.bins_values = [0, 5, 10]
        tested_graph = plot._get_chart()
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)

    def test_get_chart_horizontal(self):
        series = pd.Series(
            [5, 8, 10], index=['i1', 'i2', 'i3']
        )
        expected_x = (1, 2)
        expected_y = (']0, 5]', ']5, 10]')
        plot = DistributionBarGraph(series)
        plot.bins_values = [0, 5, 10]
        tested_graph = plot._get_chart(orientation="h")
        self.assertTupleEqual(tested_graph.x, expected_x)
        self.assertTupleEqual(tested_graph.y, expected_y)

    def test_in_bins_and_count(self):
        series = pd.Series(
            [5, 8, 10], index=['i1', 'i2', 'i3']
        )
        expected_object = pd.Series(
            [1, 2], index=[']0, 5]', ']5, 10]']
        )
        plot = DistributionBarGraph(series)
        plot.bins_values = [0, 5, 10]
        tested_object = plot.in_bins_and_count()
        pd.testing.assert_series_equal(tested_object, expected_object)
