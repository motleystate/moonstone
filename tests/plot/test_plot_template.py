from unittest import TestCase

import pandas as pd

from moonstone.plot.plot_template import (
    BarGraph
)


class TestBarGraph(TestCase):

    def test_compute_heterogeneous_bins(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = BarGraph(tested_object)
        tested_object_instance.bins_values     # call compute_heterogeneous_bins()
        expected_object = [-0.1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20.0]
        self.assertListEqual(tested_object_instance.bins_values, expected_object)

    def test_in_bins_and_count(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = BarGraph(tested_object)
        tested_object_instance.bins_values = [0, 5, 10, 15]
        tested_object_instance.in_bins_and_count()
        self.assertListEqual(tested_object_instance.xnames, [']0, 5]', ']5, 10]', ']10, 15]'])
        self.assertListEqual(tested_object_instance.yvalues, [0, 2, 1])

    def test_count(self):
        tested_object = pd.Series(
            {
                'A': 'F',
                'B': 'M',
                'C': 'M'
            })
        tested_object.name = 'sex'
        tested_object_instance = BarGraph(tested_object)
        tested_object_instance.count()
        self.assertListEqual(tested_object_instance.xnames, ['F', 'M'])
        self.assertListEqual(tested_object_instance.yvalues, [1, 2])

    def test_reset_xnames(self):
        tested_object = pd.Series(
            {
                'A': 'F',
                'B': 'M',
                'C': 'M'
            })
        tested_object.name = 'sex'
        tested_object_instance = BarGraph(tested_object)
        tested_object_instance.xnames = ['F', 'M']
        reset_xnames_dic = {'F': 'Female', 'M': 'Male'}
        tested_object_instance.reset_xnames(reset_xnames_dic)
        self.assertListEqual(tested_object_instance.xnames, ['Female', 'Male'])

    def test_BarGraph_numerical_nobinsvalues(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = BarGraph(tested_object, show=False)
        tested_object_instance.plot_one_graph("Mean distribution of read count of genes", "mean", "number of genes")

    def test_BarGraph_numerical_binsvalues(self):
        tested_object = pd.Series(
            {
                'gene_1': 10.5,
                'gene_2': 5.9,
                'gene_3': 9,
            })
        tested_object.name = 'mean read count'
        tested_object_instance = BarGraph(tested_object, show=False)
        tested_object_instance.bins_values = [0, 5, 10, 15]
        tested_object_instance.plot_one_graph("Mean distribution of read count of genes", "mean", "number of genes")

    def test_BarGraph_categorical(self):
        tested_object = pd.Series(
            {
                'A': 'F',
                'B': 'M',
                'C': 'M'
            })
        tested_object.name = 'sex'
        tested_object_instance = BarGraph(tested_object, show=False)
        tested_object_instance.plot_one_graph("Sex distribution in the samples", "sex", "number of samples")
