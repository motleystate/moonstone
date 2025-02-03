from unittest import TestCase

import numpy as np
import pandas as pd

import plotly.graph_objects as go


from moonstone.plot.graphs.scatter import GroupScatterGraph


class TestGroupScatterGraph(TestCase):

    def setUp(self):
        self.test_df = pd.DataFrame.from_dict({
            "sample1": [1.23, 4.8, "A", "y"],
            "sample2": [-3.4, 0.87, "B", "n"],
            "sample3": [0.99, -1.89, "A", "n"],
            "sample4": [-2.22, -1.09, "C", "y"],
        }, orient="index", columns=["var1", "var2", "group", "smoking"])
        self.ins = GroupScatterGraph(self.test_df)

    def test_plot_one_graph_1groupcol(self):
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "group",
            show=False
        )
        self.assertEqual(tested_fig["data"][0]["marker"]["color"], "#A63A50")
        np.testing.assert_array_equal(
            tested_fig["data"][0]["x"], np.array([1.23, 0.99]))
        np.testing.assert_array_equal(
            tested_fig["data"][0]["y"], np.array([4.8, -1.89]))
        self.assertEqual(tested_fig["data"][0]["name"], "A")

        self.assertEqual(tested_fig["data"][1]["marker"]["color"], "#FFBF00")
        self.assertEqual(tested_fig["data"][2]["marker"]["color"], "#68ace8")

    def test_plot_one_graph_2groupcol(self):
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "group",
            group_col2="smoking", show=False
        )
        self.assertEqual(tested_fig["data"][0]["marker"], go.scatter.Marker({'color': '#A63A50', 'symbol': 0}))
        self.assertEqual(tested_fig["data"][1]["marker"], go.scatter.Marker({'color': '#FFBF00', 'symbol': 0}))
        self.assertEqual(tested_fig["data"][2]["marker"], go.scatter.Marker({'color': '#68ace8', 'symbol': 0}))
        self.assertEqual(tested_fig["data"][3]["marker"], go.scatter.Marker({'color': '#A63A50', 'symbol': 1}))

    def test_symbol_scheme_under24(self):
        tested_object = ["A", "B", "C"]
        expected_object = {"A": 0, "B": 1, "C": 2}
        np.testing.assert_array_equal(self.ins._symbol_scheme(tested_object), expected_object)
