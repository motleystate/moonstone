from unittest import TestCase

import numpy as np
import pandas as pd

import plotly.graph_objects as go


from moonstone.plot.graphs.scatter import GroupScatterGraph, GroupScatter3DGraph


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
        # testing that it runs without error and checking a few variables
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "group",
            show=False, plotting_options={"layout": {"height": 500, "width": 500}}
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
        # testing that it runs without error and checking a few variables
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "group",
            group_col2="smoking", show=False
        )
        self.assertEqual(tested_fig["data"][0]["marker"], go.scatter.Marker({'color': '#A63A50', 'symbol': 0}))
        self.assertEqual(tested_fig["data"][1]["marker"], go.scatter.Marker({'color': '#68ace8', 'symbol': 0}))
        self.assertEqual(tested_fig["data"][2]["marker"], go.scatter.Marker({'color': '#A63A50', 'symbol': 1}))
        self.assertEqual(tested_fig["data"][3]["marker"], go.scatter.Marker({'color': '#FFBF00', 'symbol': 1}))

    def test_translating_forced_symbols(self):
        tested_object = {"square", "103", 0}
        expected_object = {0, 103, 1}
        self.assertSetEqual(self.ins._translating_forced_symbols(tested_object), expected_object)

    def test_translating_forced_symbols_invalidtext(self):
        tested_object = {"square", "heart"}
        with self.assertRaises(ValueError):
            self.ins._translating_forced_symbols(tested_object)

    def test_symbol_scheme_under72_woforcedsymbols(self):
        tested_object = ["A", "B", "C"]
        expected_object = {"A": 0, "B": 1, "C": 2}
        np.testing.assert_array_equal(self.ins._symbol_scheme(tested_object), expected_object)

    def test_symbol_scheme_under72_wforcedsymbols(self):
        tested_object = ["A", "B", "C", "D", "E"]
        forced_symbols = {
            "A": "square",         # = 1 = "1"
            "B": "103",            # = 103 = "cross-open"
            "C": 0                 # = "0" = "circle"
        }
        expected_object = {
            "D": 2, "E": 3,
            "A": "square", "B": "103", "C": 0,
            }
        np.testing.assert_array_equal(
            self.ins._symbol_scheme(tested_object, symbols=forced_symbols),
            expected_object)

    def test_symbol_scheme_over72(self):
        tested_object = list(range(0, 73))
        expected_object = {
            0: 0, 1: 100, 2: 200, 3: 300, 4: 1, 5: 101, 6: 201, 7: 301, 8: 2, 9: 102, 10: 202, 11: 302, 12: 3, 13: 103,
            14: 203, 15: 303, 16: 4, 17: 104, 18: 204, 19: 304, 20: 5, 21: 105, 22: 205, 23: 305, 24: 6, 25: 106,
            26: 206, 27: 306, 28: 7, 29: 107, 30: 207, 31: 307, 32: 8, 33: 108, 34: 208, 35: 308, 36: 9, 37: 109,
            38: 209, 39: 309, 40: 10, 41: 110, 42: 210, 43: 310, 44: 11, 45: 111, 46: 211, 47: 311, 48: 12, 49: 112,
            50: 212, 51: 312, 52: 13, 53: 113, 54: 213, 55: 313, 56: 14, 57: 114, 58: 214, 59: 314, 60: 15, 61: 115,
            62: 215, 63: 315, 64: 16, 65: 116, 66: 216, 67: 316, 68: 17, 69: 117, 70: 217, 71: 317, 72: 18}
        np.testing.assert_array_equal(self.ins._symbol_scheme(tested_object), expected_object)

    def test_symbol_scheme_over162(self):
        tested_object = list(range(0, 163))
        expected_object = {
            0: 0, 1: 100, 2: 200, 3: 300, 4: 1, 5: 101, 6: 201, 7: 301, 8: 2, 9: 102, 10: 202, 11: 302, 12: 3, 13: 103,
            14: 203, 15: 303, 16: 4, 17: 104, 18: 204, 19: 304, 20: 5, 21: 105, 22: 205, 23: 305, 24: 6, 25: 106,
            26: 206, 27: 306, 28: 7, 29: 107, 30: 207, 31: 307, 32: 8, 33: 108, 34: 208, 35: 308, 36: 9, 37: 109,
            38: 209, 39: 309, 40: 10, 41: 110, 42: 210, 43: 310, 44: 11, 45: 111, 46: 211, 47: 311, 48: 12, 49: 112,
            50: 212, 51: 312, 52: 13, 53: 113, 54: 213, 55: 313, 56: 14, 57: 114, 58: 214, 59: 314, 60: 15, 61: 115,
            62: 215, 63: 315, 64: 16, 65: 116, 66: 216, 67: 316, 68: 17, 69: 117, 70: 217, 71: 317, 72: 18, 73: 118,
            74: 218, 75: 318, 76: 19, 77: 119, 78: 219, 79: 319, 80: 20, 81: 120, 82: 220, 83: 320, 84: 21, 85: 121,
            86: 221, 87: 321, 88: 22, 89: 122, 90: 222, 91: 322, 92: 23, 93: 123, 94: 223, 95: 323, 96: 24, 97: 124,
            98: 224, 99: 324, 100: 25, 101: 125, 102: 26, 103: 126, 104: 27, 105: 127, 106: 28, 107: 128, 108: 29,
            109: 129, 110: 30, 111: 130, 112: 31, 113: 131, 114: 32, 115: 132, 116: 33, 117: 133, 118: 34, 119: 134,
            120: 35, 121: 135, 122: 36, 123: 136, 124: 236, 125: 336, 126: 37, 127: 137, 128: 38, 129: 138, 130: 39,
            131: 139, 132: 40, 133: 140, 134: 41, 135: 141, 136: 42, 137: 142, 138: 43, 139: 143, 140: 44, 141: 144,
            142: 45, 143: 145, 144: 46, 145: 146, 146: 47, 147: 147, 148: 48, 149: 148, 150: 49, 151: 149, 152: 50,
            153: 150, 154: 51, 155: 151, 156: 52, 157: 152, 158: 53, 159: 153, 160: 54, 161: 154,
            162: 0  # -> repeating symbols because more groups than available symbols
            }
        np.testing.assert_array_equal(self.ins._symbol_scheme(tested_object), expected_object)

    def test_symbol_scheme_nothingtodo(self):
        tested_object = ["A", "B", "C"]
        tested_symbols = {"A": 0, "B": 1, "C": 101}
        expected_object = {"A": 0, "B": 1, "C": 101}
        np.testing.assert_array_equal(self.ins._symbol_scheme(tested_object, tested_symbols), expected_object)

    def test_confidence_ellipse_95(self):
        tested_x = pd.Series(
            {'sample1': 4, 'sample2': 3, 'sample3': 6, 'sample4': 8, 'sample5': 5, 'sample6': 8, 'sample7': 9})
        tested_y = pd.Series(
            {'sample1': 9, 'sample2': 5, 'sample3': 7, 'sample4': 6, 'sample5': 6, 'sample6': 9, 'sample7': 5})
        ellipse_trace_data = self.ins._confidence_ellipse(x=tested_x, y=tested_y, q_confidence=0.95, n_points=10)
        # n_points=10 to not overwhelm this file, so it's enneagone rather than an ellipse
        np.testing.assert_array_almost_equal(
            ellipse_trace_data[0],
            np.array([11.68943495, 10.53366356, 7.32338506, 3.56072442, 1.00627238, 0.85528544, 3.17841206, 6.88863549,
                      10.24990093, 11.68943495])
        )
        np.testing.assert_array_almost_equal(
            ellipse_trace_data[1],
            np.array([6.42008503, 9.16380926, 10.7613742, 10.46526146, 8.41402548, 5.56746237, 3.25751065, 2.5650224,
                      3.81402057, 6.42008503])
        )

    def test_confidence_ellipse_50(self):
        tested_x = pd.Series(
            {'sample1': 4, 'sample2': 3, 'sample3': 6, 'sample4': 8, 'sample5': 5, 'sample6': 8, 'sample7': 9})
        tested_y = pd.Series(
            {'sample1': 9, 'sample2': 5, 'sample3': 7, 'sample4': 6, 'sample5': 6, 'sample6': 9, 'sample7': 5})
        ellipse_trace_data = self.ins._confidence_ellipse(x=tested_x, y=tested_y, q_confidence=0.5, n_points=10)
        # n_points=10 to not overwhelm this file, so it's enneagone rather than an ellipse
        np.testing.assert_array_almost_equal(
            ellipse_trace_data[0],
            np.array([8.81086031, 8.25491359, 6.71071219, 4.9008051, 3.67206796, 3.59944054, 4.71690602, 6.50158987,
                      8.11841869, 8.810860315])
        )
        np.testing.assert_array_almost_equal(
            ellipse_trace_data[1],
            np.array([6.57276992, 7.89255037, 8.66100769, 8.51857216, 7.53189095, 6.16264317, 5.05151506, 4.71841582,
                      5.31920628, 6.57276992])
        )

    def test_plot_one_graph_confidence_ellipse(self):
        tested_object = pd.DataFrame(
            [
                [1, 9.3, 2],
                [12.1, 2.2, 1],
                [0, 4.5, 1],
                [9.1, 1.1, 1],
                [3.0, 5.6, 2],
                [4.3, 11.2, 2],
                [9.4, 0.1, 1],
                [1.8, 0.8, 2],
                [11.6, 2.3, 1],
                [5.6, 6.6, 1],
                [3.2, 8.9, 2],
                [3.4, 1.2, 2]
            ],
            columns=['species1', 'species2', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8', 'sample9',
                   'sample10', 'sample11', 'sample12'],
        )
        ins = GroupScatterGraph(tested_object)
        tested_fig = ins.plot_one_graph(
            first_col="species1", second_col='species2',
            group_col='group',
            q_confel=0.95, show=False)

        self.assertEqual(len(tested_fig.data), 4)   # 2 groups + 2 ellipses (one for each group)

        # first group = '2'
        self.assertEqual(tested_fig.data[0].name, '2')
        self.assertEqual(tested_fig.data[0].marker.color, '#A63A50')
        # ... and its associated ellipse
        self.assertEqual(tested_fig.data[1].hovertext, 'Ellipse 2 [95% confidence]')
        self.assertEqual(tested_fig.data[1].line, go.scatter.Line({'color': '#A63A50', 'dash': 'dot'}))

        # second/last group = '1'
        self.assertEqual(tested_fig.data[2].name, '1')
        self.assertEqual(tested_fig.data[2].marker.color, '#FFBF00')
        # ... and its associated ellipse
        self.assertEqual(tested_fig.data[3].hovertext, 'Ellipse 1 [95% confidence]')
        self.assertEqual(tested_fig.data[3].line, go.scatter.Line({'color': '#FFBF00', 'dash': 'dot'}))

    def test_plot_one_graph_confidence_ellipse_group_col2(self):
        tested_object = pd.DataFrame(
            [
                [1, 9.3, 'F', 2],
                [12.1, 2.2, 'M', 1],
                [0, 4.5, 'F', 1],
                [9.1, 1.1, 'M', 1],
                [3.0, 5.6, 'M', 2],
                [4.3, 11.2, 'F', 2],
                [9.4, 0.1, 'F', 1],
                [1.8, 0.8, 'M', 2],
                [11.6, 2.3, 'M', 1],
                [5.6, 6.6, 'F', 1],
                [3.2, 8.9, 'F', 2],
                [3.4, 1.2, 'M', 2]
            ],
            columns=['species1', 'species2', 'sex', 'group'],
            index=['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8', 'sample9',
                   'sample10', 'sample11', 'sample12'],
        )
        ins = GroupScatterGraph(tested_object)
        tested_fig = ins.plot_one_graph(
            first_col="species1", second_col='species2',
            group_col2='group', group_col='sex',
            q_confel=0.95, show=False)

        self.assertEqual(len(tested_fig.data), 8)   # 4 groups + 4 ellipses (one for each group)

        # first group = 'F - 2' (sex = 'F'; group = 2)
        self.assertEqual(tested_fig.data[0].name, 'F - 2')
        self.assertEqual(tested_fig.data[0].marker, go.scatter.Marker({'color': '#A63A50', 'symbol': 0}))
        # ... and its associated ellipse
        self.assertEqual(tested_fig.data[1].hovertext, 'Ellipse F - 2 [95% confidence]')
        self.assertEqual(tested_fig.data[1].line, go.scatter.Line({'color': '#A63A50', 'dash': 'dot'}))

        # last group = 'M - 1'
        self.assertEqual(tested_fig.data[6].name, 'M - 1')
        self.assertEqual(tested_fig.data[6].marker, go.scatter.Marker({'color': '#FFBF00', 'symbol': 1}))
        # ... and its associated ellipse
        self.assertEqual(tested_fig.data[7].hovertext, 'Ellipse M - 1 [95% confidence]')
        self.assertEqual(tested_fig.data[7].line, go.scatter.Line({'color': '#FFBF00', 'dash': 'dash'}))


class TestGroupScatter3DGraph(TestCase):

    def setUp(self):
        self.test_df = pd.DataFrame.from_dict({
            "sample1": [1.23, 4.8, 0.34, "A", "y"],
            "sample2": [-3.4, 0.87, 2.34, "B", "n"],
            "sample3": [0.99, -1.89, -3.5, "A", "n"],
            "sample4": [-2.22, -1.09, -0.019, "C", "y"],
        }, orient="index", columns=["var1", "var2", "var3", "group", "smoking"])
        self.ins = GroupScatter3DGraph(self.test_df)

    def test_plot_one_graph_1groupcol(self):
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "var3", "group",
            show=False, plotting_options={"layout": {"height": 500, "width": 500}}
        )
        self.assertEqual(tested_fig["data"][0]["marker"]["color"], "#A63A50")
        np.testing.assert_array_equal(
            tested_fig["data"][0]["x"], np.array([1.23, 0.99]))
        np.testing.assert_array_equal(
            tested_fig["data"][0]["y"], np.array([4.8, -1.89]))
        np.testing.assert_array_equal(
            tested_fig["data"][0]["z"], np.array([0.34, -3.5]))
        self.assertEqual(tested_fig["data"][0]["name"], "A")

        self.assertEqual(tested_fig["data"][1]["marker"]["color"], "#FFBF00")
        self.assertEqual(tested_fig["data"][2]["marker"]["color"], "#68ace8")

    def test_plot_one_graph_2groupcol(self):
        tested_fig = self.ins.plot_one_graph(
            "var1", "var2", "var3", "group",
            group_col2="smoking", show=False
        )
        self.assertEqual(tested_fig["data"][0]["marker"], go.scatter3d.Marker({'color': '#A63A50', 'symbol': 'circle'}))
        self.assertEqual(tested_fig["data"][1]["marker"], go.scatter3d.Marker({'color': '#68ace8', 'symbol': 'circle'}))
        self.assertEqual(tested_fig["data"][2]["marker"], go.scatter3d.Marker({'color': '#A63A50', 'symbol': 'cross'}))
        self.assertEqual(tested_fig["data"][3]["marker"], go.scatter3d.Marker({'color': '#FFBF00', 'symbol': 'cross'}))

    def test_symbol_scheme3d_under8_woforcedsymbols(self):
        tested_object = ["A", "B", "C"]
        expected_object = {"A": "circle", "B": "cross", "C": "diamond"}
        np.testing.assert_array_equal(self.ins._symbol_scheme3d(tested_object), expected_object)

    def test_symbol_scheme3d_under8_wforcedsymbols(self):
        tested_object = ["A", "B", "C", "D", "E"]
        forced_symbols = {
            "A": "square",
            "B": "diamond",
        }
        expected_object = {
            "C": "circle", "D": "cross", "E": "circle-open",
            "A": "square", "B": "diamond",
            }
        np.testing.assert_array_equal(
            self.ins._symbol_scheme3d(tested_object, symbols=forced_symbols),
            expected_object)

    def test_symbol_scheme3d_over8(self):
        tested_object = list(range(0, 9))
        expected_object = {
            0: 'circle', 1: 'cross', 2: 'diamond', 3: 'square', 4: 'circle-open', 5: 'diamond-open', 6: 'square-open',
            7: 'x', 8: 'circle'}
        np.testing.assert_array_equal(self.ins._symbol_scheme3d(tested_object), expected_object)

    def test_symbol_scheme3d_nothingtodo(self):
        tested_object = ["A", "B", "C"]
        tested_symbols = {"A": 'circle', "B": 'cross', "C": 'square-open'}
        expected_object = {"A": 'circle', "B": 'cross', "C": 'square-open'}
        np.testing.assert_array_equal(self.ins._symbol_scheme3d(tested_object, tested_symbols), expected_object)
