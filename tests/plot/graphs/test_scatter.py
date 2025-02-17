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

    def test_dictionary_marker_symbols(self):
        expected_object = {
            "circle": 0, "circle-open": 100, "circle-dot": 200, "circle-open-dot": 300, "square": 1,
            "square-open": 101, "square-dot": 201, "square-open-dot": 301, "diamond": 2, "diamond-open": 102,
            "diamond-dot": 202, "diamond-open-dot": 302, "cross": 3, "cross-open": 103, "cross-dot": 203,
            "cross-open-dot": 303, "x": 4, "x-open": 104, "x-dot": 204, "x-open-dot": 304, "triangle-up": 5,
            "triangle-up-open": 105, "triangle-up-dot": 205, "triangle-up-open-dot": 305, "triangle-down": 6,
            "triangle-down-open": 106, "triangle-down-dot": 206, "triangle-down-open-dot": 306, "triangle-left": 7,
            "triangle-left-open": 107, "triangle-left-dot": 207, "triangle-left-open-dot": 307, "triangle-right": 8,
            "triangle-right-open": 108, "triangle-right-dot": 208, "triangle-right-open-dot": 308, "triangle-ne": 9,
            "triangle-ne-open": 109, "triangle-ne-dot": 209, "triangle-ne-open-dot": 309, "triangle-se": 10,
            "triangle-se-open": 110, "triangle-se-dot": 210, "triangle-se-open-dot": 310, "triangle-sw": 11,
            "triangle-sw-open": 111, "triangle-sw-dot": 211, "triangle-sw-open-dot": 311, "triangle-nw": 12,
            "triangle-nw-open": 112, "triangle-nw-dot": 212, "triangle-nw-open-dot": 312, "pentagon": 13,
            "pentagon-open": 113, "pentagon-dot": 213, "pentagon-open-dot": 313, "hexagon": 14, "hexagon-open": 114,
            "hexagon-dot": 214, "hexagon-open-dot": 314, "hexagon2": 15, "hexagon2-open": 115, "hexagon2-dot": 215,
            "hexagon2-open-dot": 315, "octagon": 16, "octagon-open": 116, "octagon-dot": 216, "octagon-open-dot": 316,
            "star": 17, "star-open": 117, "star-dot": 217, "star-open-dot": 317, "hexagram": 18, "hexagram-open": 118,
            "hexagram-dot": 218, "hexagram-open-dot": 318, "star-triangle-up": 19, "star-triangle-up-open": 119,
            "star-triangle-up-dot": 219, "star-triangle-up-open-dot": 319, "star-triangle-down": 20,
            "star-triangle-down-open": 120, "star-triangle-down-dot": 220, "star-triangle-down-open-dot": 320,
            "star-square": 21, "star-square-open": 121, "star-square-dot": 221, "star-square-open-dot": 321,
            "star-diamond": 22, "star-diamond-open": 122, "star-diamond-dot": 222, "star-diamond-open-dot": 322,
            "diamond-tall": 23, "diamond-tall-open": 123, "diamond-tall-dot": 223, "diamond-tall-open-dot": 323,
            "diamond-wide": 24, "diamond-wide-open": 124, "diamond-wide-dot": 224, "diamond-wide-open-dot": 324,
            "hourglass": 25, "hourglass-open": 125, "bowtie": 26, "bowtie-open": 126, "circle-cross": 27,
            "circle-cross-open": 127, "circle-x": 28, "circle-x-open": 128, "square-cross": 29,
            "square-cross-open": 129, "square-x": 30, "square-x-open": 130, "diamond-cross": 31,
            "diamond-cross-open": 131, "diamond-x": 32, "diamond-x-open": 132, "cross-thin": 33,
            "cross-thin-open": 133, "x-thin": 34, "x-thin-open": 134, "asterisk": 35, "asterisk-open": 135, "hash": 36,
            "hash-open": 136, "hash-dot": 236, "hash-open-dot": 336, "y-up": 37, "y-up-open": 137, "y-down": 38,
            "y-down-open": 138, "y-left": 39, "y-left-open": 139, "y-right": 40, "y-right-open": 140, "line-ew": 41,
            "line-ew-open": 141, "line-ns": 42, "line-ns-open": 142, "line-ne": 43, "line-ne-open": 143, "line-nw": 44,
            "line-nw-open": 144, "arrow-up": 45, "arrow-up-open": 145, "arrow-down": 46, "arrow-down-open": 146,
            "arrow-left": 47, "arrow-left-open": 147, "arrow-right": 48, "arrow-right-open": 148, "arrow-bar-up": 49,
            "arrow-bar-up-open": 149, "arrow-bar-down": 50, "arrow-bar-down-open": 150, "arrow-bar-left": 51,
            "arrow-bar-left-open": 151, "arrow-bar-right": 52, "arrow-bar-right-open": 152, "arrow": 53,
            "arrow-open": 153, "arrow-wide": 54, "arrow-wide-open": 154}
        self.assertDictEqual(self.ins.dictionary_marker_symbols, expected_object)

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
