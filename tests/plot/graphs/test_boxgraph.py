import numpy as np
from unittest import TestCase

import pandas as pd

from moonstone.plot.graphs.box import (
    GroupBoxGraph
)


class TestGroupBoxGraph(TestCase):

    def test_orientation_vertical(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M"],
                [3.0, "F"],
                [9.0, "M"],
                [6.0, "M"],
                [2.0, "F"]
            ],
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
            columns=["data", "group"],
        )
        expected_x_gpM = ['M', 'M', 'M']
        expected_y_gpM = [1.0, 9.0, 6.0]
        expected_x_gpF = ['F', 'F']
        expected_y_gpF = [3.0, 2.0]
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="group", show=False
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpM)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpM)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpF)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpF)

    def test_orientation_horizontal(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M"],
                [3.0, "F"],
                [9.0, "M"],
                [6.0, "M"],
                [2.0, "F"]
            ],
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
            columns=["data", "group"],
        )
        expected_x_gpM = [1.0, 9.0, 6.0]
        expected_y_gpM = ['M', 'M', 'M']
        expected_x_gpF = [3.0, 2.0]
        expected_y_gpF = ['F', 'F']
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="group", show=False,
            orientation="h"
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpM)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpM)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpF)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpF)
