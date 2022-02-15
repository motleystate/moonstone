import numpy as np
from unittest import TestCase

import pandas as pd

from moonstone.plot.graphs.box import (
    GroupBoxGraph
)


class TestGroupBoxGraph(TestCase):

    def test_vertical_orientation(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M"],
                [3.0, "F"],
                [9.0, "M"],
                [6.0, "M"],
                [2.0, "F"]
            ],
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
            columns=["data", "sex"],
        )
        expected_x_gpM = ['M', 'M', 'M']
        expected_y_gpM = [1.0, 9.0, 6.0]
        expected_x_gpF = ['F', 'F']
        expected_y_gpF = [3.0, 2.0]
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="sex", show=False
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpM)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpM)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpF)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpF)

    def test_horizontal_orientation(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M"],
                [3.0, "F"],
                [9.0, "M"],
                [6.0, "M"],
                [2.0, "F"]
            ],
            index=["sample1", "sample2", "sample3", "sample4", "sample5"],
            columns=["data", "sex"],
        )
        expected_x_gpM = [1.0, 9.0, 6.0]
        expected_y_gpM = ['M', 'M', 'M']
        expected_x_gpF = [3.0, 2.0]
        expected_y_gpF = ['F', 'F']
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="sex", show=False,
            orientation="h"
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpM)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpM)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpF)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpF)

    def test_with_group_col2_without_groups(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M", "A"],
                [3.0, "F", "B"],
                [9.0, "M", "A"],
                [6.0, "M", "B"],
                [2.0, "F", "A"],
                [4.0, "F", "A"],
                [2.0, "M", "B"],
                [6.0, "M", "B"],
            ],
            index=[
                "sample1", "sample2", "sample3", "sample4", "sample5",
                "sample6", "sample7", "sample8"
            ],
            columns=["data", "sex", "group"],
        )
        expected_x_gpA = ['M', 'M', 'F', 'F']
        expected_y_gpA = [1.0, 9.0, 2.0, 4.0]
        expected_x_gpB = ['F', 'M', 'M', 'M']
        expected_y_gpB = [3.0, 6.0, 2.0, 6.0]
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="sex", group_col2="group", show=False,
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpA)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpA)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpB)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpB)

    def test_with_group_col2_with_groups_and_groups2(self):
        tested_df = pd.DataFrame(
            [
                [1.0, "M", "A"],
                [3.0, "F", "B"],
                [9.0, "M", "A"],
                [6.0, "M", "B"],
                [2.0, "F", "A"],
                [4.0, "F", "A"],
                [2.0, "M", "B"],
                [6.0, "M", "B"],
                [8.0, "M", "C"],
                [5.0, "F", "C"],
                [7.0, "M", "C"],
            ],
            index=[
                "sample1", "sample2", "sample3", "sample4", "sample5",
                "sample6", "sample7", "sample8", "sample9", "sample10",
                "sample11"
            ],
            columns=["data", "sex", "group"],
        )
        groups = ["F", "M"]    # change order
        groups2 = ["A", "B"]   # don't show group "C" (+ dictate order)
        expected_x_gpA = ['F', 'F', 'M', 'M']
        expected_y_gpA = [2.0, 4.0, 1.0, 9.0]
        expected_x_gpB = ['F', 'M', 'M', 'M']
        expected_y_gpB = [3.0, 6.0, 2.0, 6.0]
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="sex", group_col2="group", 
            groups=groups, groups2=groups2,
            show=False,
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpA)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpA)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpB)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpB)

    def test_sort_groups_with_group_col2(self):
        tested_df = pd.DataFrame(
            [
                [8.0, "M", "C"],   # "C" appears before so unsorted order is ["C", "A", "B"]
                [1.0, "M", "A"],   # and "M" appears before so unsorted order is ["M", "F"]
                [3.0, "F", "B"],
                [9.0, "M", "A"],
                [6.0, "M", "B"],
                [2.0, "F", "A"],
                [4.0, "F", "A"],
                [2.0, "M", "B"],
                [6.0, "M", "B"],
                [5.0, "F", "C"],
                [7.0, "M", "C"],
            ],
            index=[
                "sample1", "sample2", "sample3", "sample4", "sample5",
                "sample6", "sample7", "sample8", "sample9", "sample10",
                "sample11"
            ],
            columns=["data", "sex", "group"],
        )
        expected_x_gpA = ['F', 'F', 'M', 'M']
        expected_y_gpA = [2.0, 4.0, 1.0, 9.0]
        expected_x_gpB = ['F', 'M', 'M', 'M']
        expected_y_gpB = [3.0, 6.0, 2.0, 6.0]
        expected_x_gpC = ['F', 'M', 'M']
        expected_y_gpC = [5., 8., 7.]
        plot = GroupBoxGraph(tested_df)
        tested_graph = plot.plot_one_graph(
            data_col="data", group_col="sex", group_col2="group", 
            sort_groups=True,
            show=True,
            )
        np.testing.assert_array_equal(tested_graph.data[0].x, expected_x_gpA)
        np.testing.assert_array_equal(tested_graph.data[0].y, expected_y_gpA)
        np.testing.assert_array_equal(tested_graph.data[1].x, expected_x_gpB)
        np.testing.assert_array_equal(tested_graph.data[1].y, expected_y_gpB)
        np.testing.assert_array_equal(tested_graph.data[2].x, expected_x_gpC)
        np.testing.assert_array_equal(tested_graph.data[2].y, expected_y_gpC)