from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.plot.graphs import ParallelCategoriesGraph, SankeyDiagram
from moonstone.plot.graphs.category_flows import MISSING_LABEL, prepare_categories


class TestPrepareCategories(TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {"age": ["young", "young", "old", "old", "old"], "city": ["Paris", "Lyon", "Paris", "Paris", np.nan]},
            index=list("abcde"),
        )

    def test_missing_shown_as_category_not_dropped_or_filled(self):
        with self.assertLogs("moonstone.plot.graphs.category_flows", level="WARNING") as logs:
            str_df, weights, cat_order, _ = prepare_categories(self.df)
        self.assertIn("1 missing", logs.output[0])
        self.assertEqual(str_df.shape, self.df.shape)
        self.assertEqual(str_df.loc["e", "city"], MISSING_LABEL)
        self.assertTrue(self.df.loc["e", "city"] != self.df.loc["e", "city"])  # input untouched (still NaN)
        self.assertEqual(list(weights), [1.0] * 5)

    def test_default_order_largest_first_then_name(self):
        _, _, cat_order, _ = prepare_categories(self.df)
        self.assertEqual(cat_order["age"], ["old", "young"])
        self.assertEqual(cat_order["city"], ["Paris", MISSING_LABEL, "Lyon"])

    def test_alphabetical_order(self):
        _, _, cat_order, _ = prepare_categories(self.df, sort_categories="alphabetical")
        self.assertEqual(cat_order["age"], ["old", "young"])
        self.assertEqual(cat_order["city"], ["(missing)", "Lyon", "Paris"])
        with self.assertRaises(ValueError):
            prepare_categories(self.df, sort_categories="random")

    def test_order_given_then_rest(self):
        _, _, cat_order, _ = prepare_categories(self.df, order={"city": ["Lyon"]})
        self.assertEqual(cat_order["city"], ["Lyon", "Paris", MISSING_LABEL])

    def test_order_with_unknown_category_warns(self):
        with self.assertLogs("moonstone.plot.graphs.category_flows", level="WARNING"):
            prepare_categories(self.df.dropna(), order={"city": ["Rome"]})

    def test_column_titles(self):
        _, _, _, titles = prepare_categories(self.df, column_titles={"age": "Age group"})
        self.assertEqual(titles, ["Age group", "city"])

    def test_numbers_are_categories_in_different_columns_do_not_clash(self):
        df = pd.DataFrame({"a": [1, 1, 2], "b": [1, 2, 2]})
        str_df, _, _, _ = prepare_categories(df)
        self.assertEqual(list(str_df["a"]), ["1", "1", "2"])

    def test_errors(self):
        with self.assertRaises(ValueError):
            prepare_categories(self.df[["age"]])
        with self.assertRaises(ValueError):
            prepare_categories(self.df, order={"nope": []})
        with self.assertRaises(ValueError):
            prepare_categories(self.df, weights=pd.Series([1, 1, 1, 1, 1]))  # wrong index
        with self.assertRaises(ValueError):
            prepare_categories(self.df, weights=pd.Series([1, 1, -1, 1, 1], index=self.df.index))


class TestSankeyDiagram(TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {
                "age": ["18-30", "18-30", "31-50", "31-50", "51+", "18-30"],
                "city": ["Paris", "Lyon", "Paris", "Paris", "Lyon", "Paris"],
                "smoker": ["yes", "no", "no", "yes", "no", "no"],
            },
            index=[f"s{i}" for i in range(6)],
        )

    @staticmethod
    def _links(fig):
        trace = fig.data[0]
        labels = [lab.split("<br>")[0] for lab in trace.node.label]
        return {(labels[s], labels[t]): v for s, t, v in zip(trace.link.source, trace.link.target, trace.link.value)}

    def test_links_count_samples_between_neighbouring_columns(self):
        fig = SankeyDiagram(self.df).plot_one_graph(show=False)
        links = self._links(fig)
        self.assertEqual(links[("18-30", "Paris")], 2)
        self.assertEqual(links[("18-30", "Lyon")], 1)
        self.assertEqual(links[("Paris", "no")], 2)
        self.assertEqual(sum(v for (a, _), v in links.items() if a in {"18-30", "31-50", "51+"}), 6)
        self.assertEqual(sum(v for (a, _), v in links.items() if a in {"Paris", "Lyon"}), 6)

    def test_same_value_in_two_columns_gives_two_nodes(self):
        df = pd.DataFrame({"a": ["x", "y", "x"], "b": ["x", "x", "y"]})
        fig = SankeyDiagram(df).plot_one_graph(show=False)
        self.assertEqual(len(fig.data[0].node.label), 4)

    def test_node_labels_show_counts(self):
        fig = SankeyDiagram(self.df).plot_one_graph(show=False)
        self.assertIn("18-30<br>3", fig.data[0].node.label)

    def test_order_gives_node_positions_top_to_bottom(self):
        fig = SankeyDiagram(self.df).plot_one_graph(show=False, order={"age": ["51+", "31-50", "18-30"]})
        trace = fig.data[0]
        y = {lab.split("<br>")[0]: yy for lab, yy, xx in zip(trace.node.label, trace.node.y, trace.node.x) if xx < 0.1}
        self.assertLess(y["51+"], y["31-50"])
        self.assertLess(y["31-50"], y["18-30"])

    def test_column_titles_added_as_annotations(self):
        fig = SankeyDiagram(self.df).plot_one_graph(show=False, column_titles={"age": "Age group"})
        self.assertEqual([a.text for a in fig.layout.annotations], ["<b>Age group</b>", "<b>city</b>", "<b>smoker</b>"])

    def test_missing_values_do_not_crash(self):
        df = self.df.copy()
        df.loc["s1", "city"] = np.nan
        with self.assertLogs("moonstone.plot.graphs.category_flows", level="WARNING"):
            fig = SankeyDiagram(df).plot_one_graph(show=False)
        self.assertIn(("18-30", MISSING_LABEL), self._links(fig))

    def test_weights(self):
        weights = pd.Series([1, 1, 1, 1, 1, 5.0], index=self.df.index)
        fig = SankeyDiagram(self.df, weights=weights).plot_one_graph(show=False)
        self.assertEqual(self._links(fig)[("18-30", "Paris")], 6)

    def test_plotting_options_and_output_file(self):
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "sankey.html")
            fig = SankeyDiagram(self.df).plot_one_graph(
                show=False, plotting_options={"layout": {"title_text": "T"}}, output_file=path)
            self.assertEqual(fig.layout.title.text, "T")
            self.assertTrue(os.path.getsize(path) > 0)


class TestParallelCategoriesGraph(TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {"a": ["x", "x", "y", "y"], "b": ["u", "v", "u", "u"], "c": [1, 2, 1, 2], "grp": ["g1", "g2", "g1", "g1"]},
            index=list("wxyz"),
        )

    def test_dimensions_and_counts(self):
        fig = ParallelCategoriesGraph(self.df[["a", "b", "c"]]).plot_one_graph(
            show=False, column_titles={"a": "A"}, order={"a": ["y", "x"]})
        trace = fig.data[0]
        self.assertEqual([d.label for d in trace.dimensions], ["A", "b", "c"])
        self.assertEqual(list(trace.dimensions[0].categoryarray), ["y", "x"])
        self.assertEqual(list(trace.counts), [1, 1, 1, 1])

    def test_color_by_column_not_plotted(self):
        fig = ParallelCategoriesGraph(self.df[["a", "b"]]).plot_one_graph(show=False, color_by=self.df["grp"])
        line = fig.data[0].line
        # g1 -> first step, g2 -> second step of a two-step colourscale
        self.assertEqual(list(line.color), [0.5, 1.5, 0.5, 0.5])
        self.assertEqual(len(line.colorscale), 4)

    def test_color_by_series_wrong_index(self):
        with self.assertRaises(ValueError):
            ParallelCategoriesGraph(self.df[["a", "b"]]).plot_one_graph(
                show=False, color_by=pd.Series(["g"] * 4))

    def test_custom_colors(self):
        fig = ParallelCategoriesGraph(self.df[["a", "b"]]).plot_one_graph(
            show=False, color_by=self.df["grp"], colors={"g2": "#ff0000"})
        self.assertIn("#ff0000", [step[1] for step in fig.data[0].line.colorscale])

    def test_colorscale_spreads_categories_along_scale_in_given_order(self):
        df = pd.DataFrame(
            {"a": ["x", "y", "x", "y"], "b": ["u", "u", "v", "v"], "age": ["<=30", ">40", "31-40", ">40"]})
        fig = ParallelCategoriesGraph(df[["a", "b"]]).plot_one_graph(
            show=False, color_by=df["age"].rename("age"), colorscale=[[0, "rgb(0,0,0)"], [1, "rgb(255,255,255)"]],
            order={"age": ["<=30", "31-40", ">40"]})
        line = fig.data[0].line
        steps = [step[1] for step in line.colorscale][::2]
        self.assertEqual(steps, ["rgb(0, 0, 0)", "rgb(128, 128, 128)", "rgb(255, 255, 255)"])
        # rows: <=30 -> first step, >40 -> last step, 31-40 -> middle step
        self.assertEqual(list(line.color), [0.5, 2.5, 1.5, 2.5])

    def test_colorscale_single_category_and_colors_override(self):
        df = pd.DataFrame({"a": ["x", "y"], "b": ["u", "v"], "g": ["only", "only"]})
        fig = ParallelCategoriesGraph(df[["a", "b"]]).plot_one_graph(show=False, color_by=df["g"], colorscale="Viridis")
        self.assertEqual(len(fig.data[0].line.colorscale), 2)
        fig = ParallelCategoriesGraph(df[["a", "b"]]).plot_one_graph(
            show=False, color_by=df["g"], colorscale="Viridis", colors={"only": "#123456"})
        self.assertEqual(fig.data[0].line.colorscale[0][1], "#123456")
