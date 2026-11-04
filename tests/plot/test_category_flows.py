from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.plot.category_flows import plot_category_flows, plot_taxonomy_flows


class TestPlotCategoryFlows(TestCase):
    def setUp(self):
        self.df = pd.DataFrame(
            {"a": ["x", "x", "y"], "b": ["u", "v", "u"], "c": ["p", "p", "q"]}, index=["s1", "s2", "s3"])

    def test_kinds(self):
        sankey = plot_category_flows(self.df, show=False)
        self.assertEqual(sankey.data[0].type, "sankey")
        parcats = plot_category_flows(self.df, kind="parcats", show=False)
        self.assertEqual(parcats.data[0].type, "parcats")

    def test_columns_selection(self):
        fig = plot_category_flows(self.df, columns=["a", "c"], kind="parcats", show=False)
        self.assertEqual([d.label for d in fig.data[0].dimensions], ["a", "c"])

    def test_color_by_column_outside_plotted_columns(self):
        fig = plot_category_flows(self.df, columns=["a", "c"], kind="parcats", color_by="b", show=False)
        self.assertEqual(len(fig.data[0].line.colorscale), 4)

    def test_bad_kind_and_color_by_with_sankey(self):
        with self.assertRaises(ValueError):
            plot_category_flows(self.df, kind="pie", show=False)
        with self.assertRaises(ValueError):
            plot_category_flows(self.df, color_by="b", show=False)
        with self.assertRaises(ValueError):
            plot_category_flows(self.df, colorscale="Viridis", show=False)

    def test_sort_categories_alphabetical(self):
        fig = plot_category_flows(self.df, kind="parcats", sort_categories="alphabetical", show=False)
        self.assertEqual(list(fig.data[0].dimensions[0].categoryarray), ["x", "y"])
        fig = plot_taxonomy_flows(
            pd.DataFrame({"s": [1, 3]}, index=pd.MultiIndex.from_tuples([("g1", "a"), ("g2", "b")], names=["g", "sp"])),
            ["g", "sp"], kind="parcats", sort_categories="alphabetical", show=False)
        self.assertEqual(list(fig.data[0].dimensions[0].categoryarray), ["g1", "g2"])


class TestPlotTaxonomyFlows(TestCase):
    def setUp(self):
        index = pd.MultiIndex.from_tuples(
            [
                ("F1", "G1", "S1"), ("F1", "G1", "S2"), ("F1", "G2", "S3"),
                ("F2", "G3", "S4"), ("F2", "G4", "S5"), ("F3", "G5", "S6"),
            ],
            names=["family", "genus", "species"],
        )
        self.counts = pd.DataFrame(
            {"sample1": [10, 10, 20, 30, 20, 10], "sample2": [0, 20, 20, 20, 20, 20]}, index=index)

    @staticmethod
    def _links(fig):
        trace = fig.data[0]
        labels = [lab.split("<br>")[0] for lab in trace.node.label]
        return {(labels[s], labels[t]): v for s, t, v in zip(trace.link.source, trace.link.target, trace.link.value)}

    def test_weights_are_mean_relative_abundance(self):
        fig = plot_taxonomy_flows(self.counts, ["species", "genus", "family"], top_n=None, show=False)
        links = self._links(fig)
        # S1: 10% in sample1, 0% in sample2 -> mean 5; S2: 10% and 20% -> mean 15
        self.assertAlmostEqual(links[("S1", "G1")], 5)
        self.assertAlmostEqual(links[("S2", "G1")], 15)
        self.assertAlmostEqual(links[("G1", "F1")], 20)
        self.assertAlmostEqual(sum(v for (a, _), v in links.items() if a.startswith("S")), 100)

    def test_total_reads(self):
        fig = plot_taxonomy_flows(
            self.counts, ["genus", "family"], weight="total_reads", top_n=None, show=False)
        self.assertEqual(self._links(fig)[("G1", "F1")], 40)

    def test_top_n_lumps_into_other(self):
        fig = plot_taxonomy_flows(self.counts, ["species", "genus"], top_n=2, show=False)
        labels = {lab.split("<br>")[0] for lab in fig.data[0].node.label}
        self.assertIn("Other", labels)
        self.assertLessEqual(len([x for x in labels if x.startswith("S")]), 2)
        # nothing is lost: all of the weight is still drawn
        trace = fig.data[0]
        self.assertAlmostEqual(sum(trace.link.value), 100)

    def test_unclassified_taxon_is_shown_as_missing(self):
        index = self.counts.index.to_frame(index=False)
        index.loc[5, "genus"] = np.nan
        counts = self.counts.copy()
        counts.index = pd.MultiIndex.from_frame(index)
        with self.assertLogs("moonstone.plot.graphs.category_flows", level="WARNING"):
            fig = plot_taxonomy_flows(counts, ["species", "genus"], top_n=None, show=False)
        self.assertIn(("S6", "(missing)"), self._links(fig))

    def test_parcats_kind(self):
        fig = plot_taxonomy_flows(self.counts, ["species", "genus"], kind="parcats", show=False)
        self.assertAlmostEqual(sum(fig.data[0].counts), 100)

    def test_errors(self):
        with self.assertRaises(ValueError):
            plot_taxonomy_flows(self.counts, ["species", "order"], show=False)
        with self.assertRaises(ValueError):
            plot_taxonomy_flows(self.counts, ["species", "genus"], weight="nope", show=False)
