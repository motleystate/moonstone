import os
from unittest import TestCase

import pandas as pd
import plotly.graph_objects as go
from pandas.testing import assert_frame_equal

from moonstone.parsers.metadata import MetadataParser, YAMLBasedMetadataParser


class TestMetadataParser(TestCase):
    def setUp(self):
        self.metadata_file = os.path.join(
            os.path.dirname(__file__), "data/metadata/metadata.tsv"
        )
        self.metadata_file_no_header = os.path.join(
            os.path.dirname(__file__), "data/metadata/metadata_noheader.tsv"
        )
        self.metadata_file_dirty = os.path.join(
            os.path.dirname(__file__), "data/metadata/dirty_metadata.tsv"
        )

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {"col_2": [13.3, 15.3, 19.1], "col_3": ["M", "F", "M"]}
        )
        expected_df.index = pd.Index(["s1", "s2", "s3"], name="col_1")
        parser = MetadataParser(self.metadata_file, index_col="col_1")
        assert_frame_equal(parser.dataframe, expected_df)

    def test_get_stats_headers(self):
        expected_list = [
            {
                "col_name": "col_2",
                "col_type": "float64",
                "python_col_type": "float",
                "n_values": 3,
                "n_uniq_values": 3,
                "mean": 15.9,
                "values_repartition": {13.3: 1, 15.3: 1, 19.1: 1},
            },
            {
                "col_name": "col_3",
                "col_type": "object",
                "python_col_type": "str",
                "n_values": 3,
                "n_uniq_values": 2,
                "values_repartition": {"M": 2, "F": 1},
            },
        ]
        parser = MetadataParser(self.metadata_file, index_col="col_1")
        self.assertListEqual(parser.get_stats(), expected_list)

    def test_get_stats_no_header(self):
        expected_list = [
            {
                "col_name": 1,
                "col_type": "float64",
                "python_col_type": "float",
                "n_values": 3,
                "n_uniq_values": 3,
                "mean": 15.9,
                "values_repartition": {13.3: 1, 15.3: 1, 19.1: 1},
            },
            {
                "col_name": 2,
                "col_type": "object",
                "python_col_type": "str",
                "n_values": 3,
                "n_uniq_values": 2,
                "values_repartition": {"M": 2, "F": 1},
            },
        ]
        parser = MetadataParser(
            self.metadata_file_no_header, no_header=True, index_col=0
        )
        self.assertListEqual(parser.get_stats(), expected_list)

    def test_parse_file_force_dtype(self):
        expected_df = pd.DataFrame(
            {"col_2": ["13.3", "15.3", "19.1"], "col_3": ["M", "F", "M"]}
        )
        expected_df.index = pd.Index(["s1", "s2", "s3"], name="col_1")
        parsing_options = {"dtype": {"col_2": "object"}}
        parser = MetadataParser(
            self.metadata_file, parsing_options=parsing_options, index_col="col_1"
        )
        assert_frame_equal(parser.dataframe, expected_df)

    def test_parse_dirty_metadata_and_clean(self):
        expected_df = pd.DataFrame(
            {
                "age": [29, 48, 36, 25],
            }
        )
        expected_df.index = pd.Index(["s1", "s2", "s3", "s4"], name="sample")
        cleaning_operations = {
            "samples": [("to_slug", {}), ("rename", {"new_name": "sample"})]
        }
        parser = MetadataParser(
            self.metadata_file_dirty,
            sep=",",
            cleaning_operations=cleaning_operations,
            index_col="sample",
        )
        assert_frame_equal(parser.dataframe, expected_df)

    def test_get_dimensions(self):
        categories = ["col_3"]
        parser = MetadataParser(self.metadata_file, index_col="col_1")
        output_dimensions = parser._get_dimensions(categories)
        self.assertEqual(len(output_dimensions), 1)
        for dim in output_dimensions:
            self.assertTrue(isinstance(dim, go.parcats.Dimension))

    def test_get_color(self):
        color_by = "col_3"
        index = pd.Series(("s1", "s2", "s3"), name="col_1")
        expected_series = pd.Series([1, 2, 1], index=index, name=color_by)
        parser = MetadataParser(self.metadata_file, index_col="col_1")
        pd.testing.assert_series_equal(parser._get_color(color_by), expected_series)


class MockedYAMLBasedMetadataParser(YAMLBasedMetadataParser):
    """
    Mocked to skip __init__ and test only private methods of the class
    """

    def __init__(self):
        pass


class TestYAMLBasedMetadataParser(TestCase):
    def setUp(self):
        # For unit tests
        self.parsing_config = [
            {"col_name": "col_1", "dtype": "object"},
            {"col_name": "col_2", "operations": [{"name": "to_slug"}]},
            {
                "col_name": "col_3",
                "operations": [{"name": "rename", "options": {"new_name": "new"}}],
            },
        ]

    def test_extract_parsing_options(self):
        expected_dict = {"dtype": {"col_1": "object"}}
        parser = MockedYAMLBasedMetadataParser()
        self.assertDictEqual(
            parser._extract_parsing_options(self.parsing_config), expected_dict
        )

    def test_extract_cleaning_operations(self):
        expected_dict = {
            "col_2": [("to_slug", {})],
            "col_3": [("rename", {"new_name": "new"})],
        }
        parser = MockedYAMLBasedMetadataParser()
        self.assertDictEqual(
            parser._extract_cleaning_operations(self.parsing_config), expected_dict
        )

    def test_parse_yaml_config(self):
        config_file = os.path.join(
            os.path.dirname(__file__), "data/metadata/config.yaml"
        )
        expected_parsing_options = {"dtype": {"age": "object"}}
        expected_cleaning_operations = {
            "samples": [("to_slug", {}), ("rename", {"new_name": "sample"})],
        }
        parser = MockedYAMLBasedMetadataParser()
        parser._parse_yaml_config(config_file)
        self.assertDictEqual(parser.parsing_options, expected_parsing_options)
        self.assertDictEqual(parser.cleaning_operations, expected_cleaning_operations)

    def test_parse_end_to_end(self):
        metadata_file_dirty = os.path.join(
            os.path.dirname(__file__), "data/metadata/dirty_metadata.tsv"
        )
        config_file = os.path.join(
            os.path.dirname(__file__), "data/metadata/config.yaml"
        )
        parser = YAMLBasedMetadataParser(
            metadata_file_dirty, config_file, sep=",", index_col="sample"
        )
        expected_df = pd.DataFrame(
            {
                "age": ["29", "48", "36", "25"],
            }
        )
        expected_df.index = pd.Index(["s1", "s2", "s3", "s4"], name="sample")
        pd.testing.assert_frame_equal(parser.metadata_parser.dataframe, expected_df)
