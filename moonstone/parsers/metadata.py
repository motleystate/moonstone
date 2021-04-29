"""Classes to handle metadata import."""
import yaml
from collections import defaultdict
from typing import Dict, List

from pandas import DataFrame, Series
import plotly.io
import plotly.graph_objects as go

from moonstone.analysis.columns_statistics import DataframeStatistics
from moonstone.parsers.transform.cleaning import DataFrameCleaner
from .base import BaseParser


class MetadataParser(BaseParser):
    """Parse metadata file and allows to apply transformations on them (cleaning...)."""

    DEFAULT_COLORSCALE = [
        [0, "rgb(166,206,227)"],
        [0.25, "rgb(31,120,180)"],
        [0.45, "rgb(178,223,138)"],
        [0.65, "rgb(51,160,44)"],
        [0.85, "rgb(251,154,153)"],
        [1, "rgb(227,26,28)"],
    ]

    def __init__(
        self,
        *args,
        index_col: str = "sample",
        cleaning_operations: dict = None,
        **kwargs
    ):
        """
        Parse metadata file and allows to apply transformations on them (cleaning...).

        Cleaning operations are based on DataFrameCleaner object that allows to perform transformation
        operations on different columns.

        Format is the following:

        .. code-block:: python

            {'col_name': [('operation1', 'operation1_options'), ('operation2', 'operation2_options')]}

        Args:
            index_col: name of the column used as dataframe index
            cleaning_operations: cleaning operations to apply to the input table
        """
        self.index_col = index_col
        self.cleaning_operations = (
            {} if cleaning_operations is None else cleaning_operations
        )
        super().__init__(*args, **kwargs)

    def _load_data(self) -> DataFrame:
        """Load and return Pandas dataframe."""
        dataframe = super()._load_data()
        df_cleaner = DataFrameCleaner(dataframe)
        for col_name, transformations in self.cleaning_operations.items():
            for transformation in transformations:
                transf_name = transformation[0]
                transf_options = transformation[1]
                df_cleaner.run_transform(col_name, transf_name, transf_options)
        return df_cleaner.df.set_index(self.index_col)

    def get_stats(self) -> List[Dict]:
        """
        Retrieve statistics about each columns.

        Returns:
            list of dict containing statistics about each column
        """
        return DataframeStatistics(self.dataframe).get_stats()

    def _get_dimensions(self, categories) -> List[go.parcats.Dimension]:
        dimensions = []
        for cat in categories:
            dimensions.append(
                go.parcats.Dimension(
                    values=self.dataframe[cat],
                    categoryorder="category ascending",
                    label=cat,
                )
            )
        return dimensions

    def _get_color(self, color_by: str) -> Series:
        cpt = 1
        color_dict = {}
        for i in self.dataframe[color_by].unique():
            color_dict[i] = cpt
            cpt += 1
        return self.dataframe[color_by].apply(lambda x: color_dict.get(x, 0))

    def visualize_categories(
        self,
        categories: list,
        color_by: str,
        colorscale: list = None,
        title: str = "Metadata categories distribution",
        output_file: str = "",
    ):
        """
        Visualize category metadata with parallel categories diagram.

        :param categories: list of column to display
        :param color_by: perform coloration on the given category
        """
        if colorscale is None:
            colorscale = self.DEFAULT_COLORSCALE
        dimensions = self._get_dimensions(categories)
        color = self._get_color(color_by)

        fig = go.Figure(
            data=[
                go.Parcats(
                    dimensions=dimensions,
                    line={"color": color, "colorscale": colorscale},
                    hoverinfo="count",
                    arrangement="freeform",
                )
            ]
        )
        fig.update_layout(title=title)
        fig.show()
        if output_file:
            plotly.io.write_html(fig, output_file)


class YAMLBasedMetadataParser:
    """Metadata Parser with operations configured in a YAML file."""

    def __init__(self, metadata_file_path, config_file_path, **kwargs):
        """Metadata Parser with operations configured in a YAML file."""
        self._parse_yaml_config(config_file_path)
        self.metadata_parser = MetadataParser(
            metadata_file_path,
            cleaning_operations=self.cleaning_operations,
            parsing_options=self.parsing_options,
            **kwargs
        )

    def _parse_yaml_config(self, config_file_path):
        with open(config_file_path, "r") as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
        self.parsing_options = self._extract_parsing_options(config["parsing"])
        self.cleaning_operations = self._extract_cleaning_operations(config["parsing"])

    def _extract_parsing_options(self, parsing_config):
        parsing_options = defaultdict(lambda: {})
        for col_parsing_config in parsing_config:
            if "dtype" in col_parsing_config.keys():
                parsing_options["dtype"][
                    col_parsing_config["col_name"]
                ] = col_parsing_config["dtype"]
        return parsing_options

    def _extract_cleaning_operations(self, parsing_config):
        cleaning_operations = defaultdict(lambda: [])
        for col_parsing_config in parsing_config:
            if "operations" in col_parsing_config.keys():
                for operation in col_parsing_config["operations"]:
                    cleaning_operations[col_parsing_config["col_name"]].append(
                        (operation["name"], operation.get("options", {}))
                    )
        return cleaning_operations
