import yaml
from collections import defaultdict
from typing import Dict, List

from pandas import DataFrame

from moonstone.analysis.columns_statistics import DataframeStatistics
from moonstone.parsers.transform.cleaning import DataFrameCleaner
from .base import BaseParser


class MetadataParser(BaseParser):

    def __init__(self, *args, cleaning_operations: dict = None, **kwargs):
        """
        Cleaning operations are based on DataFrameCleaner object that allow to perform transformation
        operations on different columns.

        Format is the following:

            {'col_name': [('operation1', 'operation1_options'), ('operation2', 'operation2_options')]}

        :param cleaning_operations: cleaning operations to apply to the input table
        """
        self.cleaning_operations = cleaning_operations
        if self.cleaning_operations is None:
            self.cleaning_operations = {}
        super().__init__(*args, **kwargs)

    def to_dataframe(self) -> DataFrame:
        dataframe = super().to_dataframe()
        df_cleaner = DataFrameCleaner(dataframe)
        for col_name, transformations in self.cleaning_operations.items():
            for transformation in transformations:
                transf_name = transformation[0]
                transf_options = transformation[1]
                df_cleaner.run_transform(col_name, transf_name, transf_options)
        return df_cleaner.df

    def get_stats(self) -> List[Dict]:
        """
        :return: list of dict containing statistics about each column
        """
        return DataframeStatistics(self.dataframe).get_stats()


class YAMLBasedMetadataParser:

    def __init__(self, metadata_file_path, config_file_path, **kwargs):
        self._parse_yaml_config(config_file_path)
        self.metadata_parser = MetadataParser(
            metadata_file_path, cleaning_operations=self.cleaning_operations,
            parsing_options=self.parsing_options, **kwargs
        )

    def _parse_yaml_config(self, config_file_path):
        with open(config_file_path, 'r') as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
        self.parsing_options = self._extract_parsing_options(config['parsing'])
        self.cleaning_operations = self._extract_cleaning_operations(config['parsing'])

    def _extract_parsing_options(self, parsing_config):
        parsing_options = defaultdict(lambda: {})
        for col_parsing_config in parsing_config:
            if 'dtype' in col_parsing_config.keys():
                parsing_options['dtype'][col_parsing_config['col_name']] = col_parsing_config['dtype']
        return parsing_options

    def _extract_cleaning_operations(self, parsing_config):
        cleaning_operations = defaultdict(lambda: [])
        for col_parsing_config in parsing_config:
            if 'operations' in col_parsing_config.keys():
                for operation in col_parsing_config['operations']:
                    cleaning_operations[col_parsing_config['col_name']].append(
                        (operation['name'], operation.get('options', {}))
                    )
        return cleaning_operations
