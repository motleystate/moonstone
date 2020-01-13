import os
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from moonstone.parsers.base import BaseParser


class TestBaseParser(TestCase):

    def setUp(self):
        self.test_file = os.path.join(os.path.dirname(__file__), "data/metadata/base.csv")
        self.parser = BaseParser(self.test_file, sep=',')

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {
                'col_1': ['s1', 's2', 's3'],
                'col_2': [13.3, 15.3, 19.1],
                'col_3': ['M', 'F', 'M']
            }
        )
        assert_frame_equal(self.parser.dataframe, expected_df)

    def test_access_df_twice(self):
        """make sure the property works"""
        self.assertIsInstance(self.parser.dataframe, pd.DataFrame)
        self.assertIsInstance(self.parser.dataframe, pd.DataFrame)
