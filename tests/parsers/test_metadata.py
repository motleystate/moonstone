import os
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from metautils.parsers.metadata import MetadataParser


class TestPicrust2PathwaysParser(TestCase):

    def setUp(self):
        self.test_file = os.path.join(os.path.dirname(__file__), "metadata.tsv")

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {
                'col_1': ['s1', 's2', 's3'],
                'col_2': [13.3, 15.3, 19.1],
                'col_3': ['M', 'F', 'M'],
                'col_4': ['T', 'F', 'F']
            }
        )
        parser = MetadataParser(self.test_file)
        assert_frame_equal(parser.dataframe, expected_df)
