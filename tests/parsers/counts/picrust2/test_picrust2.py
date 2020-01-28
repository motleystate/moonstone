import os
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from moonstone.parsers.counts.picrust2 import Picrust2PathwaysParser


class TestPicrust2PathwaysParser(TestCase):

    def setUp(self):
        self.test_file = os.path.join(os.path.dirname(__file__), "picrust2.tsv")

    def test_parse_file(self):
        expected_df = pd.DataFrame(
            {
                'sample_1': [14.3, 94.1, 9.4],
                'sample_2': [123.4, 1231.1, 15.5],
                'sample_3': [12.0, 124.2, 12.2]
            }
        )
        expected_df.index = ['pathway_1', 'pathway_2', 'pathway_4']
        expected_df.index.name = 'pathways'
        parser = Picrust2PathwaysParser(self.test_file)
        assert_frame_equal(parser.dataframe, expected_df)
