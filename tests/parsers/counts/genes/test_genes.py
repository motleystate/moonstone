import os
import tempfile
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from moonstone.parsers.counts.genes import GeneCountsParser


class TestGeneCountsParser(TestCase):

    def test_parse_file(self):
        path = os.path.join(tempfile.mkdtemp(), 'genes.tsv')
        with open(path, 'w') as handle:
            handle.write('genes\tsample_1\tsample_2\ngene_1\t3\t19\ngene_2\t9\t10\n')
        expected_df = pd.DataFrame({'sample_1': [3, 9], 'sample_2': [19, 10]},
                                   index=pd.Index(['gene_1', 'gene_2'], name='genes'))
        assert_frame_equal(GeneCountsParser(path).dataframe, expected_df)
