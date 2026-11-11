import os
import tempfile
from unittest import TestCase

import pandas as pd
from pandas.testing import assert_frame_equal

from moonstone.parsers.handleCounts import Inputs


class TestHandleCounts(TestCase):
    """Count file read by the command line: one row per taxon, one column per sample."""

    def write(self, text):
        path = os.path.join(tempfile.mkdtemp(), 'counts.csv')
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_samples_become_rows(self):
        path = self.write('OTU,s1,s2\ntaxonA,1,2\ntaxonB,3,4\n')
        expected = pd.DataFrame({'taxonA': [1.0, 2.0], 'taxonB': [3.0, 4.0]},
                                index=pd.Index(['s1', 's2'], name='sample'))
        expected.columns.name = ''
        assert_frame_equal(Inputs(path).opencounts(), expected, check_index_type=False, check_column_type=False)

    def test_unnamed_taxa_removed_and_counted(self):
        path = self.write('OTU,s1,s2\ntaxonA,1,2\n,5,6\ntaxonB,3,4\n,7,8\n')
        with self.assertLogs('moonstone.parsers.handleCounts', level='WARNING') as logs:
            df = Inputs(path).opencounts()
        self.assertEqual(list(df.columns), ['taxonA', 'taxonB'])
        self.assertEqual(list(df['taxonB']), [3.0, 4.0])
        self.assertIn('Removed 2 unnamed/unidentified variable(s)!', logs.output[0])
