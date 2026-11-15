import os
import tempfile
import time
from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone import examples
from moonstone.ui.catalog import TOOLS
from moonstone.ui.checks import check_inputs, metadata_columns, summary
from moonstone.validation import has_errors

METADATA = examples.path('metadata.csv')


class TestChecks(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def write(self, name, text):
        path = os.path.join(self.tmp, name)
        with open(path, 'w') as handle:
            handle.write(text)
        return path

    def test_examples_have_no_error(self):
        for tool in TOOLS.values():
            meta = METADATA if tool.needs_metadata else None
            findings = check_inputs(tool, examples.path(tool.example['countfile']), meta)
            self.assertFalse(has_errors(findings), (tool.key, findings))

    def test_missing_file(self):
        codes = [f.code for f in check_inputs(TOOLS['permanova'], '/nope.tsv', METADATA)]
        self.assertEqual(codes, ['file_missing'])

    def test_wrong_table_for_the_tool(self):
        findings = check_inputs(TOOLS['permanova'], examples.path('resistance_genes.csv'), METADATA)
        unreadable = [f for f in findings if f.code == 'not_readable']
        self.assertIn('no taxonomy column found', unreadable[0].message)

    def test_tab_separated_table_for_classification(self):
        codes = [f.code for f in check_inputs(TOOLS['classify'], examples.path('species_counts.tsv'), METADATA)]
        self.assertIn('separator_mismatch', codes)

    def test_french_locale_file(self):
        path = self.write('fr.csv', 'gene;S001;S002\ntetM;1,5;2,25\n')
        findings = check_inputs(TOOLS['correlation'], path, METADATA)
        self.assertIn('decimal_comma', [f.code for f in findings])
        self.assertTrue(all('fr.csv' in f.message for f in findings if f.code == 'decimal_comma'))

    def test_second_count_table_is_checked(self):
        findings = check_inputs(TOOLS['correlation'], examples.path('species_counts.tsv'), METADATA, '/nope.csv')
        self.assertEqual([f.code for f in findings if f.level == 'error'], ['file_missing'])

    def test_numeric_looking_sample_ids_stay_text(self):
        counts = self.write('c.csv', 'gene,01,02,03\ntetM,1,2,3\nermB,4,5,6\n')
        meta = self.write('m.csv', 'sample,group\n1,a\n2,b\n3,a\n')
        self.assertTrue(has_errors(check_inputs(TOOLS['correlation'], counts, meta)))
        meta = self.write('m2.csv', 'sample,group\n01,a\n02,b\n03,a\n')
        self.assertFalse(has_errors(check_inputs(TOOLS['correlation'], counts, meta)))

    def test_nothing_is_written_next_to_the_inputs(self):
        counts = self.write('c.csv', 'gene,s1,s2,s3\ntetM,1,2,3\n')
        meta = self.write('m.csv', 'sample,group\ns1,a\ns2,b\ns3,a\n')
        check_inputs(TOOLS['correlation'], counts, meta)
        self.assertEqual(sorted(os.listdir(self.tmp)), ['c.csv', 'm.csv'])

    def test_metadata_columns(self):
        columns = metadata_columns(METADATA)
        self.assertEqual(columns['Type'], ['delivery', 'inclusion'])
        self.assertNotIn('sample', columns)
        self.assertEqual(metadata_columns('/nope.csv'), {})

    def test_summary(self):
        self.assertEqual(summary(TOOLS['permanova'], examples.path('species_counts.tsv')), '60 samples')
        self.assertEqual(summary(TOOLS['permanova'], '/nope.tsv'), '')

    def test_large_table_is_checked_in_reasonable_time(self):
        rng = np.random.default_rng(0)
        table = pd.DataFrame(rng.poisson(5, (20000, 500)), columns=[f's{i}' for i in range(500)])
        table.insert(0, 'species', [f'species {i}' for i in range(20000)])
        path = os.path.join(self.tmp, 'large.tsv')
        table.to_csv(path, sep='\t', index=False)
        meta = self.write('m.csv', 'sample,group\n' + ''.join(f's{i},{i % 2}\n' for i in range(500)))
        start = time.time()
        findings = check_inputs(TOOLS['permanova'], path, meta)
        self.assertLess(time.time() - start, 60)
        self.assertFalse(has_errors(findings), findings)
