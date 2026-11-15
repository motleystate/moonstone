import contextlib
import os
import tempfile
from unittest import TestCase

import matplotlib

from moonstone import examples
from moonstone.filtering.read_threshold_analysis import analyze_normalized_reads
from moonstone.validation import check_file_format, has_errors

matplotlib.use('Agg')


class TestExamples(TestCase):

    def test_every_example_is_installed_and_clean(self):
        for name in examples.NAMES:
            with self.subTest(name=name):
                self.assertFalse(has_errors(check_file_format(examples.path(name))))

    def test_unknown_example(self):
        with self.assertRaisesRegex(FileNotFoundError, 'species_counts.tsv'):
            examples.path('nope.csv')

    def test_taxa_counts_fits_the_read_threshold_format(self):
        stats, kept = analyze_normalized_reads(examples.path('taxa_counts.csv'))
        self.assertTrue(0 < stats['species_removed'] < stats['total_taxa'])

    def test_generator_reproduces_the_committed_files(self):
        from moonstone.examples import make_examples
        with tempfile.TemporaryDirectory() as tmp, contextlib.chdir(tmp):
            make_examples.main()
            for name in examples.NAMES:
                with self.subTest(name=name), open(name, 'rb') as new, open(examples.path(name), 'rb') as old:
                    self.assertEqual(new.read(), old.read())
        self.assertTrue(os.path.exists(examples.path('metadata.csv')))
