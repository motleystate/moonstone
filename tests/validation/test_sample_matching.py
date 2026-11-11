from unittest import TestCase

import numpy as np

from moonstone.validation.sample_matching import match_samples


def codes(matching):
    return {finding.code: finding.level for finding in matching.findings()}


class TestMatchSamples(TestCase):

    def test_same_samples_same_order(self):
        matching = match_samples(['s1', 's2', 's3'], ['s1', 's2', 's3'])
        self.assertEqual(matching.matched, ['s1', 's2', 's3'])
        self.assertTrue(matching.same_order)
        self.assertEqual(matching.findings(), [])

    def test_different_order_is_matched_by_id(self):
        matching = match_samples(['s1', 's2', 's3'], ['s3', 's1', 's2'])
        self.assertEqual(matching.matched, ['s1', 's2', 's3'])
        self.assertFalse(matching.same_order)
        self.assertEqual(codes(matching), {'different_order': 'info'})

    def test_extra_metadata_samples_are_info(self):
        matching = match_samples(['s1', 's2'], ['s1', 's2', 'x1', 'x2'])
        self.assertEqual(matching.only_in_metadata, ['x1', 'x2'])
        self.assertEqual(codes(matching), {'sample_without_counts': 'info'})

    def test_same_number_of_samples_but_different_ids(self):
        # The old check compared row counts only and accepted this.
        matching = match_samples(['s1', 's2', 's3'], ['s1', 's2', 's4'])
        self.assertEqual(matching.matched, ['s1', 's2'])
        self.assertEqual(matching.only_in_counts, ['s3'])
        self.assertEqual(matching.only_in_metadata, ['s4'])
        self.assertEqual(codes(matching), {'sample_without_metadata': 'warning', 'sample_without_counts': 'info'})

    def test_no_common_sample_is_an_error(self):
        matching = match_samples(['s1', 's2'], ['p1', 'p2'])
        self.assertEqual(matching.matched, [])
        self.assertEqual(codes(matching)['no_common_samples'], 'error')

    def test_duplicated_ids_are_errors(self):
        self.assertEqual(codes(match_samples(['s1', 's1', 's2'], ['s1', 's2']))['duplicated_sample_in_counts'],
                         'error')
        self.assertEqual(codes(match_samples(['s1', 's2'], ['s1', 's2', 's2']))['duplicated_sample_in_metadata'],
                         'error')

    def test_case_and_space_differences_are_reported_not_fixed(self):
        matching = match_samples(['S1', 'S2', 'S3'], ['S1', 's2', 'S3 '])
        self.assertEqual(matching.matched, ['S1'])
        self.assertEqual(matching.near_matches, [('S2', 's2'), ('S3', 'S3 ')])
        self.assertEqual(codes(matching)['near_match'], 'warning')

    def test_ids_compared_as_text(self):
        matching = match_samples(['1', '2'], [1, 2])
        self.assertEqual(matching.matched, ['1', '2'])

    def test_blank_metadata_id_is_ignored_with_a_warning(self):
        matching = match_samples(['s1', 's2'], ['s1', 's2', np.nan])
        self.assertEqual(matching.blank_in_metadata, 1)
        self.assertEqual(matching.only_in_metadata, [])
        self.assertEqual(codes(matching), {'blank_sample_id': 'warning'})

    def test_long_lists_are_shortened(self):
        counts = [f's{i}' for i in range(30)]
        matching = match_samples(counts, ['s0'])
        message = [f.message for f in matching.findings() if f.code == 'sample_without_metadata'][0]
        self.assertIn('(29 in all)', message)
