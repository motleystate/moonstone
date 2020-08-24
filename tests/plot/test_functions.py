from unittest import TestCase

import numpy as np
from typing import List

from moonstone.plot.functions import (
    check_list_of,
    check_type,
    check_type_s
)


class TestFunctions(TestCase):

    def test_check_list_of(self):
        self.assertTrue(check_list_of(['A', 'B', 'C'], str))
        self.assertFalse(check_list_of(['A', 12, 'C'], str))
        self.assertTrue(check_list_of([12, 13.5, np.nan], (int, float)))
        self.assertFalse(check_list_of([12, 13.5, 'NaN'], (int, float)))

    def test_check_type(self):
        self.assertTrue(check_type(12, int))
        self.assertTrue(check_type(36.5, float))
        self.assertTrue(check_type('A', str))
        self.assertTrue(check_type(['A', 'B', 'C'], List[str]))

        self.assertEqual(check_type(36.5, str), 'str')
        self.assertEqual(check_type(['A', 12, 'C'], List[str]), 'List[str]')

    def test_check_type_s(self):
        self.assertTrue(check_type_s(['A', 'B', 'C'], List[str]))
        self.assertTrue(check_type_s(36.5, [int, float]))

        self.assertEqual(check_type_s('B', [int, float]), 'int or float')
        print(check_type_s(['A', 12, 'C'], List[str]))  # EUH PROBLEM !!!
