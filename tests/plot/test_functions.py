from unittest import TestCase

import numpy as np
from typing import List
import warnings

from moonstone.plot.functions import (
    check_list_type,
    check_type,
    check_type_s,
    check_types_in_plotting_options,
    add_x_to_plotting_options
)


class TestFunctions(TestCase):

    def test_check_list_type(self):
        self.assertTrue(check_list_type(['A', 'B', 'C'], str))
        self.assertFalse(check_list_type(['A', 12, 'C'], str))
        self.assertTrue(check_list_type([12, 13.5, np.nan], (int, float)))
        self.assertFalse(check_list_type([12, 13.5, 'NaN'], (int, float)))
        self.assertTrue(check_list_type([12, 13, "text"], (int, str)))

    def test_check_type(self):
        self.assertTrue(check_type(12, int))
        self.assertTrue(check_type(36.5, float))
        self.assertTrue(check_type('A', str))
        self.assertTrue(check_type(['A', 'B', 'C'], list))

        self.assertEqual(check_type(36.5, str), 'str')
        self.assertEqual(check_type(['A', 12, 'C'], List[str]), 'List[str]')

    def test_check_type_s(self):
        self.assertTrue(check_type_s(['A', 'B', 'C'], List[str]))
        self.assertTrue(check_type_s(36.5, [int, float]))

        self.assertEqual(check_type_s('B', [int, float]), 'int or float')
        self.assertEqual(check_type_s(['A', 12, 'C'], List[str]), 'List[str]')  # EUH PROBLEM !!!

    def test_check_types_in_plotting_options(self):
        good_plotting_options = {'log': True, 'colorbar': '#ee048e'}
        self.assertDictEqual(check_types_in_plotting_options(good_plotting_options), good_plotting_options)

        dirty_plotting_options = {'log': 'bad'}
        with warnings.catch_warnings(record=True) as w:
            self.assertDictEqual(check_types_in_plotting_options(dirty_plotting_options), {})
            assert len(w) == 1
            assert str(w[0].message) == 'Warning : log value in plotting_options must \
be a bool, str given. Value overridden'

        dirty_plotting_options = {'colorbar': ['red', 63]}
        with warnings.catch_warnings(record=True) as w:
            self.assertDictEqual(check_types_in_plotting_options(dirty_plotting_options), {})
            assert len(w) == 1
            assert str(w[0].message) == 'Warning : colorbar value in plotting_options must \
be a str or List[str]. Please check the type of elements in the list given. Value overridden'

        dirty_plotting_options = {'log': 'bad', 'colorbar': ['red', 63], 'tickangle': 50}
        with warnings.catch_warnings(record=True) as w:
            self.assertDictEqual(check_types_in_plotting_options(dirty_plotting_options), {'tickangle': 50})
            assert len(w) == 2
            assert str(w[0].message) == 'Warning : log value in plotting_options must \
be a bool, str given. Value overridden'
            assert str(w[1].message) == 'Warning : colorbar value in plotting_options must \
be a str or List[str]. Please check the type of elements in the list given. Value overridden'

    def test_add_x_to_plotting_options(self):
        # test with empty dictionnary
        plotting_options = {}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'colorbar', '#ee048e'),
                             {'colorbar': '#ee048e'})

        # test with not-empty dictionnary
        plotting_options = {'log': True}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'colorbar', '#ee048e'),
                             {'log': True, 'colorbar': '#ee048e'})

        # test with value already in plotting options
        plotting_options = {'colorbar': 'red'}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'colorbar', '#ee048e'),
                             {'colorbar': 'red'})
