from unittest import TestCase

import numpy as np

from moonstone.utils.plot import (
    check_list_type, add_x_to_plotting_options
)


class TestUtils(TestCase):

    def test_check_list_type(self):
        self.assertTrue(check_list_type(['A', 'B', 'C'], str))
        self.assertFalse(check_list_type(['A', 12, 'C'], str))
        self.assertTrue(check_list_type([12, 13.5, np.nan], (int, float)))
        self.assertFalse(check_list_type([12, 13.5, 'NaN'], (int, float)))
        self.assertTrue(check_list_type([12, 13, "text"], (int, str)))

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
