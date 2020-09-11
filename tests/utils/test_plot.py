from unittest import TestCase

import numpy as np

from moonstone.utils.plot import (
    check_list_type, add_x_to_plotting_options, add_default_titles_to_plotting_options
)


class TestUtils(TestCase):

    def test_check_list_type(self):
        self.assertTrue(check_list_type(['A', 'B', 'C'], str))
        self.assertFalse(check_list_type(['A', 12, 'C'], str))
        self.assertTrue(check_list_type([12, 13.5, np.nan], (int, float)))
        self.assertFalse(check_list_type([12, 13.5, 'NaN'], (int, float)))
        self.assertTrue(check_list_type([12, 13, "text"], (int, str)))

    def test_add_x_to_plotting_options(self):
        # test without option_cat dictionary
        plotting_options = {}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'traces', 'marker_color', '#ee048e'),
                             {'traces': {'marker_color': '#ee048e'}})

        # test with option_cat dictionary but without x
        plotting_options = {'traces': {'marker_opacity': 0.5}}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'traces', 'marker_color', '#ee048e'),
                             {'traces': {'marker_opacity': 0.5, 'marker_color': '#ee048e'}})

        # test with x already in plotting options[option_cat]
        plotting_options = {'traces': {'marker_color': 'red'}}
        self.assertDictEqual(add_x_to_plotting_options(plotting_options, 'traces', 'marker_color', '#ee048e'),
                             {'traces': {'marker_color': 'red'}})

    def test_add_default_titles_to_plotting_options(self):
        title = "graph about something"
        xlabel = "something"
        ylabel = "samples"

        # test with empty dictionary
        plotting_options = {}
        expected_plotting_options = {'layout': {'title_text': 'graph about something', 'title_x': 0.5},
                                     'xaxes': {'title_text': 'something'},
                                     'yaxes': {'title_text': 'samples'}}
        self.assertDictEqual(add_default_titles_to_plotting_options(plotting_options, title, xlabel, ylabel),
                             expected_plotting_options)
