from unittest import TestCase

from moonstone.utils.colors import generate_color_code


class TestGenerateColorCode(TestCase):

    def test_method(self):
        input_string = "testest"
        expected_color = "#c8b99d"
        self.assertEqual(generate_color_code(input_string), expected_color)

    def test_method_short_string(self):
        input_string = "test"
        expected_color = "#c8b979"
        self.assertEqual(generate_color_code(input_string), expected_color)

    def test_method_number_in_str(self):
        input_string = "test123"
        expected_color = "#c8b934"
        self.assertEqual(generate_color_code(input_string), expected_color)
