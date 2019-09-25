from unittest import TestCase

from metautils.utils.convert import pandas_to_python_type


class TestPandasToPythonType(TestCase):

    def test_object(self):
        self.assertEqual(pandas_to_python_type('object'), 'str (object)')

    def test_int64(self):
        self.assertEqual(pandas_to_python_type('int64'), 'int (int64)')

    def test_float64(self):
        self.assertEqual(pandas_to_python_type('float64'), 'float (float64)')

    def test_datetime64(self):
        self.assertEqual(pandas_to_python_type('datetime64'), 'datetime64')
