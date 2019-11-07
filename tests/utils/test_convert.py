from unittest import TestCase

from moonstone.utils.convert import pandas_to_python_type


class TestPandasToPythonType(TestCase):

    def test_object(self):
        self.assertEqual(pandas_to_python_type('object'), 'str')

    def test_int64(self):
        self.assertEqual(pandas_to_python_type('int64'), 'int')

    def test_float64(self):
        self.assertEqual(pandas_to_python_type('float64'), 'float')

    def test_datetime64(self):
        self.assertEqual(pandas_to_python_type('datetime64'), 'datetime64')
