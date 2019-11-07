import sys
from io import StringIO
from unittest import TestCase

from moonstone.utils.print import print_multi_col


class TestPrintMultiCol(TestCase):

    def setUp(self):
        self.capt_stdout = StringIO()
        sys.stdout = self.capt_stdout

    def tearDown(self):
        sys.stdout = sys.__stdout__

    def test_print_multi_col(self):
        elements = [
            'a:1', 'b:2', 'c:3', 'd:4'
        ]
        line_length = 10
        expected_string = "\n  a:1  b:2\n  c:3  d:4\n"
        print_multi_col(elements, line_length=line_length)
        self.assertEqual(self.capt_stdout.getvalue(), expected_string)
