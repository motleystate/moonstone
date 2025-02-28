from unittest import TestCase
import pytest

import numpy as np
import pandas as pd

from moonstone.utils.pandas.remodelling import StructureRemodelling


class TestStructureRemodelling(TestCase):
    """
    asymmetric      -- _make_series_symmetric -->        symmetric
    series          <- _make_series_asymmetric --        series
      | Ʌ                                                  | Ʌ
    * | | **             * _transform_series_to_dataframe  | | ** _transform_dataframe_to_series
      V |                                                  V |
    asymetric       -- _make_dataframe_symmetric -->    symmetric
    dataframe       <- _make_dataframe_asymmetric --    dataframe
    """

    @pytest.fixture(autouse=True)
    def _pass_fixtures(self, capsys):
        self.capsys = capsys

    def test_check_symmetric_series_asym1(self):
        # asymmetric series case 1: unilateral (ex: ("A", "B") in the index but not ("B", "A"))
        tested_object = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "C"): 4,
            ("B", "D"): 5,
            ("C", "D"): 6
        })
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as an asymmetric series\n')
        # can seems redundant but check to be sure:
        self.assertEqual(tested_instance.check_symmetric_series(tested_object), False)

    def test_check_symmetric_series_asym2(self):
        # asymmetric series case 2: different value for both sides (ex: ("A", "B") = 1 =/= ("B", "A") = 4)
        tested_object = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "A"): 4,
            ("B", "C"): 5,
            ("B", "D"): 6,
            ("C", "A"): 7,
            ("C", "B"): 8,
            ("C", "D"): 9,
            ("D", "A"): 10,
            ("D", "B"): 11,
            ("D", "C"): 12,
        })
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as an asymmetric series\n')
        self.assertEqual(tested_instance.check_symmetric_series(tested_object), False)

    def test_check_symmetric_series_sym(self):
        # symmetric series
        tested_object = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "A"): 1,
            ("B", "C"): 4,
            ("B", "D"): 5,
            ("C", "A"): 2,
            ("C", "B"): 4,
            ("C", "D"): 6,
            ("D", "A"): 3,
            ("D", "B"): 5,
            ("D", "C"): 6,
        })
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as a symmetric series\n')
        self.assertEqual(tested_instance.check_symmetric_series(tested_object), True)

    def test_check_symmetric_dataframe_asym1(self):
        # asymmetric dataframe:
        tested_object = pd.DataFrame.from_dict({
            'A': [1, 2, 3],
            'B': [np.nan, 4, 5],
            'C': [np.nan, np.nan, 6],
        }, orient="index", columns=["B", "C", "D"])
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as an asymmetric dataframe\n')
        self.assertEqual(tested_instance.check_symmetric_dataframe(tested_object), False)

    def test_check_symmetric_dataframe_asym2(self):
        # asymmetric dataframe
        tested_object = pd.DataFrame.from_dict({
            'A': [0, 1, 2, 3],
            'B': [4, 5, 6, 7],
            'C': [8, 9, 10, 11],
            'D': [12, 13, 14, 15]
        }, orient="index", columns=["A", "B", "C", "D"])
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as an asymmetric dataframe\n')
        self.assertEqual(tested_instance.check_symmetric_dataframe(tested_object), False)

    def test_check_symmetric_dataframe_sym(self):
        # symmetric dataframe
        tested_object = pd.DataFrame.from_dict({
            'A': [np.nan, 1.0, 2.0, 3.0],
            'B': [1.0, np.nan, 4.0, 5.0],
            'C': [2.0, 4.0, np.nan, 6.0],
            'D': [3.0, 5.0, 6.0, np.nan]
        }, orient="index", columns=["A", "B", "C", "D"])
        tested_instance = StructureRemodelling(tested_object)
        captured = self.capsys.readouterr()
        self.assertEqual(captured.out, 'object stored as a symmetric dataframe\n')
        self.assertEqual(tested_instance.check_symmetric_dataframe(tested_object), True)

    def test_from_unsymmetric_series(self):
        tested_asym_series = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "C"): 4,
            ("B", "D"): 5,
            ("C", "D"): 6
        })

        expected_sym_series = pd.Series({
            ("A", "A"): np.nan,
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "A"): 1,
            ("B", "B"): np.nan,
            ("B", "C"): 4,
            ("B", "D"): 5,
            ("C", "A"): 2,
            ("C", "B"): 4,
            ("C", "C"): np.nan,
            ("C", "D"): 6,
            ("D", "A"): 3,
            ("D", "B"): 5,
            ("D", "C"): 6,
            ("D", "D"): np.nan,
        })

        expected_asym_dataframe = pd.DataFrame.from_dict({
            'A': [1.0, 2.0, 3.0],
            'B': [np.nan, 4.0, 5.0],
            'C': [np.nan, np.nan, 6.0],
        }, orient="index", columns=["B", "C", "D"])

        expected_sym_dataframe = pd.DataFrame.from_dict({
            'A': [np.nan, 1.0, 2.0, 3.0],
            'B': [1.0, np.nan, 4.0, 5.0],
            'C': [2.0, 4.0, np.nan, 6.0],
            'D': [3.0, 5.0, 6.0, np.nan]
        }, orient="index", columns=["A", "B", "C", "D"])

        tested_instance = StructureRemodelling(tested_asym_series, sym=False)

        # _make_series_symmetric
        pd.testing.assert_series_equal(tested_instance.symmetric_series.sort_index(), expected_sym_series)

        # _transform_series_to_dataframe
        pd.testing.assert_frame_equal(tested_instance.asymmetric_dataframe, expected_asym_dataframe)

        # _make_dataframe_symmetric
        pd.testing.assert_frame_equal(tested_instance.symmetric_dataframe, expected_sym_dataframe)

    def test_from_symmetric_dataframe(self):
        tested_sym_dataframe = pd.DataFrame.from_dict({
            'A': [np.nan, -1, 0, np.nan],
            'B': [-1, np.nan, np.inf, (1/3)],
            'C': [0, np.inf, np.nan, 1.2*10**-6],
            'D': [np.nan, (1/3), 1.2*10**-6, np.nan]
        }, orient="index", columns=["A", "B", "C", "D"])

        expected_sym_series = pd.Series({
            ("A", "A"): np.nan,
            ("A", "B"): -1,
            ("A", "C"): 0,
            ("A", "D"): np.nan,
            ("B", "A"): -1,
            ("B", "B"): np.nan,
            ("B", "C"): np.inf,
            ("B", "D"): (1/3),
            ("C", "A"): 0,
            ("C", "B"): np.inf,
            ("C", "C"): np.nan,
            ("C", "D"): 1.2*10**-6,
            ("D", "A"): np.nan,
            ("D", "B"): (1/3),
            ("D", "C"): 1.2*10**-6,
            ("D", "D"): np.nan,
        })

        expected_asym_series = pd.Series({
            ("A", "B"): -1,
            ("A", "C"): 0,
            ("B", "C"): np.inf,
            ("B", "D"): (1/3),
            ("C", "D"): 1.2*10**-6,
        })

        expected_asym_dataframe = pd.DataFrame.from_dict({
            'A': [-1, 0, np.nan],
            'B': [np.nan, np.inf, (1/3)],
            'C': [np.nan, np.nan, 1.2*10**-6],
        }, orient="index", columns=["B", "C", "D"])

        tested_instance = StructureRemodelling(tested_sym_dataframe, sym=True)

        # _transform_dataframe_to_series
        pd.testing.assert_series_equal(tested_instance.symmetric_series, expected_sym_series)

        # _make_series_asymmetric
        pd.testing.assert_series_equal(tested_instance.asymmetric_series, expected_asym_series)

        # _make_dataframe_asymmetric
        pd.testing.assert_frame_equal(tested_instance.asymmetric_dataframe, expected_asym_dataframe)

    def test_create_symmetric_series(self):
        # check that it gives the same series from the 2 paths
        tested_asym_dataframe = pd.DataFrame.from_dict({
            'A': [1, 2, np.nan],
            'B': [np.nan, 3, 4],
            'C': [np.nan, np.nan, 5],
        }, orient="index", columns=["B", "C", "D"])

        tested_instance = StructureRemodelling(tested_asym_dataframe, sym=False)
        obj_path1 = tested_instance.symmetric_series  # asym df -> sym df -> sym ser

        tested_instance = StructureRemodelling(tested_asym_dataframe, sym=False)
        tested_instance.asymmetric_series
        obj_path2 = tested_instance.symmetric_series  # asym df -> asym ser -> sym ser

        pd.testing.assert_series_equal(obj_path1, obj_path2)

    def test_create_asymmetric_series(self):
        # check that it gives the same series from the 2 paths
        tested_sym_dataframe = pd.DataFrame.from_dict({
            'A': [0, 1.0, 2.0, 3.0],
            'B': [1.0, 0, 4.0, 5.0],
            'C': [2.0, 4.0, 0, np.nan],
            'D': [3.0, 5.0, np.nan, 0]
        }, orient="index", columns=["A", "B", "C", "D"])

        tested_instance = StructureRemodelling(tested_sym_dataframe, sym=True)
        obj_path1 = tested_instance.asymmetric_series  # sym df -> asym df -> asym ser

        tested_instance = StructureRemodelling(tested_sym_dataframe, sym=True)
        tested_instance.symmetric_series
        obj_path2 = tested_instance.asymmetric_series  # sym df -> sym ser -> asym ser

        pd.testing.assert_series_equal(obj_path1, obj_path2)

    def test_create_symmetric_dataframe(self):
        # check that it gives the same series from the 2 paths
        tested_asym_series = pd.Series({
            ("A", "B"): -1,
            ("A", "C"): 0,
            ("B", "B"): 0,
            ("B", "C"): np.inf,
            ("B", "D"): (1/3),
            ("C", "D"): 1.2*10**-6,
        })

        tested_instance = StructureRemodelling(tested_asym_series, sym=False)
        obj_path1 = tested_instance.symmetric_dataframe  # asym ser -> sym ser -> sym df

        tested_instance = StructureRemodelling(tested_asym_series, sym=False)
        tested_instance.asymmetric_dataframe
        obj_path2 = tested_instance.symmetric_dataframe  # asym ser -> asym df -> sym df

        pd.testing.assert_frame_equal(obj_path1, obj_path2)

    def test_create_asymmetric_dataframe(self):
        # check that it gives the same series from the 2 paths
        tested_sym_series = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): np.nan,
            ("B", "A"): 1,
            ("B", "C"): 3,
            ("B", "D"): 4,
            ("C", "A"): 2,
            ("C", "B"): 3,
            ("C", "D"): 5,
            ("D", "A"): np.nan,
            ("D", "B"): 4,
            ("D", "C"): 5,
        })

        tested_instance = StructureRemodelling(tested_sym_series, sym=True)
        obj_path1 = tested_instance.asymmetric_dataframe  # sym ser -> asym ser -> asym df

        tested_instance = StructureRemodelling(tested_sym_series, sym=True)
        tested_instance.symmetric_dataframe
        obj_path2 = tested_instance.asymmetric_dataframe  # sym ser -> sym df -> asym df

        pd.testing.assert_frame_equal(obj_path1, obj_path2)
    
    def test_get_from_arguments(self):
        tested_object = pd.Series({
            ("A", "B"): 1,
            ("A", "C"): 2,
            ("A", "D"): 3,
            ("B", "C"): 4,
            ("B", "D"): 5,
            ("C", "D"): 6
        })
        tested_instance = StructureRemodelling(tested_object)

        expected_object = pd.DataFrame.from_dict({
            'A': [np.nan, 1.0, 2.0, 3.0],
            'B': [1.0, np.nan, 4.0, 5.0],
            'C': [2.0, 4.0, np.nan, 6.0],
            'D': [3.0, 5.0, 6.0, np.nan]
        }, orient="index", columns=["A", "B", "C", "D"])

        pd.testing.assert_frame_equal(
            tested_instance.get_from_arguments(structure="dataframe", sym=True), 
            expected_object
        )