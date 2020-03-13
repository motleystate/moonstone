import unittest
import numpy as np

from moonstone.normalization.processed.scaling_normalization import StandardScaler


class TestStandardScaler(unittest.TestCase):

    def test_standard_scaler(self):
        raw = np.array(
            [
                [1., -1., 2.],
                [2., 0., 0.],
                [0., 1., -1]
            ]
        )
        expected = np.array(
            [
                [0., -1.224745, 1.336306],
                [1.224745, 0., -0.267261],
                [-1.224745, 1.224745, -1.069045],
            ]
        )
        scaled_x = StandardScaler(raw).scaled_x
        np.array_equal(scaled_x, expected)
