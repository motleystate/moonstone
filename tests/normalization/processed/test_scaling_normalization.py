import unittest
import pandas as pd

from moonstone.normalization.processed.scaling_normalization import StandardScaler


class MyTestCase(unittest.TestCase):

    def test_standard_scaler(self):
        df = pd.DataFrame(
            [
                [1., -1., 2.],
                [2., 0., 0.],
                [0., 1., -1]
            ],
            columns=['A', 'B', 'C']
        )
        df_expected = pd.DataFrame(
            [
                [0., -1.224745, 1.336306],
                [1.224745, 0., -0.267261],
                [-1.224745, 1.224745, -1.069045],
            ],
            columns=['A', 'B', 'C']
        )
        scaler = StandardScaler(df).scale_x
        pd.testing.assert_frame_equal(scaler, df_expected)


if __name__ == '__main__':
    unittest.main()
