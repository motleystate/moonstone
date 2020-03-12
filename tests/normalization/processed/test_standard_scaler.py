import unittest
import pandas as pd


class MyTestCase(unittest.TestCase):

    def test_standard_scaler(self):
        self.assertEqual(True, False)

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


if __name__ == '__main__':
    unittest.main()
