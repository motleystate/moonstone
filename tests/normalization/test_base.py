from unittest import TestCase

import pandas as pd

from moonstone.normalization.base import BaseNormalization


class TestBaseNormalization(TestCase):

    def test_normalized_df(self):
        df = pd.DataFrame(
            [
                [1, 2, 3],
                [4, 5, 6]
            ],
            columns=['a', 'b', 'c']
        )
        normalization = BaseNormalization(df)
        pd.testing.assert_frame_equal(normalization.normalized_df, df)
