from unittest import TestCase

import pandas as pd

from moonstone.analysis.diversity.alpha import ShannonIndex


class TestShannonIndex(TestCase):

    def test_compute_shannon_diversity(self):

        tested_object = pd.DataFrame.from_dict(
            {
                'Items1': [4,4,0,0,4],
                'Items2': [1,0,2,0,5],
                'Items3': [0,0,0,1,4],
                'Items4': [0,3,0,0,4]
            },
            orient='index', columns=['sample1', 'sample2', 'sample3', 'sample4', 'sample5'])
        tested_object_instance = ShannonIndex(tested_object)
        # expected_object = pd.Series({
        #     'sample1': 0.721928,
        #     'sample2': 0.985228
        #     'sample3': -0.000000
        #     'sample4': -0.000000
        #     'sample5': 1.992778
        # })
        tested_object_instance.compute_shannon_diversity()
        # TO FINISH
