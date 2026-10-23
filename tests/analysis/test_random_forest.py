import os
import tempfile
from unittest import TestCase, mock

import matplotlib
import numpy as np
import pandas as pd
from sklearn.ensemble import AdaBoostClassifier, GradientBoostingClassifier, RandomForestClassifier

from moonstone.analysis import randomForest

matplotlib.use('Agg')


class TestRandomForest(TestCase):

    def setUp(self):
        rng = np.random.default_rng(0)
        samples = [f's{i}' for i in range(40)]
        self.count_df = pd.DataFrame(rng.poisson(20, (40, 12)).astype(float), index=samples,
                                     columns=[f'taxon{i}' for i in range(12)])
        # A taxon absent from every sample gets zero importance.
        self.count_df['absent_taxon'] = 0.0
        self.metadata_df = pd.DataFrame({'group': ['A'] * 20 + ['B'] * 20}, index=samples)
        self.outdir = tempfile.mkdtemp()

    def make_forest(self):
        return randomForest.RandomForest(self.count_df, self.metadata_df, self.outdir, variable='group')

    def test_get_matrix(self):
        df = self.make_forest().get_matrix()
        self.assertListEqual(list(df.columns), ['group'] + list(self.count_df.columns))
        pd.testing.assert_series_equal(df['group'], self.metadata_df['group'])

    def test_forest_writes_outputs(self):
        self.make_forest().forest('rf_Allfeatures.csv')
        features = pd.read_csv(os.path.join(self.outdir, 'rf_Allfeatures.csv'), index_col=0)
        self.assertListEqual(list(features.columns), ['features', 'importance', 'std'])
        self.assertSetEqual(set(features['features']), set(self.count_df.columns))
        self.assertTrue(os.path.isfile(os.path.join(self.outdir, 'group_rfFeatures.pdf')))

    def test_forest_uses_sqrt_max_features(self):
        with mock.patch.object(randomForest, 'RandomForestClassifier', wraps=RandomForestClassifier) as rf_class:
            self.make_forest().forest('rf.csv')
        self.assertEqual(rf_class.call_args.kwargs['max_features'], 'sqrt')

    def test_gradient_boost_split_uses_gradient_boosting(self):
        # One fit on the train/test split plus one per fold (10 folds) for each boosting model.
        with mock.patch.object(AdaBoostClassifier, 'fit', autospec=True,
                               side_effect=AdaBoostClassifier.fit) as ab_fit, \
                mock.patch.object(GradientBoostingClassifier, 'fit', autospec=True,
                                  side_effect=GradientBoostingClassifier.fit) as gb_fit:
            self.make_forest().forest('rf.csv')
        self.assertEqual(ab_fit.call_count, 11)
        self.assertEqual(gb_fit.call_count, 11)
