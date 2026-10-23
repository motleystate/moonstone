import os
import tempfile
from unittest import TestCase, mock

import matplotlib
import numpy as np
import pandas as pd

from moonstone.analysis import classify

matplotlib.use('Agg')


class TestFeatureAnalysis(TestCase):

    def setUp(self):
        rng = np.random.default_rng(0)
        samples = [f's{i}' for i in range(40)]
        self.count_df = pd.DataFrame(rng.poisson(20, (40, 6)).astype(float), index=samples,
                                     columns=[f'taxon{i}' for i in range(6)])
        # taxon0 separates the two groups, so it should be ranked first.
        self.count_df.loc[samples[20:], 'taxon0'] += 30
        self.metadata_df = pd.DataFrame({'group': ['A'] * 20 + ['B'] * 20}, index=samples)
        self.outdir = tempfile.mkdtemp()

    def run_feature_analysis(self):
        ml = classify.ML(self.count_df, self.metadata_df, self.outdir)
        with mock.patch.object(classify.plt, 'plot', wraps=classify.plt.plot) as mocked_plot:
            ml.feature_analysis(variable='group')
        return mocked_plot

    def test_writes_outputs(self):
        self.run_feature_analysis()
        for output_file in ['group_variable_coefs.pdf', 'group_SVM_components.csv', 'group_rfecv.pdf']:
            with self.subTest(output_file=output_file):
                self.assertTrue(os.path.isfile(os.path.join(self.outdir, output_file)))

    def test_components_sorted_by_rank(self):
        self.run_feature_analysis()
        components = pd.read_csv(os.path.join(self.outdir, 'group_SVM_components.csv'), index_col=0)
        self.assertListEqual(list(components.columns), ['Coef'])
        self.assertSetEqual(set(components.index), set(self.count_df.columns))
        self.assertTrue(components['Coef'].is_monotonic_increasing)
        self.assertEqual(components.index[0], 'taxon0')
        self.assertEqual(components['Coef'].iloc[0], 1)

    def test_rfecv_plot_shows_mean_score_per_number_of_features(self):
        mocked_plot = self.run_feature_analysis()
        x, y = mocked_plot.call_args.args
        self.assertListEqual(list(x), list(range(1, self.count_df.shape[1] + 1)))
        self.assertEqual(len(y), self.count_df.shape[1])
        self.assertTrue(np.all((y >= 0) & (y <= 1)))
