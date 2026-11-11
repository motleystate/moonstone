import os
import tempfile
from unittest import TestCase, mock

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from sklearn.feature_selection import RFECV
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from moonstone.analysis import classification
from moonstone.analysis.classification import Classification

matplotlib.use('Agg')


def make_data(n_per_group=(20, 20), n_taxa=8, seed=0):
    """Poisson counts where taxon0 is higher in the second group and taxon1 in the third."""
    rng = np.random.default_rng(seed)
    labels = [label for label, n in zip('ABC', n_per_group) for _ in range(n)]
    samples = [f's{i}' for i in range(len(labels))]
    counts = pd.DataFrame(rng.poisson(20, (len(samples), n_taxa)).astype(float), index=samples,
                          columns=[f'taxon{i}' for i in range(n_taxa)])
    labels = pd.Series(labels, index=samples)
    counts.loc[labels == 'B', 'taxon0'] += 30
    counts.loc[labels == 'C', 'taxon1'] += 30
    metadata = pd.DataFrame({'group': labels, 'subject': [f'p{i // 2}' for i in range(len(samples))]}, index=samples)
    return counts, metadata


class TestClassification(TestCase):

    def setUp(self):
        self.counts, self.metadata = make_data()
        self.outdir = tempfile.mkdtemp()

    def make(self, **kwargs):
        kwargs.setdefault('n_estimators', 50)
        return Classification(self.counts, self.metadata, 'group', self.outdir, **kwargs)

    def test_writes_outputs(self):
        for model in ['svm', 'random_forest']:
            self.make(model=model).run(roc=True, feature_selection=True)
            for suffix in ['cv_scores.csv', 'permutation_importance.csv', 'permutation_importance.pdf', 'roc.pdf',
                           'selected_features.csv', 'rfecv.pdf']:
                with self.subTest(model=model, suffix=suffix):
                    self.assertTrue(os.path.isfile(os.path.join(self.outdir, f'group_{model}_{suffix}')))

    def test_optional_outputs_only_when_asked(self):
        self.make().run()
        self.assertTrue(os.path.isfile(os.path.join(self.outdir, 'group_svm_cv_scores.csv')))
        for suffix in ['roc.pdf', 'selected_features.csv', 'rfecv.pdf']:
            self.assertFalse(os.path.isfile(os.path.join(self.outdir, f'group_svm_{suffix}')))

    def test_informative_taxon_found(self):
        for model in ['svm', 'random_forest']:
            with self.subTest(model=model):
                self.make(model=model).run(feature_selection=True)
                importance = pd.read_csv(os.path.join(self.outdir, f'group_{model}_permutation_importance.csv'),
                                         index_col=0)
                self.assertListEqual(list(importance.columns), ['mean_importance', 'std_importance'])
                self.assertEqual(importance.index[0], 'taxon0')
                selected = pd.read_csv(os.path.join(self.outdir, f'group_{model}_selected_features.csv'), index_col=0)
                self.assertEqual(selected.loc['taxon0', 'fraction_selected'], 1.0)
                scores = pd.read_csv(os.path.join(self.outdir, f'group_{model}_cv_scores.csv'), index_col=0)
                self.assertGreater(scores['balanced_accuracy'].mean(), 0.9)

    def test_same_seed_same_result(self):
        first = self.make(model='random_forest').cross_validate()['folds']
        second = self.make(model='random_forest').cross_validate()['folds']
        pd.testing.assert_frame_equal(first, second)

    def test_scaler_fitted_on_training_folds_only(self):
        # The scaler must never see the test samples (no leak from test to training).
        with mock.patch.object(StandardScaler, 'fit', autospec=True, side_effect=StandardScaler.fit) as fit:
            folds = self.make().cross_validate()['folds']
        self.assertListEqual([call.args[1].shape[0] for call in fit.call_args_list], list(folds['n_train']))

    def test_random_forest_not_scaled(self):
        steps = dict(self.make(model='random_forest').build_model().steps)
        self.assertNotIn('scaler', steps)
        self.assertEqual(steps['classifier'].class_weight, 'balanced')

    def test_svm_balanced_without_probability(self):
        svm = self.make().build_model().named_steps['classifier']
        self.assertIsInstance(svm, SVC)
        self.assertEqual(svm.class_weight, 'balanced')
        self.assertIsNot(svm.probability, True)  # 'deprecated' (the default) from scikit-learn 1.9

    def test_feature_selection_nested_in_training_folds(self):
        with mock.patch.object(RFECV, 'fit', autospec=True, side_effect=RFECV.fit) as fit:
            folds = self.make().cross_validate(feature_selection=True)['folds']
        self.assertListEqual([call.args[1].shape[0] for call in fit.call_args_list], list(folds['n_train']))

    def test_permutation_importance_on_test_folds(self):
        with mock.patch.object(classification, 'permutation_importance',
                               wraps=classification.permutation_importance) as permute:
            folds = self.make().cross_validate()['folds']
        self.assertListEqual([call.args[1].shape[0] for call in permute.call_args_list], list(folds['n_test']))
        self.assertEqual(permute.call_args.kwargs['scoring'], 'balanced_accuracy')

    def test_folds_stratified(self):
        x, y, groups = self.make().prepare_data()
        for _, test in self.make().make_folds(y, groups, 5).split(x, y, groups):
            self.assertDictEqual(y.iloc[test].value_counts().to_dict(), {'A': 4, 'B': 4})

    def test_grouped_folds_keep_subjects_together(self):
        x, y, groups = self.make(group_column='subject').prepare_data()
        for train, test in self.make().make_folds(y, groups, 5).split(x, y, groups):
            self.assertFalse(set(groups.iloc[train]) & set(groups.iloc[test]))

    def test_folds_lowered_for_small_group(self):
        self.counts, self.metadata = make_data(n_per_group=(20, 3))
        folds = self.make().cross_validate()['folds']
        self.assertEqual(len(folds), 3)

    def test_group_of_one_sample_rejected(self):
        self.counts, self.metadata = make_data(n_per_group=(20, 1))
        with self.assertRaisesRegex(ValueError, 'at least 2'):
            self.make().cross_validate()

    def test_unbalanced_groups_scored_with_balanced_accuracy(self):
        self.counts, self.metadata = make_data(n_per_group=(30, 10))
        results = self.make().cross_validate()
        self.assertTrue(results['folds']['balanced_accuracy'].between(0, 1).all())
        # Recompute fold 1 by hand from the same split.
        x, y, groups = results['x'], results['y'], None
        train, test = next(self.make().make_folds(y, groups, 5).split(x, y))
        predicted = self.make().build_model().fit(x.iloc[train], y.iloc[train]).predict(x.iloc[test])
        self.assertAlmostEqual(results['folds'].loc[1, 'balanced_accuracy'],
                               balanced_accuracy_score(y.iloc[test], predicted))

    def test_three_groups(self):
        self.counts, self.metadata = make_data(n_per_group=(15, 15, 15))
        for model in ['svm', 'random_forest']:
            with self.subTest(model=model):
                results = self.make(model=model).run(roc=True, feature_selection=True)
                self.assertGreater(results['folds']['auc'].mean(), 0.9)
                importance = pd.read_csv(os.path.join(self.outdir, f'group_{model}_permutation_importance.csv'),
                                         index_col=0)
                self.assertSetEqual(set(importance.index[:2]), {'taxon0', 'taxon1'})

    def test_rbf_kernel(self):
        results = self.make(kernel='rbf').run(roc=True)
        self.assertGreater(results['folds']['auc'].mean(), 0.9)
        with self.assertRaisesRegex(ValueError, 'linear kernel'):
            self.make(kernel='rbf').cross_validate(feature_selection=True)

    def test_explore_kernels(self):
        table = self.make().explore_kernels()
        self.assertSetEqual(set(table.index), {'linear', 'poly', 'rbf', 'sigmoid'})
        self.assertTrue(table['mean_balanced_accuracy'].is_monotonic_decreasing)
        saved = pd.read_csv(os.path.join(self.outdir, 'group_svm_kernel_comparison.csv'), index_col=0)
        pd.testing.assert_frame_equal(saved, table, check_names=False)
        # Linear kernel row matches the main run with the same folds and seed.
        folds = self.make(kernel='linear').cross_validate()['folds']
        self.assertAlmostEqual(table.loc['linear', 'mean_balanced_accuracy'], folds['balanced_accuracy'].mean())
        with self.assertRaisesRegex(ValueError, 'SVM model only'):
            self.make(model='random_forest').explore_kernels()

    def test_auto_kernel_chosen_inside_training_folds(self):
        with mock.patch.object(classification.GridSearchCV, 'fit', autospec=True,
                               side_effect=classification.GridSearchCV.fit) as fit:
            results = self.make(kernel='auto').run(roc=True)
        folds = results['folds']
        self.assertListEqual([call.args[1].shape[0] for call in fit.call_args_list], list(folds['n_train']))
        self.assertTrue(set(folds['kernel']) <= {'linear', 'poly', 'rbf', 'sigmoid'})
        self.assertGreater(folds['balanced_accuracy'].mean(), 0.9)
        with self.assertRaisesRegex(ValueError, 'linear kernel'):
            self.make(kernel='auto').cross_validate(feature_selection=True)

    def test_missing_labels_dropped(self):
        self.metadata.loc['s0', 'group'] = np.nan
        x, y, _ = self.make().prepare_data()
        self.assertEqual(len(y), 39)
        self.assertNotIn('s0', x.index)

    def test_numeric_sample_ids_match_text_ids(self):
        # e.g. sample IDs read as integers from one file and as text from the other
        self.counts.index = range(40)
        self.metadata.index = [str(i) for i in range(40)]
        x, y, _ = self.make().prepare_data()
        self.assertEqual(len(y), 40)
        self.assertEqual(list(x.index), list(y.index))

    def test_one_group_rejected(self):
        self.metadata['group'] = 'A'
        with self.assertRaisesRegex(ValueError, 'needs at least 2 groups'):
            self.make().prepare_data()

    def test_rfecv_plot_labels_balanced_accuracy(self):
        with mock.patch.object(Axes, 'set_ylabel', autospec=True, side_effect=Axes.set_ylabel) as set_ylabel:
            self.make(rfe_step=2).run(feature_selection=True)
        self.assertIn('Inner cross-validated balanced accuracy', [call.args[1] for call in set_ylabel.call_args_list])

    def test_invalid_arguments(self):
        for kwargs in [{'model': 'tree'}, {'kernel': 'cubic'}, {'n_folds': 1}, {'group_column': 'absent'}]:
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                self.make(**kwargs)
