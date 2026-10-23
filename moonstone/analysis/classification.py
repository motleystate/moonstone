"""Supervised classification of samples from their count table.

This module answers one question: *can the microbial profile of a sample predict a metadata
variable* (for example disease status)? It replaces the older ``classify.py`` and
``randomForest.py`` modules, which scaled the data before splitting it, reported plain accuracy
on unbalanced groups, ran unseeded splits, and selected features on the same samples they were
scored on. Each of those mistakes makes a classifier look better than it really is.

How a result is produced, in plain terms:

1. **Cross-validation (CV).** The samples are cut into *k* parts called folds. The model is
   trained on *k - 1* folds and tested on the fold it has never seen; this is repeated so that
   every fold is the test fold once. Only scores on unseen samples are reported, because a model
   always looks good on the samples it was trained on.
2. **Stratified, seeded folds.** Each fold keeps the same proportion of each group as the whole
   data set (stratified), and the random assignment of samples to folds is controlled by a seed,
   so the same command gives the same numbers every time.
3. **Grouped folds (optional).** When one subject gave several samples (time points, body sites),
   those samples are kept together in the same fold. Otherwise the model can recognise the
   subject rather than the condition, and the score is inflated.
4. **Preprocessing inside the training fold.** Any step that learns something from the data
   (here, scaling each taxon to mean 0 and standard deviation 1) is fitted on the training folds
   only, then applied unchanged to the test fold. Scaling all samples first lets information
   about the test samples leak into training. A scikit-learn ``Pipeline`` enforces this order.
5. **Balanced groups.** With 80 healthy and 20 sick samples, a model that always answers
   "healthy" is 80 % accurate and useless. Class weights make errors on the smaller group count
   as much as errors on the larger group during training, and **balanced accuracy** (the
   average of the per-group hit rates) is the main score: its chance level is 1 / number of
   groups, whatever the group sizes.
6. **Feature selection inside the CV loop (nested CV).** Choosing the taxa that best separate
   the groups, and then scoring a model on those taxa with the same samples, is circular. Here
   Recursive Feature Elimination with Cross-Validation (RFECV) runs inside each training part
   only, with its own inner folds, and the outer test fold scores the result.
7. **Permutation importance.** To ask how much the model relies on a taxon, the values of that
   taxon are shuffled among the test samples and the drop in balanced accuracy is measured. A
   large drop means the model depends on that taxon; a drop near zero means it does not.
   This works for every model and, unlike the built-in random forest importance, is measured
   on unseen samples and is not biased towards taxa with many distinct values.
8. **Choosing the SVM kernel.** The kernel decides the shape of the boundary between groups.
   :meth:`Classification.explore_kernels` scores each kernel on the same folds so they can be
   compared, but the best of several scores is optimistic: some of its lead is luck on these
   particular samples. That comparison is therefore labelled exploratory. ``kernel='auto'``
   instead lets an inner cross-validation on each training part pick the kernel, and the outer
   test fold scores that choice, so the score reported for the whole procedure stays fair.
"""
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.model_selection import cross_validate as sklearn_cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

logger = logging.getLogger(__name__)

MODELS = ('svm', 'random_forest')
KERNELS = ('linear', 'poly', 'rbf', 'sigmoid')
MODEL_LABELS = {'svm': 'Support Vector Machine', 'random_forest': 'Random Forest'}


class Classification:
    """Cross-validated classification of one metadata variable from the count table.

    :param count_df: counts with one row per sample and one column per taxon (or gene).
    :param metadata_df: metadata with one row per sample, indexed like ``count_df``.
    :param variable: metadata column holding the groups to predict.
    :param outdir: directory where result files are written.
    :param model: ``'svm'`` (Support Vector Machine, SVM) or ``'random_forest'``.
    :param kernel: SVM kernel. ``'linear'`` separates groups with a flat boundary and allows
        feature selection; ``'poly'`` (polynomial), ``'rbf'`` (radial basis function) and
        ``'sigmoid'`` allow curved boundaries. ``'auto'`` picks one of these four inside each
        training part by nested cross-validation.
    :param n_folds: number of cross-validation folds. Lowered, with a warning, when the
        smallest group has fewer samples (or subjects) than folds.
    :param seed: random seed for the fold assignment, the random forest and the permutations.
    :param group_column: optional metadata column naming the subject of each sample. All
        samples of one subject are kept in the same fold.
    :param n_estimators: number of trees in the random forest.
    :param n_repeats: number of shuffles per taxon for permutation importance.
    :param rfe_step: taxa removed at each round of recursive feature elimination. 1 is the most
        precise; a larger number (or a fraction such as 0.05, meaning 5 % of the taxa) is much
        faster on tables with hundreds of taxa, especially with the random forest.
    """

    def __init__(self, count_df, metadata_df, variable, outdir, model='svm', kernel='linear', n_folds=5,
                 seed=0, group_column=None, n_estimators=500, n_repeats=10, rfe_step=1):
        if model not in MODELS:
            raise ValueError(f'model must be one of {MODELS}, not {model!r}')
        if kernel not in KERNELS + ('auto',):
            raise ValueError(f'kernel must be one of {KERNELS + ("auto",)}, not {kernel!r}')
        if n_folds < 2:
            raise ValueError('n_folds must be at least 2')
        for column in [variable, group_column]:
            if column is not None and column not in metadata_df.columns:
                raise ValueError(f'{column!r} is not a column of the metadata')
        self.count_df = count_df
        self.metadata_df = metadata_df
        self.variable = variable
        self.outdir = outdir
        self.model = model
        self.kernel = kernel
        self.n_folds = n_folds
        self.seed = seed
        self.group_column = group_column
        self.n_estimators = n_estimators
        self.n_repeats = n_repeats
        self.rfe_step = rfe_step

    @property
    def name(self):
        """Short label used in file names and printed output."""
        return f'{self.variable}_{self.model}'

    # ---------------------------------------------------------------------------------------------
    # Data
    # ---------------------------------------------------------------------------------------------
    def prepare_data(self):
        """Return the features ``x``, the groups to predict ``y`` and the subjects ``groups``.

        Samples are matched between the count table and the metadata by name. Samples with no
        value for the variable (or for the subject column) cannot be used and are dropped with
        a warning, rather than being silently treated as a group of their own.
        """
        count_df = self.count_df.copy()
        metadata_df = self.metadata_df.copy()
        if count_df.index.dtype != metadata_df.index.dtype:
            count_df.index = count_df.index.astype(str)
            metadata_df.index = metadata_df.index.astype(str)
        samples = count_df.index.intersection(metadata_df.index, sort=False)
        columns = [self.variable] + ([self.group_column] if self.group_column else [])
        metadata_df = metadata_df.loc[samples, columns].dropna()
        dropped = len(samples) - len(metadata_df)
        if dropped:
            logger.warning(f'{dropped} samples have no value for {columns} and are left out.')
            print(f'Warning: {dropped} samples have no value for {", ".join(columns)} and are left out.')

        x = count_df.loc[metadata_df.index].astype(float)
        y = metadata_df[self.variable]
        groups = metadata_df[self.group_column] if self.group_column else None
        if y.nunique() < 2:
            raise ValueError(f'{self.variable!r} needs at least 2 groups to classify, found {y.nunique()}')
        return x, y, groups

    # ---------------------------------------------------------------------------------------------
    # Model and folds
    # ---------------------------------------------------------------------------------------------
    def build_model(self, kernel=None):
        """Return an untrained model as a scikit-learn ``Pipeline``.

        A pipeline chains the preprocessing and the classifier, so that cross-validation fits
        both on the training folds only (see point 4 of the module description).

        * SVM: taxa are first standardised (mean 0, standard deviation 1), because an SVM
          compares samples by distance and would otherwise be dominated by the most abundant
          taxa. ``class_weight='balanced'`` weights each group by the inverse of its size.
        * Random forest: no scaling, because decision trees only ask "is this value above a
          threshold?", and the answer does not change when a taxon is rescaled.

        Receiver Operating Characteristic (ROC) curves use the SVM decision value (the signed
        distance to the boundary) rather than ``SVC(probability=True)``, which is deprecated in
        scikit-learn and runs an extra internal cross-validation to estimate probabilities.

        :param kernel: SVM kernel to use instead of ``self.kernel`` (``'auto'`` starts from linear;
            the kernel is then set by the inner search in :meth:`cross_validate`).
        """
        if self.model == 'svm':
            kernel = kernel or self.kernel
            return Pipeline([
                ('scaler', StandardScaler()),
                ('classifier', SVC(kernel='linear' if kernel == 'auto' else kernel, class_weight='balanced',
                                   random_state=self.seed)),
            ])
        return Pipeline([
            ('classifier', RandomForestClassifier(n_estimators=self.n_estimators, max_features='sqrt',
                                                  class_weight='balanced', random_state=self.seed)),
        ])

    def importance_getter(self):
        """Where RFECV finds the weight of each taxon in the trained pipeline.

        A linear SVM weights each taxon with a coefficient; a random forest reports how much each
        taxon improves its splits. An SVM with a curved kernel has no per-taxon weight, so it
        cannot be used for recursive feature elimination (nor can ``kernel='auto'``, which may
        pick a curved kernel).
        """
        if self.model == 'svm':
            if self.kernel != 'linear':
                raise ValueError('Feature selection with the SVM needs the linear kernel.')
            return 'named_steps.classifier.coef_'
        return 'named_steps.classifier.feature_importances_'

    def make_folds(self, y, groups, n_folds):
        """Return a stratified (and, with subjects, grouped) fold splitter.

        Every test fold must contain at least one sample of each group, so the number of folds
        cannot exceed the size of the smallest group (counted in subjects when samples are
        grouped). When it would, the number of folds is lowered with a warning; a group with a
        single sample or subject cannot be cross-validated at all.
        """
        if groups is None:
            smallest = int(y.value_counts().min())
            unit = 'samples'
        else:
            smallest = int(pd.DataFrame({'y': y, 'g': groups}).groupby('y')['g'].nunique().min())
            unit = 'subjects'
        if smallest < 2:
            raise ValueError(f'The smallest group of {self.variable!r} has {smallest} {unit}; '
                             'cross-validation needs at least 2.')
        if smallest < n_folds:
            logger.warning(f'Smallest group has {smallest} {unit}: using {smallest} folds instead of {n_folds}.')
            print(f'Warning: the smallest group has {smallest} {unit}, so {smallest} folds are used '
                  f'instead of {n_folds}.')
            n_folds = smallest
        if groups is None:
            return StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=self.seed)
        return StratifiedGroupKFold(n_splits=n_folds, shuffle=True, random_state=self.seed)

    def class_scores(self, fitted, x, classes):
        """Return one score per sample and group: higher means "more likely this group".

        Used for ROC curves and the Area Under the Curve (AUC). The random forest gives the
        fraction of trees voting for each group and the SVM its decision value. For two groups the
        SVM returns a single column (positive means the second group), which is expanded here so
        that both models return one column per group.
        """
        if self.model == 'random_forest':
            return fitted.predict_proba(x)
        scores = fitted.decision_function(x)
        if len(classes) == 2:
            scores = np.column_stack([-scores, scores])
        return scores

    # ---------------------------------------------------------------------------------------------
    # Cross-validation
    # ---------------------------------------------------------------------------------------------
    def cross_validate(self, feature_selection=False):
        """Train and test the model on each fold, and collect the results.

        :param feature_selection: run RFECV inside each training part (nested cross-validation).
        :return: dictionary with the per-fold scores, permutation importances, ROC curves and,
            with feature selection, the taxa kept in each fold.
        """
        x, y, groups = self.prepare_data()
        classes = np.array(sorted(y.unique()))
        outer = self.make_folds(y, groups, self.n_folds)
        model = self.build_model()
        getter = self.importance_getter() if feature_selection else None
        tune_kernel = self.model == 'svm' and self.kernel == 'auto'

        results = {'x': x, 'y': y, 'classes': classes, 'folds': [], 'importance': [], 'roc': [],
                   'support': [], 'rfecv_scores': []}
        for fold, (train, test) in enumerate(outer.split(x, y, groups)):
            x_train, x_test, y_train, y_test = x.iloc[train], x.iloc[test], y.iloc[train], y.iloc[test]
            groups_train = groups.iloc[train] if groups is not None else None
            if tune_kernel:
                # Nested CV for the kernel: the inner folds, made from the training part only,
                # score each kernel; the best one is refitted on the whole training part and the
                # outer test fold, which played no part in the choice, scores it.
                inner = self.make_folds(y_train, groups_train, self.n_folds)
                fitted = GridSearchCV(clone(model), {'classifier__kernel': list(KERNELS)}, cv=inner,
                                      scoring='balanced_accuracy')
                fitted.fit(x_train, y_train, groups=groups_train)
            elif feature_selection:
                # Nested CV: RFECV sees the training part only. It repeatedly drops the taxon with
                # the smallest weight and measures, with inner folds made from the training part,
                # how balanced accuracy changes; it keeps the number of taxa that scored best.
                # The outer test fold, never seen during selection, then scores the chosen model.
                inner = self.make_folds(y_train, groups_train, self.n_folds)
                fitted = RFECV(clone(model), step=self.rfe_step, cv=inner, scoring='balanced_accuracy',
                               importance_getter=getter, min_features_to_select=1)
                fitted.fit(x_train, y_train, groups=groups_train)
                results['support'].append(pd.Series(fitted.support_, index=x.columns))
                results['rfecv_scores'].append(pd.Series(fitted.cv_results_['mean_test_score'],
                                                         index=fitted.cv_results_['n_features']))
            else:
                fitted = clone(model).fit(x_train, y_train)

            predicted = fitted.predict(x_test)
            scores = self.class_scores(fitted, x_test, classes)
            fold_roc = self.roc_for_fold(y_test, scores, classes)
            results['roc'].append(fold_roc)
            fold_auc = np.mean([auc for _, auc in fold_roc.values()]) if fold_roc else np.nan
            results['folds'].append({
                'fold': fold + 1,
                'n_train': len(train),
                'n_test': len(test),
                'balanced_accuracy': balanced_accuracy_score(y_test, predicted),
                'accuracy': accuracy_score(y_test, predicted),
                'auc': fold_auc,
                'n_features': int(fitted.n_features_) if feature_selection else x.shape[1],
            })
            if self.model == 'svm':
                results['folds'][-1]['kernel'] = (fitted.best_params_['classifier__kernel'] if tune_kernel
                                                  else self.kernel)

            # Permutation importance on the test fold: shuffle one taxon at a time among the
            # unseen samples and record the loss of balanced accuracy, averaged over n_repeats
            # shuffles. With feature selection, taxa the model dropped get an importance of 0.
            # Taxa are shuffled in parallel on all processor cores; the seed keeps the result fixed.
            permuted = permutation_importance(fitted, x_test, y_test, scoring='balanced_accuracy',
                                              n_repeats=self.n_repeats, random_state=self.seed, n_jobs=-1)
            results['importance'].append(pd.Series(permuted.importances_mean, index=x.columns))
            logger.info(f'{self.name} fold {fold + 1}: balanced accuracy '
                        f'{results["folds"][-1]["balanced_accuracy"]:.3f}')
        results['folds'] = pd.DataFrame(results['folds']).set_index('fold')
        return results

    def explore_kernels(self):
        """Compare the SVM kernels on the same folds and write ``<variable>_svm_kernel_comparison.csv``.

        Each kernel goes through the same pipeline (scaling fitted on the training folds), the
        same stratified, seeded folds and the same balanced accuracy score as :meth:`run`, so the
        differences between kernels come from the kernels alone.

        These scores are **exploratory**. Picking the kernel with the highest score and quoting
        that score is optimistic, because part of the winner's lead is luck on these samples. To
        report a fair score for a model whose kernel was chosen from the data, use
        ``kernel='auto'``, which makes the choice inside each training part (nested CV).

        :return: one row per kernel, best first: mean and standard deviation over folds of the
            balanced accuracy and accuracy on the test folds.
        """
        if self.model != 'svm':
            raise ValueError('Kernel exploration applies to the SVM model only.')
        x, y, groups = self.prepare_data()
        folds = self.make_folds(y, groups, self.n_folds)
        rows = []
        for kernel in KERNELS:
            scores = sklearn_cross_validate(self.build_model(kernel), x, y, groups=groups, cv=folds,
                                            scoring=['balanced_accuracy', 'accuracy'])
            rows.append({'kernel': kernel,
                         'mean_balanced_accuracy': scores['test_balanced_accuracy'].mean(),
                         'std_balanced_accuracy': scores['test_balanced_accuracy'].std(ddof=1),
                         'mean_accuracy': scores['test_accuracy'].mean(),
                         'std_accuracy': scores['test_accuracy'].std(ddof=1)})
        table = pd.DataFrame(rows).set_index('kernel').sort_values('mean_balanced_accuracy', ascending=False)
        table.to_csv(os.path.join(self.outdir, f'{self.variable}_svm_kernel_comparison.csv'))

        print(f'\nExploratory SVM kernel comparison for {self.variable} ({folds.get_n_splits()}-fold, seed {self.seed}).')
        print('  The best score here is optimistic; use --kernel auto for a fair score of a chosen kernel.')
        for kernel, row in table.iterrows():
            print(f'  {kernel:8s} balanced accuracy {row["mean_balanced_accuracy"]:.1%} '
                  f'(sd {row["std_balanced_accuracy"]:.1%})')
        return table

    @staticmethod
    def roc_for_fold(y_test, scores, classes):
        """Return ``{group: (ROC curve on a common grid, AUC)}`` for one test fold.

        A ROC curve plots, for every possible threshold on the score, the fraction of a group's
        samples correctly detected (true positive rate) against the fraction of other samples
        wrongly flagged (false positive rate). The AUC is 1 for perfect separation and 0.5 for
        guessing. For two groups the second group (in sorted order) is the positive group. For
        three or more groups each group is compared with all the others ("one-vs-rest") and the
        AUC reported is the average over groups.
        """
        grid = np.linspace(0, 1, 101)
        curves = {}
        targets = classes[1:] if len(classes) == 2 else classes
        for label in targets:
            column = int(np.flatnonzero(classes == label)[0])
            is_label = (y_test.to_numpy() == label)
            if is_label.all() or not is_label.any():
                continue  # a grouped fold can miss a group; its curve is undefined
            fpr, tpr, _ = roc_curve(is_label, scores[:, column])
            curve = np.interp(grid, fpr, tpr)
            curve[0] = 0.0
            curves[label] = (curve, roc_auc_score(is_label, scores[:, column]))
        return curves

    # ---------------------------------------------------------------------------------------------
    # Outputs
    # ---------------------------------------------------------------------------------------------
    def run(self, roc=False, feature_selection=False):
        """Run the cross-validation, write the result files and print a summary.

        Files written to ``outdir`` (``NAME`` is ``<variable>_<model>``):

        * ``NAME_cv_scores.csv``: balanced accuracy, accuracy and AUC of each test fold.
        * ``NAME_permutation_importance.csv`` and ``.pdf``: mean and standard deviation over folds
          of the drop in balanced accuracy when each taxon is shuffled.
        * ``NAME_roc.pdf`` (with ``roc``): ROC curves on the test folds.
        * ``NAME_selected_features.csv`` and ``NAME_rfecv.pdf`` (with ``feature_selection``): how
          often each taxon was kept across folds, and the inner score against number of taxa.

        :return: the dictionary returned by :meth:`cross_validate`.
        """
        logger.info(f'Starting {self.name} classification.')
        results = self.cross_validate(feature_selection=feature_selection)
        results['folds'].to_csv(self.output_path('cv_scores.csv'))
        importance = self.write_importance(results)
        if roc:
            self.plot_roc(results)
        if feature_selection:
            self.write_feature_selection(results)
        self.print_summary(results, importance, feature_selection)
        logger.info(f'{self.name} classification completed.')
        return results

    def output_path(self, suffix):
        return os.path.join(self.outdir, f'{self.name}_{suffix}')

    def write_importance(self, results, top=30):
        importance = pd.concat(results['importance'], axis=1)
        table = pd.DataFrame({'mean_importance': importance.mean(axis=1), 'std_importance': importance.std(axis=1)})
        table = table.sort_values('mean_importance', ascending=False)
        table.index.name = 'feature'
        table.to_csv(self.output_path('permutation_importance.csv'))

        shown = table.head(top).iloc[::-1]
        fig, ax = plt.subplots(figsize=(6, max(3, 0.25 * len(shown) + 1)))
        ax.barh(shown.index.astype(str), shown['mean_importance'], xerr=shown['std_importance'], color='tab:blue')
        ax.axvline(0, color='black', lw=0.8)
        ax.set_xlabel('Drop in balanced accuracy when shuffled (mean ± sd over folds)')
        ax.set_title(f'{self.variable}: {MODEL_LABELS[self.model]} permutation importance', fontsize=10)
        fig.tight_layout()
        fig.savefig(self.output_path('permutation_importance.pdf'), format='pdf', dpi=150)
        plt.close(fig)
        return table

    def plot_roc(self, results):
        grid = np.linspace(0, 1, 101)
        classes = results['classes']
        fig, ax = plt.subplots(figsize=(6, 6))
        if len(classes) == 2:
            curves, aucs = zip(*[fold[classes[1]] for fold in results['roc'] if classes[1] in fold])
            for i, (curve, auc) in enumerate(zip(curves, aucs)):
                ax.plot(grid, curve, lw=1, alpha=0.3, label=f'Fold {i + 1} (AUC = {auc:.2f})')
            self._plot_mean_curve(ax, grid, curves, aucs, 'Mean', 'tab:blue')
            title = f'{self.variable} ({classes[1]} as positive)'
        else:
            for label, colour in zip(classes, plt.cm.tab10.colors * (len(classes) // 10 + 1)):
                curves, aucs = zip(*[fold[label] for fold in results['roc'] if label in fold])
                self._plot_mean_curve(ax, grid, curves, aucs, f'{label} vs rest', colour)
            title = f'{self.variable} (one-vs-rest)'
        ax.plot([0, 1], [0, 1], linestyle='--', color='grey', label='Chance')
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel('False positive rate')
        ax.set_ylabel('True positive rate')
        ax.set_title(f'{MODEL_LABELS[self.model]} ROC on test folds\n{title}', fontsize=11)
        ax.legend(loc='lower right', fontsize=8)
        fig.tight_layout()
        fig.savefig(self.output_path('roc.pdf'), format='pdf', dpi=150)
        plt.close(fig)

    @staticmethod
    def _plot_mean_curve(ax, grid, curves, aucs, label, colour):
        """Plot the mean of the fold curves with a band of one standard deviation."""
        mean_curve = np.mean(curves, axis=0)
        mean_curve[-1] = 1.0
        spread = np.std(curves, axis=0)
        ax.plot(grid, mean_curve, color=colour, lw=2,
                label=f'{label} (AUC = {np.nanmean(aucs):.2f} ± {np.nanstd(aucs):.2f})')
        ax.fill_between(grid, np.maximum(mean_curve - spread, 0), np.minimum(mean_curve + spread, 1),
                        color=colour, alpha=0.15)

    def write_feature_selection(self, results):
        support = pd.concat(results['support'], axis=1)
        table = pd.DataFrame({'times_selected': support.sum(axis=1).astype(int),
                              'fraction_selected': support.mean(axis=1)})
        table = table.sort_values('times_selected', ascending=False)
        table.index.name = 'feature'
        table.to_csv(self.output_path('selected_features.csv'))

        # Each outer fold has its own inner curve: score against number of taxa kept. Their mean
        # shows how many taxa are needed; it is computed on training parts only, so the honest
        # score of the selected models is the one in cv_scores.csv, not the top of this curve.
        curves = pd.concat(results['rfecv_scores'], axis=1).sort_index()
        mean, spread = curves.mean(axis=1), curves.std(axis=1).fillna(0)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(curves.index, mean, color='tab:blue', marker='.')
        ax.fill_between(curves.index, mean - spread, mean + spread, color='tab:blue', alpha=0.15)
        ax.set_xlabel('Number of features kept')
        ax.set_ylabel('Inner cross-validated balanced accuracy')
        ax.set_title(f'{self.variable}: {MODEL_LABELS[self.model]} recursive feature elimination', fontsize=10)
        fig.tight_layout()
        fig.savefig(self.output_path('rfecv.pdf'), format='pdf', dpi=150)
        plt.close(fig)

    def print_summary(self, results, importance, feature_selection):
        folds, y = results['folds'], results['y']
        model = MODEL_LABELS[self.model] + (f' ({self.kernel} kernel)' if self.model == 'svm' else '')
        grouping = f', grouped by {self.group_column}' if self.group_column else ''
        print(f'\n{model} classification of {self.variable}: {len(folds)}-fold stratified '
              f'cross-validation{grouping}, seed {self.seed}')
        counts = y.value_counts().sort_index()
        print(f'  {len(y)} samples: ' + ', '.join(f'{label} {n} ({n / len(y):.1%})' for label, n in counts.items()))
        print(f'  Balanced accuracy on test folds: {folds["balanced_accuracy"].mean():.1%} '
              f'(sd {folds["balanced_accuracy"].std():.1%}); chance level {1 / len(counts):.1%}')
        print(f'  Accuracy on test folds: {folds["accuracy"].mean():.1%} (sd {folds["accuracy"].std():.1%})')
        print(f'  Area under the ROC curve: {folds["auc"].mean():.2f} (sd {folds["auc"].std():.2f})')
        if self.kernel == 'auto' and self.model == 'svm':
            print('  Kernel chosen per fold (inner cross-validation): ' + ', '.join(folds['kernel']))
        if feature_selection:
            print('  Features kept per fold: ' + ', '.join(str(n) for n in folds['n_features']))
        top = importance.head(5)
        print('  Top permutation importance: ' +
              ', '.join(f'{feature} ({value:.3f})' for feature, value in top['mean_importance'].items()))
        print(f'  Results written to {self.output_path("*")}')
