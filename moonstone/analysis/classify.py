import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interp
from sklearn import preprocessing, svm
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.feature_selection import RFECV

from moonstone.analysis import stats


class SVM(object):
    def __init__(self, countfile, metadata):
        self.countfile = countfile
        self.metadata = metadata

    def merge(self, variable):
        dc = self.countfile
        dm = self.metadata

        """The naming/numbering of samples is going to be variable.
        Merging appear to be the best option to preserve sample designations
        and combine count data with any metadata.
        Pandas data frames requires indexes to be of the same type in order to correctly merge.

        For the moment is appear that metadata samples, if they are numbered, result in a int64 index type.
        Count data yields an 'object' type, even if samples are numbered.

        Thus we will check for mismatched index type and try to correct"""
        if not isinstance(dc.index, type(dm.index)):
            dc.set_index(np.int64(np.array(dc.index)), inplace=True)

        df = pd.merge(dm[variable], dc, left_index=True, right_index=True)
        return df  # returned to analyze function as 'df_final'

    def analyze(self, variable=""):
        if variable == 'all':
            pass  # Later to include an iterator over al clinical variables.

        df_final = SVM.merge(self, variable)

        # Setup the features and labels
        x = np.array(df_final.drop([variable], axis=1).astype(float))
        # RBF kernel assumes features centered around zero, with equal variance.
        # Metagenomic data is sparse. maxabs_scale() handles both of these conditions.
        x_maxabs_scaled = preprocessing.maxabs_scale(x)
        print("Counts standardized by Maximum Absolute Value. Check Mean & Variance near Zero:")
        stats.normalized_stats(x_maxabs_scaled)
        y = np.array(df_final[variable])

        # My variable counter
        stats.count_items(y)

        # Support Vector Machine
        print('\nRunning Support Vector Machine classifier. This will take a few moments...')
        for kernel in ["linear", "poly", "rbf", "sigmoid"]:
            clf = svm.SVC(kernel=kernel, tol=0.00001, gamma='scale')
            # Get accuracy
            score = cross_val_score(clf, x_maxabs_scaled, y, cv=10)
            print("\tAccuracy for %s kernel: %.2f%s (+/- %.2f%s)" %
                  (kernel, score.mean()*100, "%", score.std()*200, "%"))

    def roc_analysis(self, variable=""):
        df_final = SVM.merge(self, variable)

        x = np.array(df_final.drop([variable], axis=1).astype(float))
        x = preprocessing.maxabs_scale(x)
        y = np.array(df_final[variable])

        # Here we make sure 2, and only 2 categories are present.
        if len(set(y)) != 2:
            raise SystemExit("ROC Analysis is only valid for binary [2] variables")
        # My variable counter
        stats.count_items(y)

        # y = label_binarize(y, classes=[1, 2]).ravel()
        # Here is where the sample designations are interpreted. There is now no need to explicitly name
        # each variable category. The LabelEncoder will simply assign numbers to each new category.
        labeler = preprocessing.LabelEncoder()
        labeler.fit(y)
        y = labeler.transform(y)

        clf = svm.SVC(kernel='linear', gamma='auto', probability=True, tol=0.00001)
        cv = StratifiedKFold(shuffle=True, n_splits=3)
        n_samples, n_features = x.shape

        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        i = 0
        for train, test in cv.split(x, y):  # train and test are just for interpreter. x/y[test/train] made by cv.split.
            clf.fit(x[train], y[train])
            probabilities = clf.predict_proba(x[test])
            # Compute ROC curve and area the curve
            fpr, tpr, thresholds = roc_curve(y[test], probabilities[:, 1])
            tprs.append(interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            roc_auc = auc(fpr, tpr)
            aucs.append(roc_auc)
            plt.plot(fpr, tpr, lw=1, alpha=0.3, label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

            i += 1
        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Random', alpha=.8)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        plt.plot(mean_fpr, mean_tpr, color='b',
                 label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc), lw=2, alpha=.8)

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                         label=r'$\pm$ 1 std. dev.')

        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('%s\nROC Analysis' % variable, fontsize=20)
        plt.legend(loc="lower right")
        plt.show()

    def feature_analysis(self, variable=""):
        df_final = SVM.merge(self, variable)

        x = np.array(df_final.drop([variable], axis=1).astype(float))
        x = preprocessing.maxabs_scale(x)
        y = np.array(df_final[variable])

        clf = svm.SVC(kernel='linear', probability=True, tol=0.00001)
        rfecv = RFECV(estimator=clf, step=1, cv=StratifiedKFold(2), scoring='accuracy', verbose=0)
        rfecv.fit(x, y)
        print("Optimal number of features : %d" % rfecv.n_features_)
        df_rank = pd.DataFrame(rfecv.ranking_.transpose(),
                               index=df_final.drop([variable], axis=1).columns, columns=['Coef'])
        print(df_rank.sort_values(by=['Coef']))
        # df_rank.sort_values(by=['Coef']).plot.barh()
        # plt.show()

        clf.fit(x, y)
        df_coef = pd.DataFrame(clf.coef_.transpose(),
                               index=df_final.drop([variable], 1).columns,
                               columns=['Coef'])
        df_coef.sort_values(by=['Coef']).plot.barh()
        plt.show()

        # Plot number of features VS. cross-validation scores
        # plt.figure()
        # plt.xlabel("Number of features selected")
        # plt.ylabel("Cross validation score (nb of correct classifications)")
        # plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
        # plt.show()