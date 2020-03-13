import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interp
from sklearn import preprocessing, svm
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_curve, auc
from sklearn.feature_selection import RFECV
from moonstone.analysis import stats
from moonstone.normalization.processed.scaling_normalization import StandardScaler
from moonstone.utils.df_merge import MergeDF

logger = logging.getLogger(__name__)


class ML(object):
    def __init__(self, countfile, metadata, outdir):
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.countfile = countfile
        self.metadata = metadata
        self.outdir = outdir

    def svm(self, variable=""):
        logger.info(f'Starting SVM analysis with variable: {variable}')

        # Collect metadata and count data in a single DataFrame
        merger = MergeDF(self.countfile, self.metadata, variable)
        df_final = merger.merged_df

        # Setup the features and labels
        x = np.array(df_final.drop([variable], axis=1).astype(float))
        y = np.array(df_final[variable])

        # Scale the Features
        scaler = StandardScaler(x)
        x_scaled = scaler.scaled_x

        # My variable counter
        stats.count_items(y)

        # Support Vector Machine
        print('\nRunning Support Vector Machine classifier. This could take a few moments...')
        for kernel in ["linear", "poly", "rbf", "sigmoid"]:
            clf = svm.SVC(kernel=kernel, tol=0.00001, gamma='scale')
            # Get accuracy
            score = cross_val_score(clf, x_scaled, y, cv=10)
            print("\tAccuracy for %s kernel: %.2f%s (+/- %.2f%s)" %
                  (kernel, score.mean()*100, "%", score.std()*200, "%"))

    def roc_analysis(self, variable=""):
        logger.info(f'Starting ROC analysis with variable: {variable}')
        merger = MergeDF(self.countfile, self.metadata, variable)
        df_final = merger.merged_df

        x = np.array(df_final.drop([variable], axis=1).astype(float))
        x = preprocessing.maxabs_scale(x)
        y = np.array(df_final[variable])

        # Here we make sure 2, and only 2 categories are present.
        if len(set(y)) != 2:
            logger.fatal('ROC analysis attempted with more than 2 variable states!')
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
        # Save figure to PDF
        plt.savefig(self.outdir+"/"+variable+'_ROC.pdf', format='pdf', dpi=150)

    def feature_analysis(self, variable=""):
        logger.info(f'Running SVM feature importance analysis with variable: {variable}')
        merger = MergeDF(self.countfile, self.metadata, variable)
        df_final = merger.merged_df
        logger.info('Feature Analysis successfully collected merged database.')

        # Setup the data
        x = np.array(df_final.drop([variable], axis=1).astype(float))
        x = preprocessing.maxabs_scale(x)
        y = np.array(df_final[variable])

        # Calculate and plot the contributions of all variables
        logger.info('Running SVM with linear kernel and reporting coefficients')
        clf = svm.SVC(kernel='linear', probability=True, tol=0.00001)

        clf.fit(x, y)
        df_coef = pd.DataFrame(clf.coef_.transpose(),
                               index=df_final.drop([variable], 1).columns,
                               columns=['Coef'])
        df_coef.sort_values(by=['Coef']).plot.barh()
        plt.savefig(self.outdir + '/' + variable + '_variable_coefs.pdf', format='pdf', dpi=150)
        logger.info('Coefficients Figure for SVM Linear kernel saved to ' + variable + '_variable_coefs.pdf')

        # Run and plot Recursive Factor Elimination with Cross Validation
        rfecv = RFECV(estimator=clf, step=1, cv=StratifiedKFold(10), scoring='accuracy', verbose=0)

        rfecv.fit(x, y)
        logger.info("Optimal number of features : %d" % rfecv.n_features_)

        df_rank = pd.DataFrame(rfecv.ranking_.transpose(),
                               index=df_final.drop([variable], axis=1).columns, columns=['Coef'])

        logger.info('Full component list for %s being written to %s'
                    % (variable, self.outdir+'/'+variable+'_SVM_components.csv'))

        df_rank.sort_values(by=['Coef'], ascending=[bool])\
            .to_csv(path_or_buf=self.outdir+'/'+variable+'_SVM_components.csv', sep=',')

        plt.figure()
        plt.xlabel("Number of features selected")
        plt.ylabel("Cross validation score (nb of correct classifications)")
        plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_)
        plt.savefig(self.outdir + '/' + variable + '_rfecv.pdf', format='pdf', dpi=150)
        logger.info('RFECV Figure for SVM Linear kernel saved to ' + variable + '_rfecv.pdf')
        logger.info('SVM feature importance analysis for %s completed' % variable)
