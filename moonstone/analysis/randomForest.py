import logging
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.model_selection import StratifiedKFold, train_test_split, cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, GradientBoostingClassifier

from moonstone.analysis import classify

module_logger = logging.getLogger(__name__)


class RandomForest(object):
    def __init__(self, countfile, metadata, outdir, variable=""):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.countfile = countfile
        self.metadata = metadata
        self.variable = variable
        self.outdir = outdir

    # Can use the merging function in 'classify'
    def get_matrix(self):
        self.logger.info('Building the merged data-frame from count data and selected variable')
        merged_df = classify.SVM(self.countfile, self.metadata, self.outdir)
        df = merged_df.merge(self.variable)
        self.logger.info('Done. Returning merged dataframe.')
        return df

    def forest(self, filename):
        self.logger.info('Starting random forest analysis.')
        df = RandomForest.get_matrix(self)

        self.logger.info('Setting up variables.')
        x = np.array(df.drop([self.variable], axis=1))
        x = preprocessing.maxabs_scale(x)
        y = np.array(df[self.variable])
        c = Counter(y)
        c = dict(c)

        sample_count = 0
        for _, value in c.items():
            sample_count += value

        print(f"There appear to be a total of {sample_count} samples.")
        for category, value in c.items():
            print(" %s samples labeled at %s, or %2.1f%s of the total"
                  % (value, category, value/sample_count*100, "%"))

        # Setup the different classifiers. Decision tree and a 'forest' of tress are obvious choices.
        # Also included are two type of boosting.
        dt = DecisionTreeClassifier()
        rf = RandomForestClassifier(n_estimators=100, max_features="auto", random_state=33)
        ab = AdaBoostClassifier(n_estimators=100, random_state=33)
        gb = GradientBoostingClassifier(n_estimators=100, random_state=33)

        # We will also test the two main types of Cross Validation
        sf = StratifiedKFold(n_splits=10, random_state=33, shuffle=True)
        tts = train_test_split(x, y, test_size=.25, shuffle=True)

        # Decision Trees, TTS and then combined results of samples folds
        self.logger.info("Running Decision Tree Analysis.")
        x_train, x_test, y_train, y_test = tts
        dt.fit(x_train, y_train)
        y_predict_tts = dt.predict(x_test)
        print("Decision Tree:\n\tTrain/Test Split Accuracy: %2.1f%s"
              % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        i = 0
        score = 0
        for train, test in sf.split(x, y):
            dt.fit(x[train], y[train])
            y_predict = dt.predict(x[test])
            score += accuracy_score(y[test], y_predict)
            i += 1
        print("\tSample Folds Accuracy: %2.1f%s" % (score / i * 100, "%"))

        self.logger.info("Running Random Forest Analysis")
        # Random Forest, TTS and then combined results of samples folds
        # We don't need to generate the train/test data again.
        rf.fit(x_train, y_train)
        y_predict_tts = rf.predict(x_test)
        print("Random Forest:\n\tTrain/Test Split Accuracy: %2.1f%s"
              % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        scores = cross_val_score(rf, x, y, cv=10)
        print("\t10X (Stratified)KFold Accuracy: %0.2f%s (+/- %0.2f)" % (scores.mean() * 100, "%", scores.std() * 200))

        # AdaBoost, TTS and then combined results of samples folds
        self.logger.info("Running AdaBoost Analysis")
        ab.fit(x_train, y_train)
        y_predict_tts = ab.predict(x_test)
        print("AdaBoost:\n\tTrain/Test Split Accuracy: %2.1f%s" % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        i = 0
        score = 0
        for train, test in sf.split(x, y):
            ab.fit(x[train], y[train])
            y_predict = ab.predict(x[test])
            score += accuracy_score(y[test], y_predict)
            i += 1
        print("\tSample Folds Accuracy: %2.1f%s" % (score / i * 100, "%"))

        # Gradient Boost, TTS and then combined results of samples folds
        self.logger.info("Running Gradient Boost Analysis")
        ab.fit(x_train, y_train)
        y_predict_tts = ab.predict(x_test)
        print("Gradient Boost:\n\tTrain/Test Split Accuracy: %2.1f%s"
              % (accuracy_score(y_test, y_predict_tts) * 100, "%"))
        i = 0
        score = 0
        for train, test in sf.split(x, y):
            gb.fit(x[train], y[train])
            y_predict = gb.predict(x[test])
            score += accuracy_score(y[test], y_predict)
            i += 1
        print("\tSample Folds Accuracy: %2.1f%s" % (score / i * 100, "%\n"))

        # Here is the output if features in the RF classifier. Just for vanilla RF.
        # Screen output is limited to the top 100 features.
        importance = rf.feature_importances_
        indices = np.argsort(importance)[::-1]
        features = df.drop([self.variable], axis=1).columns[indices]
        standard_deviations = np.std([tree.feature_importances_ for tree in rf.estimators_], axis=0)
        print("%s informative features out of %s total..." % (np.count_nonzero(importance), importance.shape[0]))

        df_complete_list = pd.DataFrame(zip(features, importance[indices], standard_deviations[indices]),
                                        columns=('features', 'importance', 'std'))
        output_file = self.outdir+'/'+filename
        df_complete_list.to_csv(path_or_buf=output_file, sep=',')

        if np.count_nonzero(importance) > 100:
            features_to_print = 100
            print("Showing the top 100! Full list can be found in csv file.")
        else:
            features_to_print = np.count_nonzero(importance)
            print("Showing them all. Also writing to csv file.")
        for f in range(features_to_print):
            print("%d. feature %s (%f)" % (
                f + 1, features[f], importance[indices[f]]))

        plt.figure()
        plt.title(f"Random Forest Feature Importance: {self.variable}", fontsize=12)
        plt.bar(range(features_to_print), importance[indices][:features_to_print], color="r",
                yerr=standard_deviations[indices][:features_to_print])
        plt.xticks(range(features_to_print), features, rotation=90, fontsize=3)
        plt.xlim([-1, features_to_print])
        # plt.show()
        plt.savefig(self.outdir+"/"+self.variable+'_rfFeatures.pdf', format='pdf', dpi=150)
