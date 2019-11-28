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
module_logger.info('Using randomForest module.')


class RandomForest(object):
    def __init__(self, countfile, metadata, outdir, variable=""):
        self.logger = module_logger
        self.logger.info('Starting instance of RandomForest.')
        self.countfile = countfile
        self.metadata = metadata
        self.variable = variable
        self.outdir = outdir

    # Can use the merging function in 'classify'
    def get_matrix(self):
        merged_df = classify.SVM(self.countfile, self.metadata)
        df = merged_df.merge(self.variable)
        return df

    def forest(self, filename):
        df = RandomForest.get_matrix(self)

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

        # Setup the different classifiers. Decision tree and a 'forect of tress are obvious choices.
        # Also included are two type of boosting.
        dt = DecisionTreeClassifier()
        rf = RandomForestClassifier(n_estimators=100, max_features="auto", random_state=33)
        ab = AdaBoostClassifier(n_estimators=100, random_state=33)
        gb = GradientBoostingClassifier(n_estimators=100, random_state=33)

        # We will also test the two main types of Cross Validation
        sf = StratifiedKFold(n_splits=10, random_state=33, shuffle=True)
        tts = train_test_split(x, y, test_size=.25, shuffle=True)

        # Decision Trees, TTS and then combined results of samples folds
        print("Running Analysis...\nDecision Tree:")
        x_train, x_test, y_train, y_test = tts
        dt.fit(x_train, y_train)
        y_predict_tts = dt.predict(x_test)
        print("\tTrain/Test Split Accuracy: %2.1f%s" % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        i = 0
        score = 0
        for train, test in sf.split(x, y):
            dt.fit(x[train], y[train])
            y_predict = dt.predict(x[test])
            score += accuracy_score(y[test], y_predict)
            i += 1
        print("\tSample Folds Accuracy: %2.1f%s" % (score / i * 100, "%"))

        print("Random Forest:")
        # Random Forest, TTS and then combined results of samples folds
        # We don't need to generate the train/test data again.
        rf.fit(x_train, y_train)
        y_predict_tts = rf.predict(x_test)
        print("\tTrain/Test Split Accuracy: %2.1f%s" % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        scores = cross_val_score(rf, x, y, cv=10)
        print("\t10X (Stratified)KFold Accuracy: %0.2f%s (+/- %0.2f)" % (scores.mean() * 100, "%", scores.std() * 200))

        # AdaBoost, TTS and then combined results of samples folds
        print("AdaBoost:")
        ab.fit(x_train, y_train)
        y_predict_tts = ab.predict(x_test)
        print("\tTrain/Test Split Accuracy: %2.1f%s" % (accuracy_score(y_test, y_predict_tts) * 100, "%"))

        i = 0
        score = 0
        for train, test in sf.split(x, y):
            ab.fit(x[train], y[train])
            y_predict = ab.predict(x[test])
            score += accuracy_score(y[test], y_predict)
            i += 1
        print("\tSample Folds Accuracy: %2.1f%s" % (score / i * 100, "%"))

        # Gradient Boost, TTS and then combined results of samples folds
        print("Gradient Boost:")
        ab.fit(x_train, y_train)
        y_predict_tts = ab.predict(x_test)
        print("\tTrain/Test Split Accuracy: %2.1f%s" % (accuracy_score(y_test, y_predict_tts) * 100, "%"))
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
        plt.title(f"Random Forest Feature Importance: {self.variable}", fontsize=16)
        plt.bar(range(features_to_print), importance[indices][:features_to_print], color="r",
                yerr=standard_deviations[indices][:features_to_print], bottom=.26)
        plt.xticks(range(features_to_print), features, rotation=90)
        plt.xlim([-1, features_to_print])
        plt.show()
