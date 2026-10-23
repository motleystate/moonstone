"""Probe each 'to confirm' audit finding against installed libraries. Does not modify repo code."""
import os, sys, tempfile, traceback, warnings
import matplotlib; matplotlib.use("Agg")
import numpy as np, pandas as pd, scipy, sklearn
warnings.filterwarnings("ignore")
print(f"python {sys.version.split()[0]} pandas {pd.__version__} sklearn {sklearn.__version__} scipy {scipy.__version__} numpy {np.__version__}")

def run(name, fn):
    try:
        fn(); print(f"[{name}] NO ERROR")
    except BaseException as e:
        tb = traceback.extract_tb(e.__traceback__)
        where = next((f"{os.path.basename(f.filename)}:{f.lineno}" for f in reversed(tb) if "moonstone" in f.filename), f"{os.path.basename(tb[-1].filename)}:{tb[-1].lineno}")
        print(f"[{name}] {type(e).__name__} at {where}: {str(e).splitlines()[0][:200]}")

# 1. scipy.interp import (also blocks `moonstone --help`)
run("A classify.py:5 import", lambda: __import__("moonstone.analysis.classify"))

# Shim so the rest of classify.py can be reached (scratch only, not a fix)
scipy.interp = np.interp
from moonstone.analysis import classify, randomForest

rng = np.random.default_rng(0)
n, p = 40, 12
idx = [f"s{i}" for i in range(n)]
counts = pd.DataFrame(rng.poisson(20, (n, p)).astype(float), index=idx, columns=[f"taxon{j}" for j in range(p)])
meta = pd.DataFrame({"group": ["A"] * 20 + ["B"] * 20}, index=idx)
out = tempfile.mkdtemp()

# Library-level checks, independent of call order
run("B lib RandomForestClassifier(max_features='auto')",
    lambda: sklearn.ensemble.RandomForestClassifier(max_features="auto").fit(counts.values, meta.group.values))
def rfecv_attr():
    r = sklearn.feature_selection.RFECV(sklearn.svm.SVC(kernel="linear"), cv=3).fit(counts.values, meta.group.values)
    r.grid_scores_
run("C lib RFECV.grid_scores_", rfecv_attr)
run("D lib DataFrame.drop([col], 1)", lambda: counts.drop(["taxon0"], 1))
run("E lib sort_values(ascending=[bool])",
    lambda: pd.DataFrame({"Coef": [3, 1, 2]}).sort_values(by=["Coef"], ascending=[bool]))

# End-to-end through moonstone's own methods
ml = lambda: classify.ML(counts.copy(), meta.copy(), out)
run("F ML.svm", lambda: ml().svm(variable="group"))
run("G ML.roc_analysis", lambda: ml().roc_analysis(variable="group"))
run("H ML.feature_analysis", lambda: ml().feature_analysis(variable="group"))
rf = lambda: randomForest.RandomForest(counts.copy(), meta.copy(), out, variable="group")
run("I RandomForest.forest (as shipped)", lambda: rf().forest("rf.csv"))
# Bypass broken get_matrix to see whether line 58 is reached next
def forest_bypass():
    from moonstone.utils.pandas.df_merge import MergeDF
    r = rf(); r.get_matrix = lambda: MergeDF(counts.copy(), meta.copy(), "group").merged_df
    randomForest.RandomForest.get_matrix = lambda self: MergeDF(counts.copy(), meta.copy(), "group").merged_df
    r.forest("rf.csv")
run("J RandomForest.forest (get_matrix bypassed)", forest_bypass)
