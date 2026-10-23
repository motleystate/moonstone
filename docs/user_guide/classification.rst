.. _classification:

**************
Classification
**************

Classification asks whether the microbial profile of a sample can predict a metadata variable, for example
disease status. Moonstone trains a Support Vector Machine (SVM), a random forest, or both, and reports how well
they predict the variable for samples they were not trained on.

How the score is obtained
=========================

The steps below are what makes the reported score trustworthy. Each one closes a way in which a classifier can
look better than it really is.

Cross-validation
    The samples are cut into *k* parts (folds, 5 by default). The model is trained on *k - 1* folds and tested on
    the remaining one, until every fold has been the test fold once. Only scores on test folds are reported.

Stratified and seeded folds
    Each fold keeps the proportions of the groups found in the whole data set, and the random assignment of samples
    to folds is fixed by a seed (``--seed``, 0 by default), so the same command gives the same numbers.

Grouped folds
    When one subject gave several samples (time points, body sites), ``--group_column`` keeps all samples of a
    subject in the same fold. Otherwise the model can recognise the subject instead of the condition.

Scaling inside the training folds
    For the SVM, each taxon is scaled to mean 0 and standard deviation 1. The scaling is learnt from the training
    folds only and then applied to the test fold, so nothing about the test samples leaks into training. The random
    forest is not scaled, because decision trees are not affected by the scale of a taxon.

Unequal group sizes
    Errors in a small group are weighted as heavily as errors in a large group during training (class weights).
    The main score is **balanced accuracy**, the average of the per-group hit rates. Its chance level is
    1 / number of groups, whatever the group sizes, and it is printed next to the result. Plain accuracy is
    reported as well.

Small groups
    Every test fold needs at least one sample of each group, so the number of folds is lowered, with a warning,
    when the smallest group has fewer samples (or subjects) than ``--folds``. A group with a single sample cannot
    be cross-validated and stops the analysis with an error.

Area under the curve
    With ``--roc``, Receiver Operating Characteristic (ROC) curves are drawn for the test folds, and the Area Under
    the Curve (AUC) is reported: 1 means perfect separation and 0.5 means guessing. With two groups, the second
    group in alphabetical order is the "positive" group. With three or more groups, each group is compared with all
    the others ("one-vs-rest") and the AUC is averaged over groups.

Feature selection inside the cross-validation
    With ``--feature_selection``, Recursive Feature Elimination with Cross-Validation (RFECV) repeatedly drops the
    least useful taxon and keeps the number of taxa that scored best. It runs inside each training part, with its
    own inner folds, and the outer test fold scores the result (nested cross-validation). Selecting taxa on all
    samples and then scoring on the same samples would be circular. The SVM needs the linear kernel for this.

Permutation importance
    To measure how much a model relies on a taxon, the values of that taxon are shuffled among the test samples
    and the loss of balanced accuracy is recorded. A large loss means the model depends on that taxon. It is
    computed on test folds for both models, and replaces the coefficient plot and the built-in random forest
    importance of earlier versions.

Choosing the SVM kernel
=======================

The kernel sets the shape of the boundary between groups: ``linear`` (a flat boundary, the default),
``poly`` (polynomial), ``rbf`` (radial basis function) or ``sigmoid``. There are two ways to choose it.

``--explore_kernels``
    Scores the four kernels on the same folds and writes ``VARIABLE_svm_kernel_comparison.csv``, best first.
    These scores are **exploratory**: the best of four scores is optimistic, because part of its lead is luck on
    these particular samples. Use it to see whether the kernel matters, not to quote a final score.

``--kernel auto``
    Picks the kernel inside each training part with an inner cross-validation, and scores that choice on the outer
    test fold. The score reported is a fair estimate for the whole procedure, including the choice of kernel. The
    kernel chosen in each fold is printed and saved in ``cv_scores.csv``.

Command line
============

.. code-block:: bash

    # Support Vector Machine with the default linear kernel
    $ moonstone counts.csv metadata.csv output --classify disease

    # Both models with ROC curves, samples of one patient kept together
    $ moonstone counts.csv metadata.csv output --classify disease --model svm random_forest --roc \
        --group_column patient_id

    # Compare kernels (exploratory), then let nested cross-validation choose one
    $ moonstone counts.csv metadata.csv output --classify disease --explore_kernels --kernel auto

    # Nested feature selection with the linear SVM
    $ moonstone counts.csv metadata.csv output --classify disease --feature_selection

Options
"""""""

=============================  ==========================================================================
Option                         Meaning
=============================  ==========================================================================
``-c, --classify VARIABLE``    Metadata column holding the groups to predict.
``--model``                    ``svm``, ``random_forest`` or both (default ``svm``).
``--kernel``                   ``linear`` (default), ``poly``, ``rbf``, ``sigmoid`` or ``auto``.
``--explore_kernels``          Compare the four SVM kernels (exploratory scores).
``--folds N``                  Number of cross-validation folds (default 5).
``--seed N``                   Random seed for folds, random forest and permutations (default 0).
``--group_column COLUMN``      Metadata column naming the subject of each sample.
``--roc``                      Plot ROC curves on the test folds.
``--feature_selection``        Nested recursive feature elimination.
=============================  ==========================================================================

Output files
""""""""""""

Files are written to the output directory, named ``VARIABLE_MODEL_...`` where ``MODEL`` is ``svm`` or
``random_forest``.

=========================================  ==================================================================
File                                       Content
=========================================  ==================================================================
``_cv_scores.csv``                         Balanced accuracy, accuracy and AUC of each test fold; number of
                                           taxa used and, for the SVM, the kernel.
``_permutation_importance.csv``            Mean and standard deviation over folds of the loss of balanced
                                           accuracy when each taxon is shuffled, highest first.
``_permutation_importance.pdf``            Bar plot of the 30 most important taxa.
``_roc.pdf``                               ROC curves (with ``--roc``).
``_selected_features.csv``                 How many folds kept each taxon (with ``--feature_selection``).
``_rfecv.pdf``                             Inner balanced accuracy against number of taxa kept (with
                                           ``--feature_selection``). The fair score is in ``_cv_scores.csv``.
``VARIABLE_svm_kernel_comparison.csv``     Exploratory kernel comparison (with ``--explore_kernels``).
=========================================  ==================================================================

Replaced flags
""""""""""""""

The classifier flags of earlier versions were removed. Using one of them stops with an error that names its
replacement.

==================================  ==================================================
Removed flag                        Replacement
==================================  ==================================================
``-s, --svm VARIABLE``              ``--classify VARIABLE`` (add ``--explore_kernels`` to compare kernels)
``-sr, --svm_roc VARIABLE``         ``--classify VARIABLE --roc``
``-sc, --svm_classifier VARIABLE``  ``--classify VARIABLE --feature_selection``
``-rf, --random_forest VARIABLE``   ``--classify VARIABLE --model random_forest``
==================================  ==================================================

The old output files ``rf_Allfeatures.csv``, ``VARIABLE_rfFeatures.pdf``, ``VARIABLE_variable_coefs.pdf``,
``VARIABLE_SVM_components.csv`` and ``VARIABLE_ROC.pdf`` are no longer written. The decision tree, AdaBoost and
gradient boosting comparisons, and the single train/test split, were removed.

In Python
=========

.. code-block:: python

    from moonstone.analysis.classification import Classification

    # count_df: one row per sample, one column per taxon; metadata_df: one row per sample
    classification = Classification(count_df, metadata_df, "disease", "output", model="random_forest",
                                    group_column="patient_id")
    results = classification.run(roc=True)
    results["folds"]  # per-fold scores, as in disease_random_forest_cv_scores.csv

Feature selection with the random forest refits many forests and can be slow on tables with hundreds of taxa.
``Classification(..., rfe_step=0.05)`` removes 5 % of the taxa at each round instead of one taxon, which is much
faster. See :class:`moonstone.analysis.classification.Classification` for all parameters.
