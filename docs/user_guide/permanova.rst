.. _permanova:

*****************************************
Comparing groups by community (PERMANOVA)
*****************************************

Permutational Multivariate Analysis of Variance (PERMANOVA) asks whether samples from different groups (for example
treated and untreated) differ in their microbial community as a whole. It works on a beta-diversity distance
matrix, so it is available on every beta-diversity class (``BrayCurtis``, ``Jaccard``, ``WeightedUniFrac``,
``UnweightedUniFrac``).

.. code-block:: python

    from moonstone.analysis.diversity.beta import BrayCurtis

    bray_curtis = BrayCurtis(normalized_counts_df)

    # one variable
    result = bray_curtis.permanova(metadata_df, "treatment")

    # many variables, with correction for multiple testing
    screen = bray_curtis.permanova_screen(metadata_df)

How the test works
==================

PERMANOVA compares the distances between groups with the distances within groups, and summarises them in a
pseudo-F statistic (like the F of an ordinary analysis of variance). To get a p-value, the group labels are shuffled
many times (999 by default) and the test counts how often a shuffled pseudo-F is at least as large as the real one.
``R2`` is the share of the variation between samples that the groups explain.

Reading the results
===================

Dispersion check (PERMDISP)
    A significant PERMANOVA means the groups differ in their centre, in their spread, or both. The Permutational
    Analysis of Multivariate Dispersions (PERMDISP) tests the spread alone. If ``permdisp_p_value`` is also small,
    part or all of the PERMANOVA result may come from one group being more variable than the other, rather than
    from a different composition. Look at the PCoA plot (``visualize_pcoa``) before concluding.

Multiple testing
    Testing many metadata variables at p < 0.05 gives about one false positive per 20 variables with no real
    effect. ``permanova_screen`` corrects the p-values with the Benjamini-Hochberg false discovery rate (FDR) by
    default (``correction_method="fdr_bh"``); read ``p_value_corrected`` and ``permdisp_p_value_corrected``.

Reproducibility
    The shuffles are seeded (``seed=0`` by default), so the same call gives the same p-values.

Which samples are used
    Samples without metadata, samples with a missing value in the tested variable, and groups with fewer than
    ``min_group_size`` samples (2 by default) are left out, and this is logged. ``permanova_screen`` skips variables
    with fewer than 2 usable groups or more than ``max_groups`` (15 by default) distinct values, which are usually
    continuous variables, and says why in the ``skipped`` column.

Limits
======

* Samples must be independent. If a subject gave several samples (time points, body sites), the shuffles should
  keep each subject's samples together; this is not supported yet, so use one sample per subject.
* Variables are treated as categories. Continuous variables (age, for example) are not tested.
* PERMDISP uses the Principal Coordinates Analysis axes with positive eigenvalues. For non-Euclidean distances such
  as Bray-Curtis, R's ``vegan::betadisper`` also uses the negative ones, so PERMDISP values can differ slightly from R.
