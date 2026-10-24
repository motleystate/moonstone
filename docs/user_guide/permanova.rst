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

From the command line
=====================

Start from a normalized count table with one row per taxon, the taxonomy in columns named ``kingdom``, ``phylum``,
``class``, ``order``, ``family``, ``genus`` and ``species``, and one column per sample (tab- or comma-separated),
plus a metadata file with one row per sample:

.. code-block:: bash

    # every metadata variable, species level
    moonstone normalized_counts.tsv metadata.csv results --permanova

    # two variables, every taxonomic level
    moonstone normalized_counts.tsv metadata.csv results --permanova treatment site --taxa_level all

Put ``--permanova`` after the three file arguments. The results go to ``results/permanova.csv``, one row per level
and variable. Other options: ``--no_read_threshold``, ``--permutations`` (9,999 by default), ``--seed``
(0 by default) and ``--threads`` (number of tests run at the same time, one per processor core; 0, the default, uses
every core available to Moonstone). ``--permanova`` cannot be combined with ``-f``, ``-p``, ``-k`` or ``-c``.

Most of the run time is the PERMDISP dispersion check. As a guide, measured on synthetic data with one core
(unverified on real data): one variable at 9,999 permutations took about 8 s with 100 samples, 18 s with 300 and
2 minutes with 1,000. ``--threads`` runs the variables (and levels) in parallel, so 4 cores make a screen of many
variables about 3.5 times faster. The results are identical whatever the number of cores, because every test uses
the same seed whichever core runs it. From Python, pass ``n_jobs`` to ``run_permanova`` or ``permanova_screen``
(1 by default).

The same analysis from Python:

.. code-block:: python

    from moonstone.analysis.permanova import read_taxonomy_counts, run_permanova

    counts = read_taxonomy_counts("normalized_counts.tsv")
    results = run_permanova(counts, metadata_df, taxa_level="all")

For each taxonomic level:

1. The counts are summed to that level. Summing comes first, so a genus made of many rare species is judged on all
   of its reads.
2. Unless ``--no_read_threshold`` is given, only the taxa above the knee of the read-threshold score
   (prevalence x log10(mean reads + 1), see ``moonstone.filtering.read_threshold_analysis``) are kept. This
   filter never looks at the groups, so it does not bias the test.
3. PERMANOVA and PERMDISP are run for each variable.

The p-values of every test in the run (all levels and variables) are corrected together.

The ``warnings`` column
    The read-threshold knee has been checked at species level only, so every other level carries a warning
    (see https://github.com/motleystate/moonstone/issues/132). A warning is also added when the filter removes more
    than 10% of the reads, when taxa without a name at that level are pooled as ``unclassified``, and when samples
    with no reads are left out. Read the warnings before the p-values.

Levels with a single taxon
    With one taxon (often the kingdom), Bray-Curtis only compares total reads, not community composition, so the
    level is skipped.

How the test works
==================

PERMANOVA compares the distances between groups with the distances within groups, and summarises them in a
pseudo-F statistic (like the F of an ordinary analysis of variance). To get a p-value, the group labels are shuffled
many times (9,999 by default) and the test counts how often a shuffled pseudo-F is at least as large as the real
one: p = (1 + number at least as large) / (1 + number of shuffles).
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

Precision of the p-values
    A permutation p-value is an estimate. ``p_value_se`` and ``permdisp_p_value_se`` give its standard error,
    sqrt(p (1 - p) / N) for N permutations: rerunning with another seed would usually move the p-value by less than
    about twice this. With 9,999 permutations, a p-value of 0.05 has a standard error of about 0.002; with 999, about
    0.007. The smallest possible p-value is 1 / (N + 1), for example 0.0001 with 9,999 permutations.

Enough permutations for the correction
    With m tests corrected together, a single real effect can only reach a corrected p-value of m / (N + 1). If
    that is not below 0.05, nothing can be significant however strong the effect, so a warning is logged and added
    to the ``warnings`` column with the number of permutations needed. With 9,999 permutations this happens beyond
    499 tests; with 999, beyond 49 (for example 50 variables at one level, or 7 levels x 8 variables).

Reproducibility
    The shuffles are seeded (``seed=0`` by default), so the same call gives the same p-values. The seed gives
    repeatability, not accuracy: a different seed gives p-values that differ by about their standard error, and
    neither is more correct. More permutations, not a particular seed, make the p-values more accurate. Report the
    number of permutations and the seed with the results.

Which samples are used
    Samples without metadata and groups with fewer than ``min_group_size`` samples (2 by default) are left out, and
    this is logged. ``permanova_screen`` skips variables with fewer than 2 usable groups or more than ``max_groups``
    (15 by default) distinct values, which are usually continuous variables, logs a warning, and says why in the
    ``skipped`` column.

Missing values
    A missing value is treated as unknown, not as a group of its own. Each variable is tested on the samples that
    have a value for it (complete-case), so different variables can use different samples. PERMANOVA and PERMDISP
    always use the same samples. ``n_missing`` gives the number of samples left out for a missing value, and
    ``n_samples`` the number tested; the count is also logged, and printed by ``--permanova``.

    Missing means an empty cell, a cell holding only spaces, or a value pandas reads as missing (``NA``, ``N/A``,
    ``NaN``, ``null`` and similar). Words such as ``missing``, ``unknown`` or ``-`` are real values and form a
    group, so leave those cells empty instead.

    Missing values are never filled in or recoded. If a blank actually means something (for example "not
    detected"), recode it in the metadata before the analysis. When more than half of the samples have no value
    for a variable, a warning is logged and ``--permanova`` flags the variable: check that the samples with a value
    are not a special subset (for example only the patients seen at one clinic).

Limits
======

* Samples must be independent. If a subject gave several samples (time points, body sites), the shuffles should
  keep each subject's samples together; this is not supported yet, so use one sample per subject.
* Variables are treated as categories. Continuous variables (age, for example) are not tested.
* PERMDISP uses the Principal Coordinates Analysis axes with positive eigenvalues. For non-Euclidean distances such
  as Bray-Curtis, R's ``vegan::betadisper`` also uses the negative ones, so PERMDISP values can differ slightly from R.
