.. _correlation:

***************************************************
Correlation and association on compositional data
***************************************************

Relative abundances are compositional: each sample's values share a fixed total, so one taxon going up pushes the
others down. Correlations computed directly on them are partly artefacts. Moonstone works on the centred log-ratio
(CLR) transform instead (``moonstone.transformers.compositional``), with zeros replaced explicitly beforehand.

Which analysis?
===============

``--correlation``: which features rise and fall together
    A co-occurrence network: one coefficient for every pair of features, with permutation p-values and
    Benjamini-Hochberg (BH) correction. Results can be loaded into Gephi or Cytoscape.

``--associate``: which features are linked to a clinical variable
    One linear model per feature on CLR values, adjusted for covariates (``--covariates``). With covariates this
    is a regression, not a correlation: a Spearman coefficient between a taxon and a clinical variable cannot account
    for age, treatment or any other variable. The model is MaAsLin2's linear model on CLR values, without random
    effect.

``--permanova``: whether whole communities differ between groups
    See :ref:`permanova`.

From the command line
=====================

Each run needs a count table and a metadata file; ``--second_countfile`` adds a second count table. Example files
are given below. They are synthetic (simulated values, no real samples) and small enough to run every command on
this page in a few seconds.

Input files
-----------

Species count table (:download:`species_counts.tsv <examples/species_counts.tsv>`)
    One row per taxon, one column per sample. Tab- or comma-separated: Moonstone picks the separator from the
    header line.

    .. csv-table:: First rows and columns of ``species_counts.tsv`` (60 sample columns in all)
        :header: kingdom, phylum, class, order, family, genus, species, S001, S002, S003, ...

        Bacteria, Firmicutes, Bacilli, Lactobacillales, Lactobacillaceae, Lactobacillus, Lactobacillus crispatus, 582, 579, 16398, ...
        Bacteria, Firmicutes, Bacilli, Lactobacillales, Lactobacillaceae, Lactobacillus, Lactobacillus iners, 2193, 19220, 18175, ...
        Bacteria, Actinobacteria, Actinomycetia, Bifidobacteriales, Bifidobacteriaceae, Gardnerella, Gardnerella vaginalis, 19803, 5483, 2573, ...

    - **Taxonomy columns** are named exactly ``kingdom``, ``phylum``, ``class``, ``order``, ``family``, ``genus``
      and ``species`` (lower case). They can be left out from the top (a table can start at ``genus``), but the
      level given to ``--taxa_level`` must be there. A column ``NCBI_taxonomy_ID``, if present, is ignored.
    - **Every other column is a sample.** Its header is the sample ID and must be written exactly as in the
      metadata (``S001`` and ``s001`` are different samples; ``001`` stays ``001``).
    - **Values** are counts or normalised counts: numbers only, with a dot as decimal separator and no thousands
      separator. Write 0 for a taxon not found; a blank cell stops the run. Moonstone never fills in values.
    - Rows without a name at the chosen level are pooled as ``unclassified <level>``, with a warning.

Feature table without taxonomy (:download:`resistance_genes.csv <examples/resistance_genes.csv>`)
    For genes, pathways or any other features. The **first column holds the feature names** (its header can be
    anything, here ``gene``); every other column is a sample, with the same rules for IDs and values. Feature names
    must be unique. Units do not need to match the species table (here, reads per kilobase per million): each table
    is transformed on its own.

    .. csv-table:: First rows and columns of ``resistance_genes.csv``
        :header: gene, S001, S002, S003, S004, ...

        tetM, 127.2, 43.7, 19.3, 4.0, ...
        ermB, 13.7, 0.0, 41.5, 11.3, ...

Metadata (:download:`metadata.csv <examples/metadata.csv>`)
    One row per sample, comma-separated. **The first column holds the sample IDs** (its header can be anything,
    here ``sample``); the other columns are the variables.

    .. csv-table:: First rows of ``metadata.csv``
        :header: sample, patient, Type, disease, age, antibiotics

        S001, P01, inclusion, bv, 25, no
        S002, P01, delivery, bv, 25, no
        S003, P02, inclusion, healthy, 33, no
        S004, P02, delivery, healthy, 33, no

    - Here each patient was sampled twice. ``Type`` tells the visits apart, for ``--select Type=inclusion``;
      ``patient`` names the subject, for ``--group_column patient``, which stops the run if a patient still has
      more than one sample.
    - **Text columns** (``disease``, ``antibiotics``) are groups; the reference group is chosen with
      ``--reference disease=healthy``. **Numeric columns** (``age``) are numbers: a coefficient is per unit (per
      year). A column of group codes such as 1, 2, 3 must be listed in ``--categorical``.
    - **Missing values** are left blank (or ``NA``). Samples with a missing value in a variable of the model are
      left out and counted; nothing is filled in. Do not code missing values as 0 or 99.
    - Samples in the metadata but not in the count table are ignored. Count-table samples without a metadata row
      are listed in a warning and left out.

Try it
------

From the folder holding the three example files:

.. code-block:: bash

    # co-occurrence network of species, at inclusion only
    moonstone species_counts.tsv metadata.csv results --correlation --select Type=inclusion --group_column patient

    # species against resistance genes: only pairs across the two tables
    moonstone species_counts.tsv metadata.csv results --correlation --second_countfile resistance_genes.csv \
        --select Type=inclusion

    # species associated with bacterial vaginosis, adjusted for age and antibiotics, checked with ANCOM-BC
    moonstone species_counts.tsv metadata.csv results --associate disease --covariates age antibiotics \
        --reference disease=healthy --ancombc --select Type=inclusion --group_column patient

Put ``--associate`` and ``--covariates`` after the three file arguments. ``--correlation`` and ``--associate`` can
run together; they cannot be combined with ``-f``, ``-c`` or ``--permanova``.

In the simulation, *Gardnerella vaginalis*, *Fannyhessea vaginae*, *Prevotella bivia* and *Sneathia amnii* were
made more abundant, and *Lactobacillus crispatus* less abundant, with bacterial vaginosis; ``tetM`` follows
*Gardnerella vaginalis*. The runs find these: for example *G. vaginalis* and ``tetM`` correlate at r = 0.81, and
*L. crispatus* is lower with bacterial vaginosis (coefficient -3.3 on the CLR scale, q < 0.001).

The association run also reports *Escherichia coli* and *Lactobacillus iners* as lower with bacterial vaginosis,
although their abundance was simulated with no change. This is how CLR values behave: each value is relative to
the geometric mean of its sample, so when *L. crispatus* drops sharply the others shift on the CLR scale.
ANCOM-BC, which corrects for this bias, does not report them (``ancombc_q_value`` of 0.47 and 0.81 in
``association.csv``). Associations found by both methods are the robust ones.

The example files are made by ``docs/user_guide/examples/make_examples.py``.

What happens, in this order
---------------------------

1. **Read.** Counts with taxonomy are summed to ``--taxa_level``.
2. **Select samples.** Only samples present in the metadata are kept, then those matching every
   ``--select COLUMN=VALUE`` (repeatable, all must match; values are compared as text). The run stops if a value
   is not found, and reports how many samples each selection kept. With ``--group_column`` the run stops if a
   subject still has more than one sample.
3. **Filter features**, on the selected samples: features present in at least ``--min_prevalence`` of them (10% by
   default, as in MaAsLin2), then optionally the ``--top`` N with the highest mean relative abundance. The filter
   never looks at the metadata. Filtering comes before the transform because the CLR of a feature depends on the
   other features kept. With two tables, each is filtered on its own.
4. **Transform and test.** Zeros are replaced by multiplicative replacement (each zero becomes 65% of the smallest
   non-zero proportion of its sample) and the CLR is applied, to each table on its own.
5. **Correct** the p-values over all tests of the analysis (``--correction``: ``fdr_bh`` by default,
   ``bonferroni`` or ``none``) and report results below ``--alpha`` (0.05 by default).

One sample per subject
    The analyses are cross-sectional and assume independent samples. If subjects were sampled more than once (for
    example at inclusion and at delivery), choose one time point with ``--select``. Paired analyses are not
    available yet.

Correlation options
-------------------

``--correlation_method``
    ``spearman`` (default), ``pearson`` or ``kendall`` on CLR values; ``rho`` (proportionality); ``sparcc``
    (SparCC, Friedman and Alm 2012, which needs read counts, not proportions). Across two tables only
    ``spearman``, ``pearson`` and ``kendall`` are available. Kendall is slow on large tables.
``--permutations``
    999 by default. Each permutation shuffles every feature across samples independently and transforms the table
    again; the observed coefficient is compared with the mean and spread of the shuffled ones. This is needed
    because the CLR itself makes features slightly negatively correlated, so textbook p-values give too many
    false positives (in simulations, 30% of datasets with no real correlation had a "significant" pair instead of
    5%). ``--seed`` (0 by default) makes the run repeatable.

As a guide, with 120 features and 1,000 samples (synthetic data), 999 permutations took about 30 s for Spearman and
9 s for Pearson; SparCC is slower (about 17 s for 200 permutations).

Association options
-------------------

``--covariates``
    Variables the model adjusts for; they are not reported.
``--reference VARIABLE=LEVEL``
    The level other groups are compared with (repeatable). By default, the first value in sorted order.
``--categorical``
    Numeric variables to treat as groups, for example ``1/2/3`` codes. Text and true/false columns are always
    groups. A numeric variable with 3 to 5 whole-number values gives a warning.
``--standardize``
    Scale continuous variables to a standard deviation of 1, so coefficients are per standard deviation (MaAsLin2's
    default). By default, coefficients are per unit (for example per year of age).
``--ancombc``
    Run ANCOM-BC (scikit-bio) on the same model as a sensitivity check. Its log-fold changes are not on the CLR
    scale: compare directions and significance, not values.

Samples with a missing value in any variable of the model are left out (complete cases); the number is reported.

Output files
------------

``correlation_edges.csv``
    One row per pair: ``feature_a``, ``feature_b``, ``method``, ``r``, ``z`` (distance from the permutation null
    in standard deviations), ``p_value``, ``q_value``, ``n_samples``, ``prevalence_a``, ``prevalence_b``.
``correlation_gephi_edges.csv``
    The significant pairs, positive and negative, with the columns Gephi and Cytoscape import directly
    (``Source``, ``Target``, ``Weight`` = absolute r, ``Sign``).
``correlation_matrix.csv``
    The r matrix, for a heatmap.
``association.csv``
    One row per feature and variable level: ``coef``, ``std_error``, ``ci_low`` and ``ci_high`` (95% confidence
    interval), ``t_value``, ``p_value``, ``q_value``, ``n_samples``, ``n_dropped``, ``prevalence``, and with
    ``--ancombc`` the ANCOM-BC q-value and log-fold change (also in ``ancombc.csv``).
``methods.txt``
    A paragraph per analysis for the methods section of a paper: Moonstone version, samples and selections,
    filters, zero handling, transform, method, permutations and seed, correction and threshold. Check it against
    what you report.

From Python
===========

.. code-block:: python

    from moonstone.analysis.association import feature_association
    from moonstone.analysis.correlation import feature_correlation, filter_edges, to_gephi

    edges = feature_correlation(species_df)                      # features as rows, samples as columns
    network = to_gephi(filter_edges(edges, alpha=0.05))
    results = feature_association(species_df, metadata_df, 'disease', covariates=['age'])

What is not in Moonstone
========================

FastSpar
    FastSpar (Watts et al., 2019; installable from bioconda) is a fast C++ implementation of SparCC. Moonstone's
    SparCC is in Python and does not call FastSpar. For very large tables FastSpar is much faster; its results
    should be close to Moonstone's but are not identical (random draws differ).
Models that stay in R
    Hurdle and zero-inflated mixed models (``glmmTMB``), MaAsLin2's compound Poisson models (CPLM), ``emmeans``
    contrasts and ``DHARMa`` diagnostics have no good Python equivalent and are not ported.
Random effects
    Association models have no random effect per subject; analyses are cross-sectional.

The association model is checked against R's ``lm()`` on the same synthetic data (``tests/reference``).
