.. _changelog:

**********
Changelogs
**********

Summary of developments of moonstone library.


Unreleased
==========

* Requires Python 3.13 and pandas >= 3.0 (support for Python 3.9 to 3.12 and for pandas 2 is dropped).
  Minimum scikit-bio raised to 0.6.0. (https://github.com/motleystate/moonstone/pull/118)
* Compatibility with pandas 3 (https://github.com/motleystate/moonstone/pull/118)
* Command-line tool starts again with recent SciPy versions (``scipy.interp`` replaced by ``numpy.interp``);
  ``moonstone --help`` now exits with status 0. New dependency: ``filetype>=1.2.0``. (https://github.com/motleystate/moonstone/pull/117)
* Install instructions updated for Python 3.13; install from GitHub until a release newer than 1.0.0 is on PyPI.
* **Breaking change:** the command-line flags ``-s/--svm``, ``-sr/--svm_roc``, ``-sc/--svm_classifier`` and
  ``-rf/--random_forest`` are removed and stop with an error naming their replacement, ``--classify VARIABLE``
  (with ``--roc``, ``--feature_selection`` or ``--model random_forest``). Output file names change, and the decision
  tree, AdaBoost and gradient boosting runs are removed. See :ref:`classification`.
  (https://github.com/motleystate/moonstone/pull/123)
* Fix: ``-f/--filtering`` on the command line now removes taxa whose mean read count is below the threshold,
  as documented. It used to remove samples instead, which usually stopped the run with a sample-number mismatch.
* Fix: Faith's phylogenetic diversity and UniFrac now sum a multi-level taxonomy table to its last level as intended.
  Before, such a table raised a misleading "INCOMPLETE TREE" error, or gave 0 for every sample with
  ``force_computation=True``.
* Fix: ``GeometricMeanNormalization`` raises an error naming samples without reads, instead of looping forever.
* Fix: the command-line count-file warning gives the real number of unnamed taxa removed.
* Fixes to FASTQ read downsizing (``moonstone.normalization.reads``):

  * Reads are drawn without replacement, so each output pair has exactly the requested number of different reads.
    Before, repeated draws and an off-by-one error (the first read could never be picked) could give fewer reads,
    for example 987 when downsizing 1,000 reads to 1,000. **With the same seed, different reads are now selected.**
  * Asking for more reads than a file has raises an error; ``DownsizeDir`` leaves such pairs out with a warning.
  * ``DownsizePair`` no longer fails on uncompressed FASTQ without ``read_info``, and finds its input in ``in_dir``.
  * Files that are not FASTQ (a README, an image) are skipped with a warning instead of stopping the run, and files
    without a pair are reported.
  * Pairs are matched on the whole first read header, so a header that starts another sample's header no longer
    pairs files from different samples.
  * The reads-per-sample plot summary works with pandas 3.
* **Breaking change:** the command-line flags ``-p/--pca_plot`` and ``-k/--k_means`` are removed, with the module
  ``moonstone.analysis.clustering`` behind them. Both were broken: ``--pca_plot`` drew an empty figure and saved nothing,
  and ``--k_means`` wrote a metadata file with no rows. The flags now stop with an error saying they were removed.
* Classifiers rewritten as one module, ``moonstone.analysis.classification``, replacing ``classify.py`` and
  ``randomForest.py``. New options ``--model``, ``--kernel``, ``--explore_kernels``, ``--folds``, ``--seed`` and
  ``--group_column``. Support Vector Machine scaling now happens inside each training fold, folds are stratified and
  seeded, groups are class-weighted and scored with balanced accuracy, feature selection is nested inside the
  cross-validation, and feature importance is permutation importance on the test folds. ``--explore_kernels``
  compares kernels (exploratory) and ``--kernel auto`` chooses the kernel by nested cross-validation. Minimum
  scikit-learn raised to 1.5. (https://github.com/motleystate/moonstone/pull/123)
* New ``permanova`` and ``permanova_screen`` methods on all beta-diversity classes: Permutational Multivariate
  Analysis of Variance (PERMANOVA) with a dispersion check (PERMDISP), seeded permutations, and false discovery rate
  correction across variables. See :ref:`permanova`. (https://github.com/motleystate/moonstone/issues/129)
* New command-line flag ``--permanova [VARIABLE ...]`` with ``--taxa_level`` (``species`` by default, or ``all``),
  ``--no_read_threshold`` and ``--permutations``, and the matching function
  ``moonstone.analysis.permanova.run_permanova``. Reads a normalized count table with taxonomy columns, sums it to
  each level, filters taxa with the read-threshold knee, and corrects all tests of a run together. Existing flags
  are unchanged. (https://github.com/motleystate/moonstone/issues/131)
* PERMANOVA missing values: samples with a missing value in the tested variable are left out of that variable's
  test only (complete-case), for PERMANOVA and PERMDISP alike. Empty and whitespace-only values now count as missing
  instead of forming a group. New ``n_missing`` column in the results of ``permanova``, ``permanova_screen``,
  ``run_permanova`` and ``permanova.csv``; the count is logged as a warning and printed by ``--permanova``. Skipped
  variables are logged as warnings, and so is a variable with no value for more than half of the samples (missing
  values are never filled in).
* PERMANOVA defaults to 9,999 permutations instead of 999, in ``permanova``, ``permanova_screen``,
  ``run_permanova`` and ``--permutations``. New ``p_value_se`` and ``permdisp_p_value_se`` columns give the standard
  error of each permutation p-value, sqrt(p (1 - p) / N). A warning (logged, and in the ``warnings`` column) says when
  there are too few permutations for any test to pass the multiple-testing correction. New ``--threads`` flag
  (every available core by default) and ``n_jobs`` parameter run the tests in parallel with identical results.
  New dependencies listed explicitly (already installed with scikit-learn): ``joblib`` and ``threadpoolctl``.
* ``analyse_groups(pval_to_display="all")`` no longer crashes when ``groups`` is not given, or when a sample's
  group is blank in the metadata. The significance brackets are placed on the groups that are plotted, so a group
  that is only in the metadata no longer adds an empty box. With ``sort_groups=True``, the brackets are now placed
  after the groups are sorted; before, they could mark the wrong pairs of groups.
  (https://github.com/motleystate/moonstone/issues/101)
* ``analyse_groups(pval_to_display=...)`` now draws the significance brackets when ``group_col2`` is used (it was
  ignored with a message before), for the three choices ``"all"``, ``"same group_col values"`` and
  ``"same group_col or group_col2 values"``. Each bracket ends on its box inside the group cluster. P-values are
  still drawn on static images only. With ``mode="violin"``, the ``group_col2`` violins are now drawn side by side
  instead of on top of each other. (https://github.com/motleystate/moonstone/issues/101)


v1.1
====

v1.1.0
------

* Compatibility with Python 3.10. Upgrade of statsmodels requirement to version=1.14
* class ``ScatterTrendlines``: input counts dataframe; give the two species to plot against each other to ``plot_one_graph`` and the number of trendlines to generate (https://github.com/motleystate/moonstone/pull/104)
* 3 new functions in dict_operations: ``flatten_dict``, ``super_pop`` and ``filter_dict`` (https://github.com/motleystate/moonstone/pull/104)(https://github.com/motleystate/moonstone/pull/92)
* reorganization Chi2 and other statistical tests (https://github.com/motleystate/moonstone/pull/92)
* ``SeriesBinning``: to put series into bins (https://github.com/motleystate/moonstone/pull/92)

v1.0
====

v1.0.0
------

* upgrade to python >= 3.8. Version requirements have been updated for numpy, scikit-bio, scikit-learn, pandas, statsmodels, hdmedians, scipy and plotly
* ``mode`` argument for ``plot_most_abundant_taxa()`` and ``plot_most_prevalent_taxa()`` that allows to plot the data as boxplot or violin plot, as well as a bargraph.
* ``parsers``:
  * More formats accepted in parsers: ".xls", ".xlsx", ".odt", ".ods", ".odf", ".xlsb"
  * ``Metaphlan3Parser``: ``keep_NCBI_tax_col`` argument and possibility to compute dataframe with less taxonomical level than expected
* ``adapt_phylogenetic_tree_to_counts_df`` function in utils/phylogenetic_tree_editing.py

v0.7
====

v0.7.0
------

* https://github.com/motleystate/moonstone/pull/78
* ``Phylogenetic Diversity`` (https://github.com/motleystate/moonstone/pull/81):
  * new argument ``force_computation`` allows to force diversity computation even when species are missing in the Tree
  * only accepting trees in skbio.TreeNode
* ``PlotTaxonomyCounts`` (https://github.com/motleystate/moonstone/pull/83):
  * new arguments: ``color_df`` and ``sep_series``/``sep_how`` added to ``plot_sample_composition_most_abundant_taxa``
  * new ``plot_most_abundant_taxa`` method and renaming of the old one to ``plot_most_prevalant_taxa``
  * ``mode`` argument for ``plot_most_abundant_taxa()`` and ``plot_most_prevalent_taxa()`` to plot the data as a boxplot or a violin plot, as well as a bargraph (https://github.com/motleystate/moonstone/pull/89)
* ``plot_complex_graph`` in the ``MatrixBarGraph`` class to display some metadata below the main graph (https://github.com/motleystate/moonstone/pull/83)
* 2 new arguments, ``group_col2`` and ``groups2``, for the ``analyse_groups`` method from the ``DiversityBase`` class and its descendants (https://github.com/motleystate/moonstone/pull/86)
* ``orientation`` argument for the ``plot_one_graph`` method of the ``GroupBaseGraph`` class and its descendants (https://github.com/motleystate/moonstone/pull/86)

v0.6
====

v0.6.0
------

* https://github.com/motleystate/moonstone/pull/66
* https://github.com/motleystate/moonstone/pull/69
* Build in visualization for taxonomy counts https://github.com/motleystate/moonstone/pull/72
* Faith's PD and UniFrac https://github.com/motleystate/moonstone/pull/75

v0.5
====

v0.5.0
------

* Visualize metadata cat (https://github.com/motleystate/moonstone/pull/62)
* parser for metaphlan3 (https://github.com/motleystate/moonstone/pull/60)
* pvalue correction in analyse groups (https://github.com/motleystate/moonstone/pull/57)
* add ``HeatmapGraph`` (https://github.com/motleystate/moonstone/pull/56)
* Statistical tests for group comparison (https://github.com/motleystate/moonstone/pull/53)
* Mann - Whitney u test in alpha diversity (https://github.com/motleystate/moonstone/pull/50)
* Feature box plot alpha diversity (https://github.com/motleystate/moonstone/pull/46)
* More filtering methods : classes ``NaNPercentageFiltering`` and ``NumberOfDifferentValuesFiltering`` (https://github.com/motleystate/moonstone/pull/44)

v0.4
====

v0.4.0
------

* new features for ``AlphaDiversity`` class (https://github.com/motleystate/moonstone/pull/42):

  * ``visualize`` method now has two `mode`: ``histogram`` or ``violin``
  * ``visualize_groups`` method allows visualization with violin plot based on metadata grouped by a chosen column

* There is now definition of ``index_col`` for ``MetadataParser`` to define which column to use as the index for the metadata (https://github.com/motleystate/moonstone/pull/42)
* class ``GenesToTaxonomy`` to go from gene counts to taxonomy df (https://github.com/motleystate/moonstone/pull/36)
* add ``TaxonomyRandomSelection`` to allow random selection on multiindexed dataframe (https://github.com/motleystate/moonstone/pull/32)
* Classe ``AlphaDiversity`` and child classes ``ShannonIndex`` and ``SimpsonInverseIndex`` (https://github.com/motleystate/moonstone/pull/29)
* add ``RandomSelection`` to ``normalization`` module (https://github.com/motleystate/moonstone/pull/28)

v0.3
====

v0.3.0
------

* Refactor use of plotting_options : iteration through dict.keys() + exported in ``BaseGraph`` (Parent Class)
* ``plotting_options``, ``output_file`` and show dealt with in ``plot_one_graph`` (instead of at instantiation)
* Refactoring of ``BarGraph``
* Relocation of functions : ``check_list_of`` and ``add_x_to_plotting_options`` in `utils/plot.py`
* Relocation of method to bin series in ``utils/pandas/series.py``
* Transition to autosummary in available classes
* Listing of plot classes in available classes

v0.2
====

v0.2.0
------

* Add base for each module (with ``visualize()`` method and ``data_report`` property)
* Use the new base for ``MeanFiltering`` module.
* Refactoring of ``Filtering`` class into several classes:

  * ``NoCountsFiltering`` that filters on rows or columns with no counts at all
  * ``NamesFiltering`` that filters on a given list of row or columns names

    * Either keep the names
    * Or exclude them
  * ``TaxonomyNamesFiltering`` that filters on a given list of index names at a chosen level

    * Either keep the names
    * Or exclude them
* Add base module for plots
* Add class ``BaseGraph`` and child classes ``Histogram`` and ``BarGraph``

v0.1
====

v0.1.0
------

* First release of the work done on moonstone.
* Contains command line ``moonstone`` to run built-in analysis. See ``moonstone --help`` for more information.
* Starting modules to build your own analysis:

  * Parsers

    * Metadata
    * Counts

      * Qiime2
      * Kraken2
      * Picrust2
      * Metaphlan2
  * Normalization

    * GeometricMean
    * TotalCounts
    * StandardScalar
  * Filtering
