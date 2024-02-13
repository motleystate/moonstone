.. _changelog:

**********
Changelogs
**********

Summary of developments of moonstone library.


v1.1
====

v1.1.0
------

* Compatibility with Python 3.10. Upgrade of statsmodels requirement to version=1.14
* class `ScatterTrendlines`: input counts dataframe; give the two species to plot against each other to `plot_one_graph` and the number of trendlines to generate (https://github.com/motleystate/moonstone/pull/104)
* 3 new functions in dict_operations: `flatten_dict`, `super_pop` and `filter_dict` (https://github.com/motleystate/moonstone/pull/104)(https://github.com/motleystate/moonstone/pull/92)
* reorganization Chi2 and other statistical tests (https://github.com/motleystate/moonstone/pull/92)
* `SeriesBinning`: to put series into bins (https://github.com/motleystate/moonstone/pull/92)

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
* class ``GenesToTaxonomy``to go from gene counts to taxonomy df (https://github.com/motleystate/moonstone/pull/36)
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
