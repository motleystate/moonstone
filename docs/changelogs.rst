.. _changelog:

**********
Changelogs
**********

Summary of developments of moonstone library.

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
