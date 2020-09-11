.. _changelog:

**********
Changelogs
**********

Summary of developments of moonstone library.

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
