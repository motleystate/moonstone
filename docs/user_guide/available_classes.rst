.. _av_classes:

*****************
Available classes
*****************

Here are listed the different classes available by module.

.. Note::
    We recommend to check the API documentation for more detailed documentation about
    each class

.. _av_parsers:

Parsers
=======

How it works?
"""""""""""""

.. code-block:: python

    from moonstone.parsers import YourFavouriteParser

    parser = YourFavouriteParser("/path/to/the/file")
    df = parser.dataframe

List
""""

Classic and simple parsers:

.. currentmodule:: moonstone.parsers.metadata
.. autosummary::
   :nosignatures:

    MetadataParser
    YAMLBasedMetadataParser

Counts parsers (``from moonstone.parsers.counts``):

.. currentmodule:: moonstone.parsers.counts
.. autosummary::
   :nosignatures:

    genes.GeneCountsParser
    picrust2.Picrust2PathwaysParser

Taxonomy counts parsers (``from moonstone.parsers.counts.taxonomy``):

.. currentmodule:: moonstone.parsers.counts.taxonomy
.. autosummary::
   :nosignatures:

    kraken2.SunbeamKraken2Parser
    metaphlan.Metaphlan2Parser
    metaphlan.Metaphlan3Parser
    qiime.Qiime2Parser

.. Note::
    More details on API documentation for :ref:`api_parsers`.

.. _av_filtering:

Filtering
=========

How it works?
"""""""""""""

.. code-block:: python

    from moonstone.filtering import YourFavouriteFiltering

    filtering_instance = YourFavouriteFiltering(your_df)
    df = filtering_instance.filtered_df

List
""""

.. currentmodule:: moonstone.filtering
.. autosummary::
   :nosignatures:

    basics_filtering.NoCountsFiltering
    basics_filtering.NamesFiltering
    basics_filtering.NaNPercentageFiltering
    basics_filtering.NumberOfDifferentValuesFiltering
    taxonomy_filtering.TaxonomyNamesFiltering
    taxonomy_filtering.TaxonomyMeanFiltering
    mean_filtering.MeanFiltering

.. Note::
    More details on API documentation for :ref:`api_filtering`.

.. _av_normalization:

Normalization
=============

How it works?
"""""""""""""

.. code-block:: python

    from moonstone.normalization import YourFavouriteNormalization

    normalization = YourFavouriteNormalization(your_df)
    df = normalization.normalized_df

List
""""

.. currentmodule:: moonstone.normalization
.. autosummary::
   :nosignatures:

    counts.geometric_mean.GeometricMeanNormalization
    counts.total_counts.TotalCountsNormalization
    counts.random_selection.RandomSelection
    counts.random_selection.TaxonomyRandomSelection
    processed.scaling_normalization

.. _av_plot:

Plot
====

How it works?
"""""""""""""

.. code-block:: python

    from moonstone.plot import PlotCountsStats

    # instantiation
    plot_instance = PlotCountsStats(df)

    # call one (or more) plotting method(s)
    plot_instance.your_favorite_plot()
    plot_instance.another_of_your_favorite_plot()

List
""""

.. currentmodule:: moonstone.plot
.. autosummary::
   :nosignatures:

   counts.PlotCountsStats
   counts.PlotTaxonomyCounts
   metadata.PlotMetadataStats

.. Note::
    More details on API documentation for :ref:`api_plot`.

Analysis
========

How it works?
"""""""""""""

.. code-block:: python

    from moonstone.analysis import DifferentialAnalysis

    # instantiation
    analysis_instance = DifferentialAnalysis(df, metadata_df)
    analysis_instance.differential_analysis_by_feature(features, type_of_features, test_to_use, correction_method_used)

.. Note::
    The way analysis instances work will change using a defined common way of performing analysis.

List
""""

.. currentmodule:: moonstone.analysis
.. autosummary::
   :nosignatures:

   differential_analysis.DifferentialAnalysis
   diversity.alpha.ShannonIndex
   diversity.alpha.SimpsonInverseIndex
   diversity.alpha.Chao1Index
   diversity.alpha.FaithsPhylogeneticDiversity
   diversity.beta.BrayCurtis
   diversity.beta.WeightedUniFrac
   diversity.beta.UnweightedUniFrac

.. Note::
    More details on API documentation for :ref:`api_analysis`.