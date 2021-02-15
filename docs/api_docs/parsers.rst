.. _api_parsers:

*******
Parsers
*******

How it works?
=============

.. code-block:: python

    from moonstone.parsers import YourFavouriteParser

    parser = YourFavouriteParser("/path/to/the/file")
    df = parser.dataframe

---------

Counts
======

Simple Counts
"""""""""""""

.. currentmodule:: moonstone.parsers.counts

.. autosummary::
    :toctree: stubs

    GeneCountsParser
    Picrust2PathwaysParser

Taxonomy Counts
"""""""""""""""

.. currentmodule:: moonstone.parsers.counts.taxonomy

.. autosummary::
    :toctree: stubs

    SunbeamKraken2Parser
    Metaphlan2Parser
    Metaphlan3Parser
    Qiime2Parser

Metadata
========

.. currentmodule:: moonstone.parsers.metadata

.. autosummary::
    :toctree: stubs

    MetadataParser
    YAMLBasedMetadataParser
