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

+----------------------------+--------------------------------------------------------------------------------+
| Class name                 | Description                                                                    |
+============================+================================================================================+
| ``MetadataParser``         | Parse metadata file and allows to apply transformations on them (cleaning...). |
+----------------------------+--------------------------------------------------------------------------------+

Counts parsers (``from moonstone.parsers.counts``):

+----------------------------+----------------------------------------------------------------+
| Class name                 | Description                                                    |
+============================+================================================================+
| ``Picrust2PathwaysParser`` | Predicted sample pathway abundances output file from Picrust2_ |
+----------------------------+----------------------------------------------------------------+

.. _Picrust2: https://github.com/picrust/picrust2/wiki

Taxonomy counts parsers (``from moonstone.parsers.counts.taxonomy``):

+--------------------------+---------------------------------------------------------+
| Class name               | Description                                             |
+==========================+=========================================================+
| ``SunbeamKraken2Parser`` | output from Kraken2_ merge table from Sunbeam_ pipeline |
+--------------------------+---------------------------------------------------------+
| ``Metaphlan2Parser``     | output from Metaphlan2_ merged table                    |
+--------------------------+---------------------------------------------------------+
| ``Qiime2Parser``         | parse output csv data obtained by Qiime2_               |
+--------------------------+---------------------------------------------------------+

.. _Sunbeam: https://github.com/sunbeam-labs/sunbeam
.. _Kraken2: https://ccb.jhu.edu/software/kraken2/
.. _Metaphlan2: https://github.com/biobakery/MetaPhlAn
.. _Qiime2: https://qiime2.org/

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

+----------------------------+-------------------------------------------------------------------------------+
| Class name                 | Description                                                                   |
+============================+===============================================================================+
| ``NoCountsFiltering``      | Remove rows (default) or columns with no counts.                              |
+----------------------------+-------------------------------------------------------------------------------+
| ``NamesFiltering``         | Filtering based on row (default) or column names.                             |
+----------------------------+-------------------------------------------------------------------------------+
| ``TaxonomyNamesFiltering`` | Filtering a Taxonomy multiindexed dataframe on index names at a chosen level. |
+----------------------------+-------------------------------------------------------------------------------+
| ``MeanFiltering``          | Remove items below a given mean threshold.                                    |
+----------------------------+-------------------------------------------------------------------------------+

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

+--------------------------------+------------------------------------------------------+
| Class name                     | Description                                          |
+================================+======================================================+
| ``GeometricMeanNormalization`` | normalization based on the one performed by DeSeq2_. |
+--------------------------------+------------------------------------------------------+
| ``TotalCountsNormalization``   | normalization based on total counts.                 |
+--------------------------------+------------------------------------------------------+

.. _DeSeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html