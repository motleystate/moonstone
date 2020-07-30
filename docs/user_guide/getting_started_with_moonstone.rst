******************************
Getting started with moonstone
******************************

Moonstone can be used in two different ways:

By command-line
===============

.. Warning::
    This feature is going to change in the future and is still under construction

A list of built-in scripts for analysis is available and can be used by command line.
They are built from the different moonstone modules described below.

Please refer to the help section of the command-line tool for more information:

.. code-block:: bash

    $ moonstone --help

In Python
=========

Moonstone provides a list of modules to perform the different steps of your analysis. In this part,
you will find a detailed description of each module.

.. Warning::
    This part is under construction and will be completed while updating each modules

Parsers
"""""""

A crucial step of every analysis is to import the data. Here the aim is to be able to import data from many
different sources and format them in a homogeneous configuration using pandas dataframes.

.. Note::
    The parsers module can be described as the ETL_ procedure of moonstone.

.. _ETL: https://en.wikipedia.org/wiki/Extract,_transform,_load

.. image:: /img/countparser.png
  :alt: Countparser example

What format Moonstone can handle?
'''''''''''''''''''''''''''''''''

Detailed information about the different parsers can be found in the API documentation about :ref:`api_parsers`.
In brief, parsers are available for:

.. Note::
    Table example represent moonstone data format once the data is parsed.

- Simple count matrices (columns being samples, rows being items)

+--------+----------+----------+----------+
| index  | sample-1 | sample-2 | sample-3 |
+========+==========+==========+==========+
| gene-1 | 22       | 0        | 35       |
+--------+----------+----------+----------+
| gene-2 | 35       | 29       | 56       |
+--------+----------+----------+----------+

- Taxonomy counts matrices

+------------------------------------+---------------------+
|                indexes             |       columns       |
+----------+------------+------------+----------+----------+
| kingdom  | phylum     | class      | sample-1 | sample-2 |
+==========+============+============+==========+==========+
| Bacteria | Firmicutes | Clostridia | 35       | 56       |
+----------+------------+------------+----------+----------+

.. Note::
    Multi-indexes are used to represent taxonomy counts in order to facilitate grouping at the chosen level.

- Metadata (table with each row containing information about one sample)

+----------+-----+-----+--------+
| index    | age | sex | smoker |
+==========+=====+=====+========+
| sample-1 | 22  | f   | y      |
+----------+-----+-----+--------+
| sample-2 | 35  | m   | n      |
+----------+-----+-----+--------+

Cleaning and Transformers
'''''''''''''''''''''''''

It is not rare to have the necessity to clean up a bit the data you wish to analyse before analysis. This is often
the case for metadata that are often manually obtained through a spreadsheet without any specific constraints
leading to mistakes and inconsistency in the matrix.

As part of the **Parsers** module, some cleaning and transforming operations can be applied.
Here are some examples of operations you might need to do:

- Transform a ``string`` to snakecase_ format (particularly useful to make sure sample IDs are identical between counts and metadata)
- Remove trailing spaces from a ``string``.
- Group under a common name different values representing the same thing (e.g. M, m and male for sex).

.. _snakecase: https://en.wikipedia.org/wiki/Snake_case

.. Note::
    More details about available operations can be found in the API documentation about :ref:`api_parsers`.

Filtering
"""""""""

Normalization
"""""""""""""

Analysis
""""""""

Visualization
"""""""""""""
