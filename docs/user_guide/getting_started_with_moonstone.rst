******************************
Getting started with moonstone
******************************

How to run moonstone?
=====================

Moonstone can be used in two different ways:

By command-line
###############

.. Note::
    This feature is going to change in the future and is still under construction

A list of built-in scripts for analysis is available and can be used by command line.
They are built from the different moonstone modules.

Please refer to the help section of the command-line tool for more information:

.. code-block:: bash

    $ moonstone --help

In Python
#########

Moonstone provides a list of modules to perform the different steps of your analysis and are described below.

.. Note::
    This part is under construction and will be completed while updating each modules

Parsers
"""""""

A crucial step of every analysis is to import the data. Here the aim is to be able to import data from many
different source and format them in a common format using pandas dataframes.

What kind of data does moonstone handle?
''''''''''''''''''''''''''''''''''''''''

Moonstone parsers can handle count data for items (such as genes, functions...) or taxonomy (setting up a more
complexe dataframe structure with multi-indexing). There is also a parser for metadata files that are crucial for analysis.

Cleaning and Transforms
'''''''''''''''''''''''

Filtering
"""""""""

Normalization
"""""""""""""

Analysis
""""""""

Visualization
"""""""""""""
