************
Get started with moonstone
************



How to run moonstone?
============

More detailed documentation is available on the wiki_

.. _wiki: https://virtualenv.pypa.io/en/latest/(https://gitlab.pasteur.fr/metagenomics/data-analysis/-/wikis/Home).

moonstone can be run in a python script or by command-line in a terminal.


In a python script
########

.. code-block:: python

    import moonstone
    help(moonstone)


By command-line
########

.. code-block:: bash

    (moonstone) $ moonstone --help
    usage: moonstone [-h] [-f filtering] [-p] [-k clusters] [-s variable]
                    [-sr variable] [-sc variable] [-rf variable]
                    countfile metadata

    Microbiota Analysis Scripts for Machine LEarning

    positional arguments:
    countfile     Normalized count file input
    metadata      Clinical data input file

    optional arguments:
    -h, --help    show this help message and exit
    -f filtering  Minimum mean reads per variable: use a float >0
    -p            Generates PCA plot of data
    -k clusters   Runs K-Means clustering. Provide number of clusters
    -s variable   Run SVMachine Classification
    -sr variable  SVM Classifier with ROC Analysis
    -sc variable  SVM and output Classifier Components
    -rf variable  Random Forest Analysis, using supplied variable



What kind of data does moonstone handle?
============

moonstone can handles count data and taxonomy data files, and optionally their corresponding metadata files.
moonstone can also take `pandas dataframe`_ with the data as some modules' input (in python script exclusively).

.. _pandas dataframe : https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html