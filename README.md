# Moonstone

![Python package](https://github.com/motleystate/moonstone/workflows/Python%20package/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/moonstone/badge/?version=latest)](https://moonstone.readthedocs.io/en/latest/?badge=latest)

This repository aims to gather different scripts and items to help performing data analysis from metagenomics data.

It also contains work for the `moonstone` library that will be moved to github in the future. This library aims to provide a way of performing analysis of metagenomics counts from raw data to statistical analysis.

Thus, in moonstone you will find:

* parsers for common file types for metagenomics counts
* modules for cleaning and filtering
* normalization modules
* analysis modules
* visualization modules

The idea is to keep track of every steps applied to a raw data to be able to easily share the analysis and reproduce it.

## Install moonstone

Set up a Python 3 virtual env, for instance:

```bash
python3 -m virtualenv moonstone
source moonstone/bin/activate
(moonstone) $
```

Then simply clone and install the library:

```bash
(moonstone) $ git clone https://gitlab.pasteur.fr/metagenomics/data-analysis.git
(moonstone) $ cd data-analysis
(moonstone) $ pip install .  # you can use -e option to make an editable install
```

You can also install dependencies required for development : `pip install -r requirements-dev.txt`.

--------

## Quickstart with moonstone

More detailed documentation is available on the [wiki](https://gitlab.pasteur.fr/metagenomics/data-analysis/-/wikis/Home).

Moonstone is directly callable from your terminal:

```bash
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
```

--------

## Development guide

To work on a new feature, related or not to an issue, we open a new branch from the development (currently `master`) branch to work on it.

> *If an issue is opened, we recommend to name your branch with the issue number at the begining, e.g. 5-my-new-feature*

Once done, work is reviewed through a merge request to keep developers up to date on the development of moonstone.
