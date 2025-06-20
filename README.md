# Moonstone

![Python package](https://github.com/motleystate/moonstone/workflows/Python%20package/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/moonstone/badge/?version=latest)](https://moonstone.readthedocs.io/en/latest)
[![codecov](https://codecov.io/gh/motleystate/moonstone/branch/master/graph/badge.svg)](https://codecov.io/gh/motleystate/moonstone)

Moonstone aims to provide a way of performing analysis of metagenomics counts from raw data to statistical analysis and visualization of the results.

Thus, in moonstone you will find:

* parsers for common file types for metagenomics counts
* modules for cleaning and filtering your data
* normalization modules
* analysis modules
* plot modules

The main idea is to keep track of every steps applied to a raw data to be able to easily share the analysis and reproduce it.

Please check our [documentation](https://moonstone.readthedocs.io/en/latest/?badge=latest) for more information.

## Install moonstone

You can install `moonstone` either as a **user** (to use the tool as-is) or as a **developer** (to work directly on the code). We recommend using [Mamba](https://github.com/mamba-org/mamba) (a faster drop-in replacement for Conda) to manage environments.

### 🧪 User Installation

If you only want to use `moonstone` and don't need to edit the source code:

```bash
# Create and activate a Python 3.10 environment
mamba create -n moonstone python=3.10 -y
conda activate moonstone

# Install the latest published version of moonstone from PyPI
pip install moonstone
```

### 🛠️ Developer Installation

To work directly with the code and have your changes reflected automatically, follow these steps:

```bash
# Clone the repository**  
git clone git@github.com:motleystate/moonstone.git
cd moonstone

# Create a development environment
mamba create -n moonstone-dev python=3.10 -y
conda activate moonstone-dev

# Install Moonstone in editable mode
pip install -e .
```

You can also install dependencies required for development : `pip install -r requirements-dev.txt`.

--------

## Quickstart with moonstone on the command-line

More detailed documentation is available on the [documentation](https://moonstone.readthedocs.io/en/latest/?badge=latest).

Moonstone is directly callable from your terminal for available built-in analysis scripts:

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

Once done, work is reviewed through a pull request.
