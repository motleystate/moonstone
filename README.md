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

Moonstone requires **Python 3.13**. If you only want to use `moonstone` and don't need to edit the source code:

```bash
# Create and activate a Python 3.13 environment
mamba create -n moonstone python=3.13 -y
conda activate moonstone

# Install the latest version of moonstone from GitHub
pip install "git+https://github.com/motleystate/moonstone.git"

# Check the installation
moonstone --help
```

> The version on the Python Package Index (PyPI, `pip install moonstone`) is 1.0.0, which predates Python 3.13 support and does not install on Python 3.13. Use the GitHub install above until a newer release is published on PyPI.

### 🛠️ Developer Installation

To work directly with the code and have your changes reflected automatically, follow these steps:

```bash
# Clone the repository
git clone https://github.com/motleystate/moonstone.git
cd moonstone

# Create a development environment
mamba create -n moonstone-dev python=3.13 -y
conda activate moonstone-dev

# Install Moonstone in editable mode, with the development tools
pip install -e . -r requirements-dev.txt

# Run the test suite
pytest
```

Without Mamba or Conda, a standard virtual environment works too: `python3.13 -m venv .venv && source .venv/bin/activate`, then run the same `pip install` command.

--------

## Quickstart with moonstone on the command-line

More detailed documentation is available on the [documentation](https://moonstone.readthedocs.io/en/latest/?badge=latest).

Moonstone is directly callable from your terminal for available built-in analysis scripts:

```bash
(moonstone) $ moonstone --help
usage: moonstone countfile, metadata file, output directory [options]

Microbiota Analysis Scripts for Machine LEarning

positional arguments:
  countfile             Normalized count file input
  metadata              Clinical data input file
  outdir                Output files directory

options:
  -h, --help            show this help message and exit
  -os, --skip_prompt    Suppress output directory exists prompt.
  -f, --filtering min_mean_read
                        Minimum mean reads per variable: use a float >0
  -c, --classify variable
                        Cross-validated classification of a metadata variable
                        from the counts
  --model {svm,random_forest} [{svm,random_forest} ...]
                        Classifier(s) for --classify: Support Vector Machine
                        and/or random forest (default: svm)
  --kernel {linear,poly,rbf,sigmoid,auto}
                        Support Vector Machine kernel; auto picks one by
                        nested cross-validation (default: linear)
  --explore_kernels     With --classify, compare Support Vector Machine
                        kernels (exploratory scores)
  --folds FOLDS         Number of cross-validation folds (default: 5)
  --seed SEED           Random seed for folds and models (default: 0)
  --group_column column
                        Metadata column naming the subject of each sample;
                        keeps a subject in one fold
  --roc                 With --classify, plot Receiver Operating
                        Characteristic curves
  --feature_selection   With --classify, select features by nested recursive
                        feature elimination
```

For example, `moonstone counts.csv metadata.csv output --classify disease --model svm random_forest --roc`
runs a 5-fold stratified cross-validation of both classifiers. Scores are balanced accuracy and the area under the
Receiver Operating Characteristic (ROC) curve, measured on samples the model did not see during training.

--------

## Development guide

To work on a new feature, related or not to an issue, we open a new branch from the development (currently `master`) branch to work on it.

> *If an issue is opened, we recommend to name your branch with the issue number at the begining, e.g. 5-my-new-feature*

Once done, work is reviewed through a pull request.
