# Moonstone

![Python package](https://github.com/motleystate/moonstone/workflows/Python%20package/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/moonstone/badge/?version=latest)](https://moonstone.readthedocs.io/en/latest)
[![codecov](https://codecov.io/gh/motleystate/moonstone/branch/master/graph/badge.svg)](https://codecov.io/gh/motleystate/moonstone)

**M**etagenomics **O**bservation, **O**rganization, **N**ormalization, **S**tatistics, **T**esting, **O**rdination, **N**etworks and **E**xploration.

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

# Install the latest release of moonstone from PyPI
pip install moonstone

# Check the installation
moonstone --help
```

> To try changes that are not released yet, install from GitHub instead: `pip install "git+https://github.com/motleystate/moonstone.git"`. Releases from 1.5.0 onwards require Python 3.13; on older Python versions, pip installs moonstone 1.0.0.

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

## Graphical interface

`moonstone-ui` opens Moonstone in your browser: one page per tool, with a form built from the command-line options,
checks of your files before a run, example files, a tutorial and the documentation. Each page shows the `moonstone`
command it runs, so anything done in the interface can be repeated in a terminal or on a cluster.

```bash
pip install "moonstone[ui]"   # adds Streamlit
moonstone-ui                  # start it from the directory where results should go
```

The Tutorials page has worked examples for every tool and for the Python classes, run live on the example files
installed with Moonstone. Start there after a long break. More in the
[documentation](https://moonstone.readthedocs.io/en/latest/user_guide/interface.html).

On a cluster, `moonstone-ui --remote` serves the interface on a private socket, with no network port open, and you
reach it from your own computer through an SSH tunnel
([documentation](https://moonstone.readthedocs.io/en/latest/user_guide/remote.html)).

--------

## From sequencing reads to counts

`moonstone-reads` runs the [Moonbeam](https://github.com/motleystate/moonbeam2) pipeline from inside Moonstone:
adapter trimming, host-read removal, Kraken2 classification and Bracken, with the same results as Moonbeam on the
same inputs. It needs Kraken2, Bracken, Trimmomatic and KrakenTools; `environment-hpc.yml` creates an environment
with everything.

```bash
# 1. reads to Bracken-corrected raw counts, submitted to Slurm
moonstone-reads run READS_FOLDER RUN_FOLDER --db DATABASE --host-db HOST_DATABASE --slurm
moonstone-reads status RUN_FOLDER

# 2. the same counts, in the table layout the analyses read
moonstone-reads table RUN_FOLDER/output/mapping/all_samples.bracken.tsv counts_raw.csv

# 3. normalization, a separate step on any Moonstone count table
python -m moonstone.normalization.table counts_raw.csv counts_geomean.csv --method geometric_mean
```

The same three steps are tools of the interface. More in the
[documentation](https://moonstone.readthedocs.io/en/latest/user_guide/reads.html).

--------

## Development guide

To work on a new feature, related or not to an issue, we open a new branch from the development (currently `master`) branch to work on it.

> *If an issue is opened, we recommend to name your branch with the issue number at the begining, e.g. 5-my-new-feature*

Once done, work is reviewed through a pull request.
