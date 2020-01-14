# Metagenomics Data Analysis

[![pipeline status](https://gitlab.pasteur.fr/metagenomics/data-analysis/badges/master/pipeline.svg)](https://gitlab.pasteur.fr/metagenomics/data-analysis/commits/master)
[![coverage report](https://gitlab.pasteur.fr/metagenomics/data-analysis/badges/master/coverage.svg)](https://gitlab.pasteur.fr/metagenomics/data-analysis/commits/master)

This repository aims to gather different scripts and items to help performing data analysis from metagenomics data.

## Install

### Set up virtualenv

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

## Run moonstone

Moonstone is directly callable from your terminal:

```
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

## Parsing and cleaning based on YAML config file

It is now possible to perform parsing and cleaning of an input table thanks to a config YAML file.

For the moment, only documentation available is the test available at : `tests/parsers/test_metadata:TestYAMLBasedMetadataParser.test_parse_end_to_end`:

### Input table

```csv
samples,age
s1,29
S2 ,48
s3,36
s4 ,25
```

### Config File

```yaml
parsing:
  - col_name: samples
    operations:
    - name: to_slug
    - name: rename
      options:
        new_name: 'sample'
  - col_name: age
    dtype: object
```

### Usage

```python
metadata_file_dirty = os.path.join(os.path.dirname(__file__), "data/metadata/dirty_metadata.tsv")
config_file = os.path.join(os.path.dirname(__file__), "data/metadata/config.yaml")
parser = YAMLBasedMetadataParser(metadata_file_dirty, config_file, sep=",")
expected_df = pd.DataFrame(
    {
        'sample': ['s1', 's2', 's3', 's4'],
        'age': ['29', '48', '36', '25'],
    }
)
pd.testing.assert_frame_equal(parser.metadata_parser.dataframe, expected_df)
```

Better documentation coming soon...

--------

## Development guide

To work a a new feature, related or not to an issue, we open a new branch from the development (currently `master`) branch to work on it.

> *If an issue is opened, we recommend to name your branch with the issue number at the begining, e.g. 5-my-new-feature*

Once done, work is reviewed through a merge request to keep developers up to date on the development of moonstone.
