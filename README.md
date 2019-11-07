# Metagenomics Data Analysis

[![pipeline status](https://gitlab.pasteur.fr/metagenomics/data-analysis/badges/master/pipeline.svg)](https://gitlab.pasteur.fr/metagenomics/data-analysis/commits/master)
[![coverage report](https://gitlab.pasteur.fr/metagenomics/data-analysis/badges/master/coverage.svg)](https://gitlab.pasteur.fr/metagenomics/data-analysis/commits/master)

This repository aims to gather different scripts and items to help performing data analysis from metagenomics data.

## Install

### Set up virtualenv

Set up a Python 3 virtual env, for instance:

```bash
python3 -m virtualenv metautils
source metautils/bin/activate
(moonstone) $
```

Then simply clone and install the library:

```bash
(moonstone) $ git clone https://gitlab.pasteur.fr/metagenomics/data-analysis.git
(moonstone) $ cd data-analysis
(moonstone) $ pip install .
```

You can also install dependencies required for devlopment : `pip install -r requirements-dev.txt`.
