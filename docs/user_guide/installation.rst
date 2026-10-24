************
Installation
************

Pre-requisites
==============

moonstone requires Python 3.13.

.. Note::
    We highly recommend the use of a virtual environment such as `virtualenv`_, `pyenv`_ or `conda`_.

.. _virtualenv: https://virtualenv.pypa.io/en/latest/
.. _pyenv: https://github.com/pyenv/pyenv
.. _conda: http://docs.readthedocs.io/en/latest/conda.html

For example, with conda (or mamba):

.. code-block:: bash

    conda create -n moonstone python=3.13 -y
    conda activate moonstone


Installation procedure
======================

Pip
---

Install the latest release of moonstone from the Python Package Index (PyPI):

.. code-block:: bash

    pip install moonstone

To try changes that are not released yet, install from the Github repository instead:

.. code-block:: bash

    pip install "git+https://github.com/motleystate/moonstone.git"

.. Note::
    Releases from 1.5.0 onwards require Python 3.13. On older Python versions, pip installs moonstone 1.0.0.

Check the installation with:

.. code-block:: bash

    moonstone --help

Manually
--------

.. Note::
    This is particularly useful when you wish to install a version under development from
    any branches of the Github repository.

Clone the repository and install moonstone with the following commands:

.. code-block:: bash

    git clone https://github.com/motleystate/moonstone.git
    cd moonstone
    pip install .     # you can use -e option to make an editable install

To also install the development tools and run the test suite:

.. code-block:: bash

    pip install -e . -r requirements-dev.txt
    pytest

.. Note::
    Installing moonstone will install the following dependencies:
        - pandas (>=3.0.0)
        - numpy (>=1.26.4)
        - scipy (>=1.9.0, <1.16.0)
        - scikit-learn (>=1.3.1)
        - scikit-bio (>=0.6.0)
        - statsmodels (>=0.14.0)
        - matplotlib
        - plotly
        - python-slugify
        - pyaml
        - hdmedians
        - cython
        - filetype

Uninstallation procedure
=========================

You can remove moonstone with the following command:

.. code-block:: bash

    pip uninstall moonstone

.. Note::
    This will not uninstall dependencies. To do so you can make use of the pip-autoremove
    tool `pip-autoremove`_ or set up your environment with pipenv_ or poetry_...

.. _pip-autoremove: https://github.com/invl/pip-autoremove
.. _pipenv: https://github.com/pypa/pipenv
.. _poetry: https://python-poetry.org/docs/
