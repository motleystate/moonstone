************
Installation
************

Pre-requisites
==============

moonstone currently works with version of Python 3.6 and 3.7.

.. Note::
    We highly recommend the use of a virtual environment such as `virtualenv`_, `pyenv`_ or `conda`_.

.. _virtualenv: https://virtualenv.pypa.io/en/latest/
.. _pyenv: https://github.com/pyenv/pyenv
.. _conda: http://docs.readthedocs.io/en/latest/conda.html


Installation procedure
======================

Pip
---

You can use pip to install moonstone of the latest stable version published on pypi:

.. code-block:: bash

    pip install moonstone

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

.. Note::
    Installing moonstone will install the following dependencies :
        - pandas (==1.0.1)
        - scikit-learn (==0.21.3)
        - matplotlib
        - plotly
        - statsmodels
        - python-slugify
        - pyaml
        - numpy (==1.18.1)

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
