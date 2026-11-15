.. _interface:

*******************
Graphical interface
*******************

``moonstone-ui`` opens Moonstone in the browser: every tool on one page, with its documentation, a tutorial,
example files, checks of your own files, and the command it runs. It is meant for coming back to Moonstone after
months away without having to remember a flag.

The interface does not replace the command line: it **builds a** ``moonstone`` **command, shows it, and runs it**.
A result obtained in the interface is the result of that command, and the command can be copied into a terminal or
a cluster job.

Install and start
=================

The interface needs `Streamlit <https://streamlit.io>`_, which is not installed by default:

.. code-block:: bash

    pip install "moonstone[ui]"

Start it from the directory where the results should go:

.. code-block:: bash

    moonstone-ui

It runs on your machine only (address ``127.0.0.1``) and follows the light or dark setting of your system. Extra
arguments are passed to ``streamlit run``, for example ``moonstone-ui --server.port 8502``.

The pages
=========

Home
    One card per tool: the question it answers and the files it needs.

Tutorials
    Worked examples for every tool and for the classes of the Python API (parsers, filtering, normalization,
    transforms, diversity, plots). Each piece of code is **run on the example files when the tutorial is opened**,
    and its real result is shown under it, so a tutorial cannot describe something the code no longer does.
    Start with *Your first analysis, start to finish*.

Check my files
    Checks a file against the format Moonstone expects (UTF-8, comma or tab separators, decimal dots, ISO dates),
    and, if you name the tool, that this tool can read it and that its samples match the metadata. A file saved
    by Excel in a French locale can be repaired: the repair is written to a **new file**, with the list of every
    value that changed. Your file is never modified.

Read threshold, Classification, PERMANOVA, Correlation, Association
    One page per tool, with four tabs:

    * **Run**: the form. Each option shows its command-line flag and the help text of ``moonstone --help``.
      Metadata variables are chosen from the columns of your metadata file. *Use the example files* fills the form
      with files and options that work. The files and options are checked as you go, and *Run* stays disabled
      while a check reports an error. The command is displayed above the button.
    * **Tutorial**: what the tool is for, how to run it, and how to read each column of its output.
    * **Inputs and examples**: the rules for the input files, a preview of a correct file, and what a run writes.
    * **Documentation**: the options with their defaults, and the documentation of the functions behind the tool.

Previous runs
    Every run made from the interface, with its status, command, log and results.

Python API
    Every public class and function, with its signature and documentation, searchable by name or by word.

Where the results go
====================

Each run is written to a new directory under the *run directory* shown in the sidebar (``moonstone_runs`` in the
directory where ``moonstone-ui`` was started, by default). A directory is never reused and nothing is overwritten.
Besides the files of the tool, a run directory holds:

* ``command.txt``: the command, the version of Moonstone and the start time;
* ``run.log``: everything the command printed;
* ``exit_code.txt``: 0 when the run finished, another number when it failed (absent if it was interrupted).

A run occupies its page until it finishes. For a long analysis, copy the command to a terminal or a cluster job
instead.

Example files
=============

The examples are installed with Moonstone, in ``moonstone/examples`` (all values are simulated):

.. code-block:: python

    from moonstone import examples

    examples.NAMES                           # the files available
    examples.path('species_counts.tsv')      # where one of them is

For developers: adding a tool
=============================

A tool is declared once, in ``moonstone/ui/catalog.py``: how it is run, which options of its parser the form
shows, how its count table is read, and an example. Its page, form, checks and tests follow from that entry.
Option names, defaults and help texts are never written in the interface: they are read from the tool's
``argparse`` parser (``moonstone.main.build_parser``).

Two tests keep the interface in step with the code: ``tests/ui/test_catalog.py`` fails when a flag is added to
the command line without being placed on a form, and ``tests/ui/test_tutorials.py`` runs every tutorial.
