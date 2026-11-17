.. _reads:

*******************************
From sequencing reads to counts
*******************************

``moonstone-reads`` takes a folder of shotgun sequencing reads and gives the number of reads of each taxon in
each sample, ready for the Moonstone analyses. It is the `Moonbeam <https://github.com/motleystate/moonbeam2>`_
pipeline moved into Moonstone: for the same reads, databases and settings it gives the same files.

There are three steps, each usable on its own, from a terminal or from the interface (:ref:`interface`):

1. **Reads to counts** (``moonstone-reads run``): the pipeline. Its result is the Bracken-corrected raw counts.
2. **Counts to Moonstone table** (``moonstone-reads table``): the same counts in Moonstone's table layout.
3. **Normalize** (``python -m moonstone.normalization.table``): a separate, optional step, on any Moonstone table.

Installation
============

The pipeline needs Snakemake and four programs that are not Python packages: Kraken2, Bracken, Trimmomatic and
KrakenTools. The file ``environment-hpc.yml`` of the repository creates one environment with all of them, Moonstone
and its interface:

.. code-block:: bash

    mamba env create -f environment-hpc.yml
    conda activate moonstone-hpc

In an environment that already has the four programs, ``pip install "moonstone[reads]"`` adds the rest. The table
and normalization steps need no extra.

The pipeline is meant for a cluster with Slurm: the databases and the reads are large. It also runs on one
machine, without ``--slurm``.

What the pipeline does
======================

====  =====================  ==========================================================================
Step  Program                What happens
====  =====================  ==========================================================================
1     Trimmomatic            Sequencing adapters are clipped from the reads.
2     Kraken2 (host)         Reads are classified against a database of host genomes.
3     KrakenTools            Reads classified as a chosen host are removed.
4     Kraken2                The remaining reads are classified against your database.
5     Bracken                Read counts are re-estimated at species level.
6     kraken-biom            One table for all samples, and a run report.
====  =====================  ==========================================================================

Steps 1, 2-3 and 5 can each be switched off (``--no-trimming``, ``--remove-host`` with no name, ``--no-bracken``).

Input
=====

Reads
    One folder; its sub-folders are searched. Files end with ``.fastq.gz``, ``.fq.gz``, ``.fastq`` or ``.fq``.
    Paired-end files are recognised by the endings ``_1``/``_2``, ``_R1``/``_R2`` or ``_R1_001``/``_R2_001`` before
    the extension, and the sample name is what comes before that ending. Single-end files (``--single-end``) have
    no such ending. Paired and single-end reads cannot be mixed in one run. ``--match`` keeps only the files whose
    name matches a regular expression. Your files are never renamed or moved: the run folder gets links to them.

Databases
    A Kraken2 database to classify the reads (``--db``), holding the Bracken k-mer file of your read length
    (``database150mers.kmer_distrib`` for 150 bases, the default of ``--read-len``), and, to remove host reads, a
    Kraken2 database of the host genomes (``--host-db``).

Problems are reported before anything runs or is written: a mate without its pair, a sample found twice, a folder
that is not a Kraken2 database, a missing k-mer file, a missing program.

Running
=======

.. code-block:: bash

    # on a cluster: submit and return at once
    moonstone-reads run READS_FOLDER RUN_FOLDER --db DATABASE --host-db HOST_DATABASE --slurm
    moonstone-reads status RUN_FOLDER
    moonstone-reads cancel RUN_FOLDER

    # on one machine, on 8 cores
    moonstone-reads run READS_FOLDER RUN_FOLDER --db DATABASE --host-db HOST_DATABASE --cores 8

    # write the run folder and list the steps, without running them
    moonstone-reads run READS_FOLDER RUN_FOLDER --db DATABASE --host-db HOST_DATABASE --dry-run

With ``--slurm``, one small job runs the workflow and submits one job per step and sample. The run does not depend
on your terminal or on the interface: close them and ask for the status later.

Options
"""""""

================================  ====================================================================================
Option                            Meaning
================================  ====================================================================================
``--remove-host [NAME ...]``      Hosts whose reads are removed: ``human`` (taxid 9606), ``mouse`` (10090); both by
                                  default. Given with no name, host removal is skipped and no host database is needed.
``--keep-trimmed-reads``          Keep the trimmed reads. This allows mapping with another database without trimming
                                  again. By default they are deleted once no step needs them.
``--keep-kraken-output``          Keep the per-read Kraken2 files, which are large. The per-sample reports are always
                                  kept.
``--single-end``, ``--match``     See *Input*.
``--no-trimming``                 Skip adapter trimming.
``--no-bracken``                  Skip Bracken: Kraken2 counts only.
``--read-len``                    Read length for Bracken (default 150).
``--bracken-threshold``           Bracken: minimum number of reads for a species to be re-estimated (default 5).
``--min-reads``                   Report: warn when fewer reads reach the classification (default 100000).
``--min-classified-pct``          Report: warn when less than this percentage is classified (default 20).
``--slurm``, ``--cores``          Where the steps run.
``--dry-run``                     List the steps without running them.
================================  ====================================================================================

Settings of your cluster
""""""""""""""""""""""""

Values that belong to one cluster are read from a site file, ``~/.config/moonstone/site.yaml`` (or the file named by
the environment variable ``MOONSTONE_SITE``). Every key is optional:

.. code-block:: yaml

    db: /shared/databases/kraken2/my_database        # default of --db
    host_db: /shared/databases/kraken2/hosts          # default of --host-db
    submit: sbatch                                    # the command that submits a job
    head_options: -p my_partition -c 1 --mem 4G       # its options for the job that runs the workflow
    profile: /path/to/slurm-profile                   # a Snakemake Slurm profile folder
    adapter_dir: /path/to/trimmomatic/adapters        # only if Trimmomatic was not installed with conda

The profile sets the partition and the resources of the step jobs. Copy the one installed with Moonstone
(``moonstone/reads/workflow/profiles/slurm``), add ``slurm_partition: my_partition`` under ``default-resources``,
and name the copy in the site file.

What a run writes
=================

In the run folder:

* ``output/pipeline_report.md``: versions, settings, and for each sample the reads, the percentage of host reads and
  of classified reads, with warnings for samples with few reads or few classified reads. Read it first.
* ``output/mapping/all_samples.bracken.tsv``: the Bracken-corrected raw counts of all samples. This is the result.
* ``output/mapping/all_samples.kraken.tsv``: the Kraken2 counts, before the Bracken correction.
* Per-sample Kraken2 and Bracken reports under ``output/mapping/``.
* ``samples.tsv`` and ``config.yaml``: the read files and the settings used.

From the pipeline's table to a Moonstone table
==============================================

The pipeline's table has one lineage per row in a single column, and names each sample column after a report file
(``S1.brackenCorrected_taxa``). The analyses read a table with the taxonomy in columns ``kingdom`` ... ``species`` and
sample names that match the metadata:

.. code-block:: bash

    moonstone-reads table RUN_FOLDER/output/mapping/all_samples.bracken.tsv counts_raw.csv

No count is changed: ``counts_raw.csv`` holds the raw, Bracken-corrected counts. It is read by PERMANOVA,
correlation, association and read threshold. A note beside it (``counts_raw.csv.methods.txt``) says what was done.

.. code-block:: python

    from moonstone.reads.table import read_moonbeam_table, to_moonstone_table

    counts = read_moonbeam_table("all_samples.bracken.tsv")          # as a dataframe
    to_moonstone_table("all_samples.bracken.tsv", "counts_raw.csv")   # as a file

Normalization
=============

Normalization is a step of its own, on any Moonstone count table, so the raw counts stay available and the method
can change without running the pipeline again. A method must be chosen; there is no default:

.. code-block:: bash

    python -m moonstone.normalization.table counts_raw.csv counts_geomean.csv --method geometric_mean
    python -m moonstone.normalization.table counts_raw.csv counts_relative.csv --method relative_abundance

``geometric_mean``
    Median of ratios, as in DESeq2: each sample is divided by one scaling factor. Values stay on a count scale, so
    the table suits every Moonstone analysis.

``relative_abundance``
    Each count as a percentage of its sample (each sample sums to 100). It suits plots, correlation (except SparCC)
    and association. It does not suit PERMANOVA or read threshold, which expect counts.

The note written beside the table states the parameters actually used. On sparse data the geometric mean lowers by
itself the share of samples in which a taxon must be present to enter the scaling factors (80% at first), until
every sample has a factor; the note gives the value it ended with and the number of taxa that defined the factors.

.. code-block:: python

    from moonstone.normalization.table import normalize

    normalized = normalize(counts_df, "geometric_mean")      # rows: features, columns: samples

Things to know
==============

Host removal is by exact taxid
    Only reads classified exactly as a chosen host are removed. A read that matches human and mouse equally well
    is classified above both (Euarchontoglires, taxid 314146) and is kept. On two control samples this was 0.1% to
    0.9% of the reads.

Disk space
    By default the trimmed reads and the per-read Kraken2 files are deleted once no step needs them, and a run
    leaves only reports and tables. With both kept, a run needs several times the size of the reads.

Differences from Moonbeam
    The workflow is Moonbeam at commit ``7a0f069``. In Moonstone the reads are found for you whatever their naming
    among the three conventions above, ``keep_intermediates`` is replaced by the two keep switches, the host-removal
    script is used from the KrakenTools package, and the report records the Moonstone version. Moonbeam remains a
    separate repository.
