"""
From sequencing reads to count tables: adapter trimming, host-read removal, Kraken2 classification and Bracken.

The workflow in ``workflow/`` is the Moonbeam pipeline (https://github.com/motleystate/moonbeam2, commit 7a0f069)
moved into Moonstone. It needs the ``reads`` extra (Snakemake) and Kraken2, Bracken, Trimmomatic and KrakenTools
on PATH (see environment-hpc.yml).
"""
import os

WORKFLOW_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'workflow')
