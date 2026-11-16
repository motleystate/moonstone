"""``moonstone-reads table``: the table written by the reads workflow, in Moonstone's layout. Counts are not changed."""
import argparse
import os
import sys

import pandas as pd

from moonstone import __version__
from moonstone.parsers import SunbeamKraken2Parser

LINEAGE_COLUMNS = ('taxonomy', 'Consensus Lineage')       # as headed by the reads workflow, and by Sunbeam
# kraken-biom names a sample column after its report file: S1.brackenCorrected_taxa (Bracken) or S1.taxa (Kraken2)
SAMPLE_SUFFIXES = ('.brackenCorrected_taxa', '.taxa')
ID_COLUMN = SunbeamKraken2Parser.new_otu_id_name


class _WorkflowTableParser(SunbeamKraken2Parser):
    """SunbeamKraken2Parser reading the lineage from the column that the table has."""

    def __init__(self, path, lineage_column):
        self.taxa_column = lineage_column
        super().__init__(path)


def read_moonbeam_table(path: str) -> pd.DataFrame:
    """
    Read a table of the reads workflow (``all_samples.bracken.tsv`` or ``all_samples.kraken.tsv``).

    :return: counts indexed by taxonomy (kingdom ... species), with the NCBI_taxonomy_ID column first, then one
        column per sample; sample names are those of the read files (the ``.brackenCorrected_taxa`` or ``.taxa``
        ending is removed).
    :raise ValueError: no lineage column, or two samples with the same name.
    """
    with open(path) as handle:
        handle.readline()
        header = handle.readline().rstrip('\n').split('\t')
    lineage = next((column for column in LINEAGE_COLUMNS if column in header), None)
    if lineage is None:
        raise ValueError(f"{path}: not a table of the reads workflow (no '{LINEAGE_COLUMNS[0]}' column on its second "
                         'line).')
    table = _WorkflowTableParser(path, lineage).dataframe

    def clean(name):
        return next((name[:-len(suffix)] for suffix in SAMPLE_SUFFIXES if name.endswith(suffix)), name)

    table.columns = [clean(str(column)) for column in table.columns]
    repeated = sorted(table.columns[table.columns.duplicated()])
    if repeated:
        raise ValueError(f"{path}: sample name(s) used more than once: {', '.join(repeated)}.")
    return table


def to_moonstone_table(path: str, output: str) -> dict:
    """
    Write the counts of a reads-workflow table in Moonstone's layout: comma-separated, the taxonomy in columns
    kingdom ... species, then NCBI_taxonomy_ID, then one column per sample. The counts are not changed (they are
    the raw, Bracken-corrected counts when the input is the Bracken table). The result is read by PERMANOVA,
    correlation, association and read threshold; normalize it with :mod:`moonstone.normalization.table`.

    A methods note is written to ``<output>.methods.txt``. Neither file is ever overwritten.

    :return: {'samples': number of samples, 'taxa': number of rows}
    """
    table = read_moonbeam_table(path)
    samples = [column for column in table.columns if column != ID_COLUMN]
    counts = table[samples]
    if (counts == counts.round()).all().all():
        table = table.astype({sample: 'int64' for sample in samples})       # 582.0 -> 582: the same count
    summary = {'samples': len(samples), 'taxa': len(table)}
    with open(output, 'x') as handle, open(output + '.methods.txt', 'x') as methods:
        table.reset_index().to_csv(handle, index=False)
        methods.write(
            f'Counts of {os.path.abspath(path)}, written in the Moonstone table layout by Moonstone {__version__} '
            f'(moonstone.reads.table): {summary["samples"]} samples, {summary["taxa"]} taxa. The counts are those of '
            'the input, not normalized. Sample names are those of the read files. Genus and species are joined in '
            'the species column.\n')
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='moonstone-reads table',
        description="Write a table of the reads workflow in Moonstone's layout. The counts are not changed.")
    parser.add_argument('table', help='all_samples.bracken.tsv (or all_samples.kraken.tsv) of a reads run')
    parser.add_argument('output', help='Table to write (comma-separated); or an existing folder, where it is named '
                                       'counts_raw.csv. An existing table is never overwritten')
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    output = os.path.join(args.output, 'counts_raw.csv') if os.path.isdir(args.output) else args.output
    try:
        summary = to_moonstone_table(args.table, output)
    except FileExistsError as error:
        print(f'Error: {error.filename} already exists; it is never overwritten.', file=sys.stderr)
        return 1
    except (OSError, ValueError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    print(f'{summary["samples"]} samples, {summary["taxa"]} taxa written to {output} (counts unchanged).')
    return 0


if __name__ == '__main__':
    sys.exit(main())
