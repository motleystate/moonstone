"""
Normalize a Moonstone count table, file to file.

    python -m moonstone.normalization.table counts.csv normalized.csv --method geometric_mean

This step is separate from the reads workflow on purpose: it takes any Moonstone count table (features in rows,
samples in columns), so the choice of normalization can change without touching how the counts were obtained.
"""
import argparse
import os
import sys

import pandas as pd

from moonstone import __version__
from moonstone.normalization.counts.geometric_mean import GeometricMeanNormalization

METHODS = ('geometric_mean', 'relative_abundance')
# What each output is for. Shown in the help, written in the methods note.
SUITS = {
    'geometric_mean': 'Values stay on a count scale: suitable for every Moonstone analysis.',
    'relative_abundance': 'Values are percentages: suitable for plots, correlation (except SparCC) and association; '
                          'not for PERMANOVA or read threshold, which expect counts.',
}
ID_COLUMNS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'NCBI_taxonomy_ID')


def _normalize(counts: pd.DataFrame, method: str) -> tuple:
    """(normalized table, sentence describing what was done)."""
    if method not in METHODS:
        raise ValueError(f"Unknown normalization {method!r}. Available: {', '.join(METHODS)}.")
    empty = [str(sample) for sample in counts.columns[(counts == 0).all()]]
    if empty:
        raise ValueError(f"Sample(s) without any read cannot be normalized: {', '.join(empty)}. Remove them first.")
    if method == 'relative_abundance':
        return counts * 100 / counts.sum(), ('Each count was divided by the total of its sample and multiplied by 100 '
                                             '(relative abundance in percent; each sample sums to 100).')
    normalization = GeometricMeanNormalization(counts)
    normalized = normalization.normalized_df
    # the class lowers its zero filter by itself until every sample has a scaling factor: report the value used
    reference = len(counts) - len(normalization.removed_zero_df)
    return normalized, ('Each sample was divided by its scaling factor, the median of the ratios of its counts to the '
                        'geometric mean of each feature (median of ratios, as in DESeq2; '
                        'moonstone GeometricMeanNormalization). The scaling factors were computed on the features '
                        f'non-zero in at least {normalization.zero_threshold}% of the samples: {reference} of '
                        f'{len(counts)}. All features were then divided by them.'), normalization.zero_threshold


def normalize(counts: pd.DataFrame, method: str) -> pd.DataFrame:
    """
    Normalize a table with features in rows and samples in columns.

    :param method: ``geometric_mean``: median of ratios to the geometric mean of each feature, as in DESeq2
        (:class:`~moonstone.normalization.counts.geometric_mean.GeometricMeanNormalization`); values stay on a
        count scale. ``relative_abundance``: percentage of the sample total (each sample sums to 100).
    :return: a new table with the same rows and columns.
    """
    return _normalize(counts, method)[0]


def read_counts(path: str) -> pd.DataFrame:
    """
    Read a Moonstone count table (tab- or comma-separated). The columns among kingdom ... species and
    NCBI_taxonomy_ID identify the rows, otherwise the first column does; every other column is a sample.
    """
    with open(path) as handle:
        separator = '\t' if '\t' in handle.readline() else ','
    table = pd.read_csv(path, sep=separator)
    identifiers = [column for column in table.columns if column in ID_COLUMNS] or [table.columns[0]]
    counts = table.set_index(identifiers)
    text = [str(column) for column in counts.columns if not pd.api.types.is_numeric_dtype(counts[column])]
    if text:
        raise ValueError(f"{path}: sample columns must be numeric, but these are not: {', '.join(text[:5])}.")
    blank = [str(column) for column in counts.columns[counts.isna().any()]]
    if blank:
        raise ValueError(f"{path}: blank cell(s) in sample(s) {', '.join(blank[:5])}. Fill them in (with 0 if the "
                         'feature was not found).')
    return counts


def normalize_table(path: str, output: str, method: str) -> dict:
    """
    Normalize the table in ``path`` and write it to ``output`` (tab-separated if it ends with .tsv, else
    comma-separated), with the same identifying columns, so it is read by the same analyses. A methods note
    is written to ``<output>.methods.txt``. Neither file is ever overwritten.

    :return: {'samples': n, 'features': n, 'method': method} and, for the geometric mean, 'zero_threshold': the
        percentage of samples in which a feature had to be non-zero to enter the scaling factors.
    """
    counts = read_counts(path)
    normalized, done, *threshold = _normalize(counts, method)
    summary = {'samples': counts.shape[1], 'features': counts.shape[0], 'method': method}
    if threshold:
        summary['zero_threshold'] = threshold[0]
    with open(output, 'x') as handle, open(output + '.methods.txt', 'x') as methods:
        normalized.reset_index().to_csv(handle, index=False, sep='\t' if output.endswith('.tsv') else ',')
        methods.write(f'{os.path.abspath(path)} ({summary["samples"]} samples, {summary["features"]} features) was '
                      f'normalized with Moonstone {__version__} (moonstone.normalization.table, method '
                      f'{method.replace("_", " ")}). {done} {SUITS[method]}\n')
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='python -m moonstone.normalization.table',
        description='Normalize a Moonstone count table (features in rows, samples in columns).')
    parser.add_argument('table', help='Count table: taxonomy columns (kingdom ... species) or a first column naming '
                                      'the features, then one column of counts per sample')
    parser.add_argument('output', help='Table to write; or an existing folder, where it is named counts_<method>.csv')
    parser.add_argument('--method', choices=METHODS, required=True,
                        help='geometric_mean: median of ratios, as in DESeq2. ' + SUITS['geometric_mean']
                             + ' relative_abundance: percentage of each sample. ' + SUITS['relative_abundance'])
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    output = os.path.join(args.output, f'counts_{args.method}.csv') if os.path.isdir(args.output) else args.output
    try:
        summary = normalize_table(args.table, output, args.method)
    except FileExistsError as error:
        print(f'Error: {error.filename} already exists; it is never overwritten.', file=sys.stderr)
        return 1
    except (OSError, ValueError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    print(f'{summary["samples"]} samples, {summary["features"]} features normalized ({args.method}) to {output}.')
    print(SUITS[args.method])
    return 0


if __name__ == '__main__':
    sys.exit(main())
