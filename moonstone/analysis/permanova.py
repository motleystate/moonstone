"""
PERMANOVA from a Moonstone-normalized count table with taxonomy, at one or all taxonomic levels.

For each taxonomic level the steps are:

1. Sum the counts to that level (summing first, so reads from rare species still count towards their genus).
2. Optionally keep only the taxa above the knee of the read-threshold score
   (:mod:`moonstone.filtering.read_threshold_analysis`). The knee has been validated at species level only, so
   every other level gets a warning in the results.
3. Run PERMANOVA with a PERMDISP dispersion check for each metadata column
   (:meth:`moonstone.analysis.diversity.beta.BetaDiversity.permanova_screen`).

The p-values of all tests in one run (every level x variable) are then corrected together.

From Python::

    from moonstone.analysis.permanova import read_taxonomy_counts, run_permanova

    counts = read_taxonomy_counts("normalized_counts.tsv")
    results = run_permanova(counts, metadata_df, taxa_level="all")
"""
import logging

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from moonstone.analysis.diversity.beta import DEFAULT_PERMUTATIONS, BrayCurtis, too_few_permutations_warning
from moonstone.filtering.read_threshold_analysis import combined_score_threshold, compute_taxa_metrics

logger = logging.getLogger(__name__)

TAXA_LEVELS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
NON_SAMPLE_COLUMNS = TAXA_LEVELS + ['NCBI_taxonomy_ID']
KNEE_VALIDATED_LEVEL = 'species'
MIN_TAXA_FOR_KNEE = 3
MAX_READS_PCT_REMOVED = 10

RESULT_COLUMNS = [
    'level', 'variable', 'n_samples', 'n_missing', 'n_groups', 'pseudo_F', 'R2', 'p_value', 'p_value_se',
    'p_value_corrected', 'permdisp_F', 'permdisp_p_value', 'permdisp_p_value_se', 'permdisp_p_value_corrected',
    'permutations', 'skipped',
    'taxa_total', 'taxa_kept', 'reads_pct_kept', 'warnings',
]


def read_taxonomy_counts(path: str) -> pd.DataFrame:
    """
    Read a normalized count table with one row per taxon and one column per sample, where the taxonomy is given
    in columns named after the levels (kingdom, phylum, class, order, family, genus, species). Tab- or
    comma-separated. An NCBI_taxonomy_ID column, if present, is ignored.

    Returns:
        Dataframe indexed by the taxonomy levels found (a MultiIndex), with one numeric column per sample.
    """
    with open(path) as handle:
        sep = '\t' if '\t' in handle.readline() else ','
    df = pd.read_csv(path, sep=sep)
    levels = [lvl for lvl in TAXA_LEVELS if lvl in df.columns]
    if not levels:
        raise ValueError(f"{path}: no taxonomy column found. Expected columns among: {', '.join(TAXA_LEVELS)}.")
    sample_columns = [col for col in df.columns if col not in NON_SAMPLE_COLUMNS]
    non_numeric = [col for col in sample_columns if not pd.api.types.is_numeric_dtype(df[col])]
    if non_numeric:
        raise ValueError(f"{path}: sample columns must be numeric, but these are not: {', '.join(non_numeric[:5])}"
                         f"{' ...' if len(non_numeric) > 5 else ''}.")
    return df.set_index(levels)[sample_columns].fillna(0)


def sum_to_level(counts_df: pd.DataFrame, level: str) -> tuple:
    """
    Sum the counts of a taxonomy-indexed table to one level.

    Taxa without a name at that level are pooled as 'unclassified <level>'.

    Returns:
        The summed dataframe (one row per taxon of that level) and a list of warnings.
    """
    if level not in counts_df.index.names:
        raise ValueError(f"Level '{level}' not in the count table. Available: {', '.join(counts_df.index.names)}.")
    names = counts_df.index.get_level_values(level)
    warnings = []
    n_unnamed = int(pd.isna(names).sum())
    if n_unnamed:
        names = names.fillna(f'unclassified {level}')
        warnings.append(f"{n_unnamed} row(s) without a {level} name pooled as 'unclassified {level}'")
    summed = counts_df.groupby(np.asarray(names)).sum()
    summed.index.name = level
    return summed, warnings


def read_threshold_filter(level_df: pd.DataFrame, level: str) -> tuple:
    """
    Keep the taxa above the knee of the read-threshold score (prevalence x log10(mean reads + 1)).

    Returns:
        The filtered dataframe, a dict with taxa_total, taxa_kept and reads_pct_kept, and a list of warnings.
    """
    warnings = []
    if level != KNEE_VALIDATED_LEVEL:
        warnings.append(f"read-threshold knee validated at {KNEE_VALIDATED_LEVEL} level only, "
                        f"not at {level} level (issue #132)")
    taxa_total = level_df.shape[0]
    if taxa_total < MIN_TAXA_FOR_KNEE:
        warnings.append(f"read-threshold not applied: only {taxa_total} taxa")
        filtered = level_df
    else:
        keep, _ = combined_score_threshold(compute_taxa_metrics(level_df))
        filtered = level_df.loc[keep[keep].index]
        if filtered.empty:
            warnings.append("read-threshold not applied: no taxon above the knee")
            filtered = level_df
    total_reads = level_df.to_numpy().sum()
    summary = {
        'taxa_total': taxa_total,
        'taxa_kept': filtered.shape[0],
        'reads_pct_kept': 100 * filtered.to_numpy().sum() / total_reads if total_reads > 0 else np.nan,
    }
    if 100 - summary['reads_pct_kept'] > MAX_READS_PCT_REMOVED:
        warnings.append(f"read-threshold removed {100 - summary['reads_pct_kept']:.1f}% of reads "
                        f"(kept {summary['taxa_kept']} of {taxa_total} taxa)")
    return filtered, summary, warnings


def _drop_empty_samples(level_df: pd.DataFrame) -> tuple:
    """Bray-Curtis is undefined for a sample with no reads: leave such samples out."""
    empty = level_df.columns[level_df.sum(axis=0) == 0]
    warnings = []
    if len(empty) > 0:
        warnings.append(f"{len(empty)} sample(s) with no reads left out: {', '.join(map(str, empty))}")
        level_df = level_df.drop(columns=empty)
    return level_df, warnings


def run_permanova(
    counts_df: pd.DataFrame, metadata_df: pd.DataFrame, columns: list = None, taxa_level: str = 'species',
    read_threshold: bool = True, correction_method: str = 'fdr_bh', max_groups: int = 15,
    permutations: int = DEFAULT_PERMUTATIONS, seed: int = 0, min_group_size: int = 2, n_jobs: int = 1,
) -> pd.DataFrame:
    """
    PERMANOVA with a PERMDISP dispersion check on Bray-Curtis distances, for one or all taxonomic levels.

    Args:
        counts_df: normalized counts, one row per taxon and one column per sample, indexed by the taxonomy levels
          (as returned by :func:`read_taxonomy_counts`).
        metadata_df: metadata indexed by sample. Samples are matched by ID; unmatched samples are left out.
        columns: metadata columns to test. Default: all.
        taxa_level: a level of the taxonomy index ('species' by default), or 'all' for every level present.
        read_threshold: keep only the taxa above the read-threshold knee at each level before testing.
        correction_method: 'fdr_bh' (default), 'bonferroni', or None. Applied once over every test of the run.
        max_groups, permutations, seed, min_group_size, n_jobs: see ``BetaDiversity.permanova_screen``.
          The seed makes a run repeatable; it does not make the p-values more accurate (more permutations do).

    Returns:
        One row per level and variable, with the columns of ``RESULT_COLUMNS``. ``warnings`` lists anything the
        reader should know about that level (separated by '; ').
    """
    available = [lvl for lvl in TAXA_LEVELS if lvl in counts_df.index.names]
    if taxa_level == 'all':
        levels = available
    elif taxa_level in available:
        levels = [taxa_level]
    else:
        raise ValueError(f"taxa_level '{taxa_level}' not in the count table. Use one of: "
                         f"{', '.join(available + ['all'])}.")

    samples = [s for s in counts_df.columns if s in metadata_df.index]
    n_unmatched = counts_df.shape[1] - len(samples)
    if n_unmatched:
        logger.warning("%i count sample(s) have no metadata and are left out.", n_unmatched)
    if not samples:
        raise ValueError("No sample is found in both the count table and the metadata.")
    counts_df = counts_df[samples]

    tables = []
    for level in levels:
        level_df, warnings = sum_to_level(counts_df, level)
        if read_threshold:
            level_df, summary, filter_warnings = read_threshold_filter(level_df, level)
            warnings += filter_warnings
        else:
            summary = {'taxa_total': level_df.shape[0], 'taxa_kept': level_df.shape[0], 'reads_pct_kept': 100.0}
        level_df, empty_warnings = _drop_empty_samples(level_df)
        warnings += empty_warnings
        for warning in warnings:
            logger.warning("%s: %s", level, warning)

        if level_df.shape[0] < 2:
            # With one taxon, Bray-Curtis only compares total reads, not community composition.
            screen = pd.DataFrame({'skipped': 'fewer than 2 taxa, composition cannot be compared'},
                                  index=pd.Index(columns or list(metadata_df.columns)))
        else:
            screen = BrayCurtis(level_df).permanova_screen(
                metadata_df, columns=columns, correction_method=None, max_groups=max_groups,
                permutations=permutations, seed=seed, min_group_size=min_group_size, n_jobs=n_jobs,
            )
        screen = screen.rename_axis('variable').reset_index()
        screen.insert(0, 'level', level)
        for key, value in summary.items():
            screen[key] = value
        screen['warnings'] = '; '.join(warnings)
        tables.append(screen)

    results = pd.concat(tables, ignore_index=True)
    if correction_method is not None:
        warning = too_few_permutations_warning(int(results['p_value'].notna().sum()), permutations)
        if warning:
            logger.warning("%s.", warning)
            results['warnings'] = [f"{w}; {warning}" if w else warning for w in results['warnings']]
    for pval_col in ['p_value', 'permdisp_p_value']:
        corrected = pd.Series(np.nan, index=results.index)
        tested = results[pval_col].dropna()
        if correction_method is not None and not tested.empty:
            corrected[tested.index] = multipletests(tested, method=correction_method)[1]
        results[f'{pval_col}_corrected'] = corrected
    return results.reindex(columns=RESULT_COLUMNS)
