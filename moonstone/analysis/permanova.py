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
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

from moonstone.analysis.association import (
    _as_list, _categorical_columns, _complete_cases, _is_categorical, _numeric_column
)
from moonstone.analysis.diversity.beta import (
    DEFAULT_PERMUTATIONS, BrayCurtis, permutation_p_value_se, resolve_n_jobs, run_permdisp, too_few_permutations_warning
)
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
# Columns of a multivariable run: one row per level and model term, with the degrees of freedom and term type added.
MODEL_RESULT_COLUMNS = RESULT_COLUMNS[:2] + ['term_type', 'df'] + RESULT_COLUMNS[2:]
NOT_APPLICABLE = 'not applicable'
PERMDISP_COLUMNS = ['permdisp_F', 'permdisp_p_value', 'permdisp_p_value_se', 'permdisp_p_value_corrected']
PERMUTATION_CHUNK = 100
# adonis2's tolerance when comparing a permuted pseudo-F with the observed one
F_TOLERANCE = np.sqrt(np.finfo(float).eps)


def read_taxonomy_counts(path: str) -> pd.DataFrame:
    """
    Read a normalized count table with one row per taxon and one column per sample, where the taxonomy is given
    in columns named after the levels (kingdom, phylum, class, order, family, genus, species). Tab- or
    comma-separated. An NCBI_taxonomy_ID column, if present, is ignored. Blank count cells stop the run: fill them
    in (with 0 if the taxon was not found) beforehand.

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
    blank = df[sample_columns].isna()
    if blank.any().any():
        samples = [str(col) for col in sample_columns if blank[col].any()]
        raise ValueError(f"{path}: {int(blank.sum().sum())} blank count cell(s) in sample(s) {', '.join(samples[:5])}"
                         f"{' ...' if len(samples) > 5 else ''}. Fill them in (with 0 if the taxon was not "
                         "found) before running.")
    return df.set_index(levels)[sample_columns]


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


def _contributions(gower: np.ndarray, basis: np.ndarray, perms: np.ndarray, parallel: bool = False) -> np.ndarray:
    """q'Gq for each column q of basis, with the rows of basis shuffled by each row of perms (perms x columns)."""
    out = np.empty((len(perms), basis.shape[1]))
    with threadpool_limits(1 if parallel else None):  # one test per core: libraries' own threads would compete
        for i, perm in enumerate(perms):
            q = basis[perm]
            out[i] = np.einsum('ij,ij->j', q, gower @ q)
    return out


def _permutation_chunks(rng, blocks: list, n: int, permutations: int):
    """Shuffles of range(n) that only move a sample within its block, PERMUTATION_CHUNK at a time."""
    for start in range(0, permutations, PERMUTATION_CHUNK):
        perms = np.tile(np.arange(n), (min(PERMUTATION_CHUNK, permutations - start), 1))
        for perm in perms:
            for block in blocks:
                perm[block] = rng.permutation(block)
        yield perms


def permanova_model(
    distmat, metadata_df: pd.DataFrame, terms: list, categorical: list = None, strata: str = None,
    by: str = 'margin', permutations: int = DEFAULT_PERMUTATIONS, seed: int = 0, dispersion: bool = True,
    dispersion_test: str = 'median', n_jobs: int = 1,
) -> pd.DataFrame:
    """
    Multivariable PERMANOVA: a distance-based linear model, as ``vegan::adonis2`` with the terms added in the order
    given. The data (not the residuals) are shuffled, as in ``adonis2``.

    Args:
        distmat: scikit-bio DistanceMatrix between samples.
        metadata_df: metadata indexed by sample. A sample with a missing value in a term or in strata is left out
          (never filled in) and counted in ``n_missing``.
        terms: metadata columns of the model. Text, boolean and ``categorical`` columns are groups; other columns
          are numbers (one slope). No interactions.
        categorical: numeric columns to treat as groups (for example 1/2/3 codes).
        strata: metadata column of blocks (for example sequencing batch): samples are only shuffled within their
          block.
        by: 'margin' tests each term after all the others, 'terms' tests them one after the other in the order
          given, each after the terms before it.
        permutations, seed: as in ``BetaDiversity.permanova``. The p-values do not depend on ``n_jobs`` or on the
          order of the samples.
        dispersion, dispersion_test: also run PERMDISP (Permutational Analysis of Multivariate Dispersions) on each
          categorical term, on the same samples. It tests one term at a time and does not adjust the spread for
          the other terms (neither does ``vegan``). Numeric terms get NaN.
        n_jobs: processes sharing the permutations and the PERMDISP runs (0 or -1: every core).

    Returns:
        One row per term: term_type, df, n_samples, n_missing, n_groups (categorical terms), pseudo_F, R2, p_value,
        p_value_se, the PERMDISP fields and permutations.
    """
    if by not in ('terms', 'margin'):
        raise ValueError(f"by must be 'terms' or 'margin', not {by!r}.")
    terms, categorical = _as_list(terms), _as_list(categorical)
    if not terms:
        raise ValueError("Name at least one variable to test.")
    model_columns = terms + _as_list(strata)
    if len(set(model_columns)) != len(model_columns):
        raise ValueError("A variable is listed twice (in the variables, the covariates or the strata).")
    unknown = [c for c in model_columns + categorical if c not in metadata_df.columns]
    if unknown:
        raise ValueError(f"Not metadata columns: {unknown}.")
    if [c for c in categorical if c not in terms]:
        raise ValueError("categorical names columns that are not in the model.")

    ids = [str(i) for i in distmat.ids]
    meta = _complete_cases(metadata_df.apply(lambda c: c.mask(c.astype(str).str.strip().eq(''))), ids, model_columns)
    meta = meta.sort_index()  # a fixed order, so the same seed gives the same p-values whatever the input order
    n = len(meta)
    n_missing = len(ids) - n
    position = {name: i for i, name in enumerate(ids)}
    keep = [position[s] for s in meta.index]
    squared = -0.5 * distmat.data[np.ix_(keep, keep)] ** 2
    gower = squared - squared.mean(0) - squared.mean(1)[:, None] + squared.mean()  # Gower-centred matrix
    total = np.trace(gower)

    kinds, blocks = [], []
    for term in terms:
        is_group = _is_categorical(meta[term], term in categorical)
        kinds.append('categorical' if is_group else 'numeric')
        columns = ([v for v, _ in _categorical_columns(term, meta[term], {})[0].values()] if is_group
                   else [_numeric_column(term, meta[term], {}, True)])
        blocks.append(np.column_stack(columns))
    df = np.array([b.shape[1] for b in blocks])
    df_residual = n - 1 - df.sum()
    if df_residual < 1:
        raise ValueError(f"{n} samples are too few for a model with {df.sum()} parameters.")

    def orthonormal(*parts):
        """Orthonormal basis of the columns (after the intercept, which Gower centring makes irrelevant)."""
        q, r = np.linalg.qr(np.column_stack([np.ones(n), *parts]))
        return q[:, 1:], np.abs(np.diag(r)).min()

    sequential, smallest = orthonormal(*blocks)
    if smallest < 1e-8:
        raise ValueError("The terms are collinear (one is a combination of the others): drop one.")
    bases = [sequential]
    if by == 'margin':  # a term's margin test is its sequential test when it comes last
        bases += [orthonormal(*blocks[:t], *blocks[t + 1:], blocks[t])[0][:, -df[t]:] for t in range(len(terms))]
    basis = np.hstack(bases)
    starts = np.cumsum(df) - df
    width = sequential.shape[1]

    def statistics(contributions):
        residual = total - contributions[:, :width].sum(1)
        used = contributions[:, width:] if by == 'margin' else contributions[:, :width]
        sum_squares = np.add.reduceat(used, starts, axis=1)
        return sum_squares, sum_squares / df / (residual / df_residual)[:, None]

    observed_ss, observed_f = (a[0] for a in statistics(_contributions(gower, basis, np.arange(n)[None])))
    blocks_of_samples = list(pd.Series(np.arange(n)).groupby(
        meta[strata].to_numpy() if strata else np.zeros(n)).indices.values())
    if max(len(b) for b in blocks_of_samples) < 2:
        raise ValueError(f"Every {strata} block has a single sample: no permutation is possible.")
    n_jobs = resolve_n_jobs(n_jobs)
    chunks = _permutation_chunks(np.random.default_rng(seed), blocks_of_samples, n, permutations)
    if n_jobs == 1:
        parts = [_contributions(gower, basis, perms) for perms in chunks]
    else:
        parts = Parallel(n_jobs=n_jobs)(delayed(_contributions)(gower, basis, perms, True) for perms in chunks)
    permuted_f = statistics(np.vstack(parts))[1]
    p_values = (1 + (permuted_f >= observed_f - F_TOLERANCE).sum(0)) / (permutations + 1)

    permdisp = {term: (np.nan, np.nan) for term in terms}
    if dispersion:
        groups = [t for t, kind in zip(terms, kinds) if kind == 'categorical']
        reduced = distmat.filter(list(meta.index))
        args = [(reduced, meta[[t]].astype(str), t, dispersion_test, permutations, seed) for t in groups]
        runs = [run_permdisp(*a) for a in args] if n_jobs == 1 else Parallel(n_jobs=n_jobs)(
            delayed(run_permdisp)(*a) for a in args)
        permdisp.update(zip(groups, runs))
    return pd.DataFrame({
        'term_type': kinds, 'df': df, 'n_samples': n, 'n_missing': n_missing,
        'n_groups': [meta[t].nunique() if k == 'categorical' else np.nan for t, k in zip(terms, kinds)],
        'pseudo_F': observed_f, 'R2': observed_ss / total, 'p_value': p_values,
        'p_value_se': [permutation_p_value_se(p, permutations) for p in p_values],
        'permdisp_F': [permdisp[t][0] for t in terms], 'permdisp_p_value': [permdisp[t][1] for t in terms],
        'permdisp_p_value_se': [permutation_p_value_se(permdisp[t][1], permutations) for t in terms],
        'permutations': permutations,
    }, index=pd.Index(terms, name='variable'))


def run_permanova(
    counts_df: pd.DataFrame, metadata_df: pd.DataFrame, columns: list = None, taxa_level: str = 'species',
    read_threshold: bool = True, correction_method: str = 'fdr_bh', max_groups: int = 15,
    permutations: int = DEFAULT_PERMUTATIONS, seed: int = 0, min_group_size: int = 2, n_jobs: int = 1,
    covariates: list = None, strata: str = None, categorical: list = None, by: str = 'margin',
) -> pd.DataFrame:
    """
    PERMANOVA with a PERMDISP dispersion check on Bray-Curtis distances, for one or all taxonomic levels.

    Args:
        counts_df: normalized counts, one row per taxon and one column per sample, indexed by the taxonomy levels
          (as returned by :func:`read_taxonomy_counts`).
        metadata_df: metadata indexed by sample. Samples are matched by ID; unmatched samples are left out.
        columns: metadata columns to test. Default: all (with covariates or strata, name them).
        taxa_level: a level of the taxonomy index ('species' by default), or 'all' for every level present.
        read_threshold: keep only the taxa above the read-threshold knee at each level before testing.
        correction_method: 'fdr_bh' (default), 'bonferroni', or None. Applied once over every test of the run.
        max_groups, permutations, seed, min_group_size, n_jobs: see ``BetaDiversity.permanova_screen``.
          The seed makes a run repeatable; it does not make the p-values more accurate (more permutations do).
        covariates, strata, categorical, by: with covariates or strata, ``columns`` and ``covariates`` become the terms
          of one model fitted with :func:`permanova_model` (every term is tested and reported, by default each one
          after all the others), with permutations restricted within ``strata``. ``max_groups`` and
          ``min_group_size`` then do not apply.

    Returns:
        One row per level and variable (or model term), with the columns of ``RESULT_COLUMNS`` (``MODEL_RESULT_COLUMNS``
        for a model). PERMDISP is ``'not applicable'`` for numeric terms. ``warnings`` lists anything the
        reader should know about that level (separated by '; ').
    """
    model = bool(covariates or strata)
    terms = list(columns or []) + _as_list(covariates)
    if model and not columns:
        raise ValueError("Name the variable(s) to test when using covariates or strata.")
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
                                  index=pd.Index(terms if model else columns or list(metadata_df.columns)))
        elif model:
            screen = permanova_model(
                BrayCurtis(level_df).beta_diversity, metadata_df, terms, categorical, strata, by, permutations, seed,
                n_jobs=n_jobs)
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
    if model and 'term_type' in results:
        numeric = results['term_type'] == 'numeric'
        results[PERMDISP_COLUMNS] = results[PERMDISP_COLUMNS].astype(object)
        results.loc[numeric, PERMDISP_COLUMNS] = NOT_APPLICABLE
    return results.reindex(columns=MODEL_RESULT_COLUMNS if model else RESULT_COLUMNS)
