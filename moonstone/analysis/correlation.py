"""
Correlation between features (taxa, genes...) of compositional count tables.

Relative abundances are compositional, so correlations computed on them directly are partly artefacts (one feature
going up pushes the others down). Every method here works on log-ratios instead:

- ``"spearman"`` (default), ``"pearson"`` and ``"kendall"``: the usual coefficients, computed on the Centred
  Log-Ratio (CLR) transform of each table (:func:`moonstone.transformers.compositional.clr`).
- ``"rho"``: proportionality (Lovell et al., 2015; Erb and Notredame, 2016), 1 - var(clr_a - clr_b) /
  (var(clr_a) + var(clr_b)). It is 1 when two features keep a constant ratio.
- ``"sparcc"``: SparCC (Sparse Correlations for Compositional data; Friedman and Alm, 2012), which estimates the
  correlation of the underlying absolute abundances, assuming most pairs are not correlated. It needs count-like
  values, because it resamples each sample from a Dirichlet distribution built on its counts. For large tables,
  FastSpar (a C++ implementation, installable from bioconda) is much faster; it is not called by Moonstone.

p-values come from a permutation null, for every method. The CLR itself makes features slightly negatively
correlated (about -1 / (number of features - 1)), so the textbook p-values of a coefficient computed on CLR values
give too many false positives: in simulations with 20 independent features, 30% of datasets had at least one
"significant" pair at a 5% false discovery rate, instead of 5%. Here each feature is shuffled across samples
independently, the table is transformed again and the coefficient recomputed, ``permutations`` times (as CoNet's
ReBoot does, Faust et al., 2012). The observed coefficient is compared with the mean and standard deviation of these
shuffled coefficients (a z-score), which brings the false-positive rate back to 5% in the same simulations and, unlike
counting exceedances, can give p-values much smaller than 1 / permutations. For SparCC each permutation uses 5
Dirichlet draws, as FastSpar does. Kendall is computed pair by pair and is slow on large tables.

With ``second_counts`` (for example species and resistance genes), each table is transformed on its own and only
pairs across the two tables are tested: two tables measured in different units do not share one total, so they must
not be merged into one composition. Only the CLR coefficients are available across tables.

p-values are corrected for multiple testing over all tested pairs (Benjamini-Hochberg by default).

Analyses are cross-sectional: one sample per subject. With ``metadata`` and ``group_column`` (for example a patient
ID), the function stops if a subject has more than one sample; it never averages or picks samples.

From Python::

    from moonstone.analysis.correlation import feature_correlation, filter_edges

    edges = feature_correlation(normalized_counts_df)
    significant = filter_edges(edges, alpha=0.05)
"""
import logging
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from moonstone.transformers.compositional import _check_table, _clr_array, clr

logger = logging.getLogger(__name__)

METHODS = ('spearman', 'pearson', 'kendall', 'rho', 'sparcc')
CLR_METHODS = ('spearman', 'pearson', 'kendall', 'rho')
CROSS_TABLE_METHODS = ('spearman', 'pearson', 'kendall')
DEFAULT_PERMUTATIONS = 999
MIN_PERMUTATIONS = 100
MAX_SHUFFLE_TRIES = 100
EDGE_COLUMNS = ['feature_a', 'feature_b', 'method', 'r', 'z', 'p_value', 'q_value', 'n_samples',
                'prevalence_a', 'prevalence_b']
MIN_SAMPLES = 3
# SparCC defaults follow Friedman and Alm (2012) and FastSpar
SPARCC_DRAWS = 20
SPARCC_PERMUTATION_DRAWS = 5
SPARCC_EXCLUSION_ITERATIONS = 10
SPARCC_THRESHOLD = 0.1
SPARCC_MIN_PARTNERS = 3


def _feature_labels(index: pd.Index) -> list:
    """Readable, unique feature names: the last level of a taxonomy MultiIndex when unique, else all levels."""
    if isinstance(index, pd.MultiIndex):
        last = [str(t[-1]) for t in index]
        if len(set(last)) == len(last):
            return last
        return [';'.join(str(level) for level in t) for t in index]
    labels = [str(i) for i in index]
    if len(set(labels)) != len(labels):
        raise ValueError("Feature names (row labels) must be unique.")
    return labels


def _align_samples(counts: pd.DataFrame, second_counts: pd.DataFrame):
    """Keep the samples present in both tables, matched by name, in the order of the first table."""
    first = [str(s) for s in counts.columns]
    second = [str(s) for s in second_counts.columns]
    for name, ids in (('counts', first), ('second_counts', second)):
        if len(set(ids)) != len(ids):
            raise ValueError(f"Sample names are used more than once in {name}.")
    common = [s for s in first if s in set(second)]
    if not common:
        raise ValueError("No sample name is found in both tables.")
    dropped = (len(first) - len(common), len(second) - len(common))
    if any(dropped):
        logger.warning("%d sample(s) of counts and %d of second_counts are not in the other table and are left out.",
                       *dropped)
    counts = counts.set_axis(first, axis=1)[common]
    second_counts = second_counts.set_axis(second, axis=1)[common]
    return counts, second_counts


def check_one_sample_per_subject(samples, metadata: pd.DataFrame, group_column: str):
    """
    Stop if a subject (a value of group_column) has more than one of the given samples.

    Args:
        samples: sample names analysed together.
        metadata: one row per sample, sample names as index.
        group_column: metadata column identifying the subject, for example a patient ID.
    """
    if group_column not in metadata.columns:
        raise ValueError(f"group_column {group_column!r} is not a metadata column.")
    meta = metadata.set_axis([str(i) for i in metadata.index], axis=0)
    samples = [str(s) for s in samples]
    missing = [s for s in samples if s not in meta.index]
    if missing:
        raise ValueError(f"{len(missing)} sample(s) have no metadata row, so their subject cannot be checked, for "
                         f"example {missing[:5]}.")
    groups = meta.loc[samples, group_column]
    if groups.isna().any():
        blank = groups[groups.isna()].index.tolist()
        raise ValueError(f"{len(blank)} sample(s) have no value in {group_column!r}, for example {blank[:5]}.")
    repeated = groups[groups.duplicated(keep=False)]
    if not repeated.empty:
        subjects = repeated.unique().tolist()
        raise ValueError(
            f"{len(subjects)} subject(s) in {group_column!r} have more than one sample, for example {subjects[:5]}. "
            "The analysis needs one sample per subject: select one time point or condition first."
        )


def _pearson_matrix(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson correlation between every row of x and every row of y; NaN for constant rows."""
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    with np.errstate(invalid='ignore', divide='ignore'):
        xs = xc / np.linalg.norm(xc, axis=1, keepdims=True)
        ys = yc / np.linalg.norm(yc, axis=1, keepdims=True)
        r = xs @ ys.T
    return np.clip(r, -1, 1)


def _kendall_matrix(x: np.ndarray, y: np.ndarray, upper_only: bool) -> np.ndarray:
    """Kendall tau between rows of x and rows of y (upper triangle only for one table)."""
    tau = np.full((x.shape[0], y.shape[0]), np.nan)
    for i in range(x.shape[0]):
        for j in range(i + 1 if upper_only else 0, y.shape[0]):
            tau[i, j] = stats.kendalltau(x[i], y[j]).statistic
    return tau


def _rho_matrix(x: np.ndarray) -> np.ndarray:
    """Proportionality rho between rows of a CLR table: 2 cov(a, b) / (var(a) + var(b))."""
    cov = np.cov(x)
    var = np.diag(cov)
    with np.errstate(invalid='ignore', divide='ignore'):
        return 2 * cov / (var[:, None] + var[None, :])


def _sparcc_basis(log_fractions: np.ndarray, exclusion_iterations: int, threshold: float) -> np.ndarray:
    """One SparCC estimate from log-fractions (features x samples)."""
    d = log_fractions.shape[0]
    cov = np.cov(log_fractions)
    var = np.diag(cov)
    variation = var[:, None] + var[None, :] - 2 * cov
    included = ~np.eye(d, dtype=bool)
    for iteration in range(exclusion_iterations + 1):
        partners = included.sum(axis=1)
        m = np.diag(partners.astype(float)) + included
        t = np.where(included, variation, 0).sum(axis=1)
        basis_var = np.linalg.solve(m, t)
        basis_var[basis_var <= 0] = 1e-4
        std = np.sqrt(basis_var)
        rho = (basis_var[:, None] + basis_var[None, :] - variation) / (2 * std[:, None] * std[None, :])
        np.fill_diagonal(rho, 1)
        if iteration == exclusion_iterations:
            break
        # exclude the strongest remaining pair, unless one of its features would keep too few partners
        candidates = np.where(included & (partners[:, None] > SPARCC_MIN_PARTNERS)
                              & (partners[None, :] > SPARCC_MIN_PARTNERS), np.abs(rho), 0)
        i, j = np.unravel_index(np.argmax(candidates), candidates.shape)
        if candidates[i, j] <= threshold:
            break
        included[i, j] = included[j, i] = False
    return np.clip(rho, -1, 1)


def _sparcc(values: np.ndarray, rng: np.random.Generator, draws: int, exclusion_iterations: int,
            threshold: float) -> np.ndarray:
    """SparCC correlation matrix: median over Dirichlet draws of the fractions of each sample (features x samples)."""
    estimates = []
    for _ in range(draws):
        gamma = rng.gamma(values + 1)
        log_fractions = np.log(gamma / gamma.sum(axis=0))
        estimates.append(_sparcc_basis(log_fractions, exclusion_iterations, threshold))
    return np.median(estimates, axis=0)


def _shuffle_features(values: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle each feature (row) across samples independently, redrawing if a sample ends up with only zeros."""
    for _ in range(MAX_SHUFFLE_TRIES):
        shuffled = rng.permuted(values, axis=1)
        if (shuffled > 0).any(axis=0).all():
            return shuffled
    raise ValueError("The table is too sparse to permute: shuffled samples keep ending up with only zeros. Filter out "
                     "rare features first.")


def _permutation_null(statistic, values: list, rng: np.random.Generator, permutations: int):
    """
    Mean and standard deviation of the statistic when every feature of every table is shuffled independently across
    samples and the tables are transformed again (renormalisation, as in CoNet's ReBoot, Faust et al., 2012).
    """
    nulls = np.array([statistic(*[_shuffle_features(v, rng) for v in values]) for _ in range(permutations)])
    with warnings.catch_warnings():  # pairs never computed (lower triangle) are all-NaN slices
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(nulls, axis=0), np.nanstd(nulls, axis=0, ddof=1)


def feature_correlation(
    counts: pd.DataFrame, second_counts: pd.DataFrame = None, method: str = 'spearman',
    zeros: str = 'multiplicative', delta: float = None, pseudocount: float = None,
    correction_method: str = 'fdr_bh', permutations: int = DEFAULT_PERMUTATIONS, seed: int = 0,
    metadata: pd.DataFrame = None, group_column: str = None,
    sparcc_draws: int = SPARCC_DRAWS, sparcc_exclusion_iterations: int = SPARCC_EXCLUSION_ITERATIONS,
    sparcc_threshold: float = SPARCC_THRESHOLD,
) -> pd.DataFrame:
    """
    Correlation between every pair of features, with p-values corrected for multiple testing.

    Filter the tables first (for example on prevalence): the CLR depends on which features are kept, and every
    extra feature adds tests to the correction.

    Args:
        counts: table with features as rows and samples as columns (counts, normalised counts or relative
          abundances; SparCC needs counts).
        second_counts: optional second table with the same layout; then only pairs across the two tables are tested.
          Samples are matched by name.
        method: 'spearman' (default), 'pearson', 'kendall', 'rho' or 'sparcc'.
        zeros, delta, pseudocount: zero handling for the CLR methods, see
          :func:`moonstone.transformers.compositional.clr`. Not used by SparCC.
        correction_method: 'fdr_bh' (Benjamini-Hochberg, default), 'bonferroni', or None for no correction.
        permutations: number of permutations used to estimate the null distribution of every coefficient
          (at least 100, because fewer give a noisy null and false positives; 999 by default).
        seed: random seed, for repeatable permutations and SparCC draws.
        metadata, group_column: if both are given, stop when a subject has more than one sample.
        sparcc_draws, sparcc_exclusion_iterations, sparcc_threshold: SparCC settings (Dirichlet draws, rounds of
          excluding strongly correlated pairs, and the correlation above which a pair is excluded). Permutations use
          5 draws each, as FastSpar does.

    Returns:
        One row per pair of features: feature_a, feature_b, method, r, z (distance of r from the permutation null,
        in standard deviations), p_value, q_value (NaN when not corrected), n_samples, prevalence_a, prevalence_b
        (share of samples where the feature is non-zero). r, z and p_value are NaN for a constant feature.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}.")
    if second_counts is not None and method not in CROSS_TABLE_METHODS:
        raise ValueError(f"Only {CROSS_TABLE_METHODS} can correlate two tables; {method!r} assumes features of one "
                         "composition.")
    if method == 'sparcc' and (delta is not None or pseudocount is not None or zeros != 'multiplicative'):
        raise ValueError("zeros, delta and pseudocount apply to the CLR methods; SparCC handles zeros itself.")
    if correction_method not in (None, 'fdr_bh', 'bonferroni'):
        raise ValueError(f"correction_method must be 'fdr_bh', 'bonferroni' or None, got {correction_method!r}.")
    if permutations < MIN_PERMUTATIONS:
        raise ValueError(f"At least {MIN_PERMUTATIONS} permutations are needed to estimate the null distribution, "
                         f"got {permutations}.")

    if second_counts is not None:
        counts, second_counts = _align_samples(counts, second_counts)
    if (metadata is None) != (group_column is None):
        raise ValueError("Give both metadata and group_column to check for repeated subjects, or neither.")
    if group_column is not None:
        check_one_sample_per_subject(counts.columns, metadata, group_column)

    values_a = _check_table(counts)
    values_b = _check_table(second_counts) if second_counts is not None else None
    n = values_a.shape[1]
    if n < MIN_SAMPLES:
        raise ValueError(f"At least {MIN_SAMPLES} samples are needed, got {n}.")
    labels_a = _feature_labels(counts.index)
    labels_b = _feature_labels(second_counts.index) if second_counts is not None else labels_a
    rng = np.random.default_rng(seed)

    raw = [values_a.to_numpy()] + ([values_b.to_numpy()] if values_b is not None else [])
    if method == 'sparcc':
        totals = raw[0].sum(axis=0)
        if np.allclose(totals, 1, atol=1e-3) or totals.max() <= 1:
            logger.warning("Every sample sums to about 1: these look like proportions. SparCC needs count-like values "
                           "(reads); on proportions its estimates are close to meaningless.")

    clr(values_a, zeros=zeros, delta=delta, pseudocount=pseudocount)  # validates the zero arguments once

    def transform(table):
        return _clr_array(table, zeros, delta, pseudocount)

    def statistic(*tables, draws=SPARCC_PERMUTATION_DRAWS):
        if method == 'sparcc':
            return _sparcc(tables[0], rng, draws, sparcc_exclusion_iterations, sparcc_threshold)
        transformed = [transform(table) for table in tables]
        x, y = transformed[0], transformed[-1]
        if method == 'rho':
            return _rho_matrix(x)
        if method == 'kendall':
            return _kendall_matrix(x, y, upper_only=len(tables) == 1)
        if method == 'spearman':
            x, y = stats.rankdata(x, axis=1), stats.rankdata(y, axis=1)
        return _pearson_matrix(x, y)

    r = statistic(*raw, draws=sparcc_draws)
    null_mean, null_sd = _permutation_null(statistic, raw, rng, permutations)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = (r - null_mean) / null_sd
    z[~np.isfinite(z)] = np.nan
    p = 2 * stats.norm.sf(np.abs(z))

    if values_b is None:
        rows, cols = np.triu_indices(len(labels_a), k=1)
    else:
        rows, cols = np.indices((len(labels_a), len(labels_b))).reshape(2, -1)
    prevalence_a = (values_a > 0).mean(axis=1).to_numpy()
    prevalence_b = (values_b > 0).mean(axis=1).to_numpy() if values_b is not None else prevalence_a
    edges = pd.DataFrame({
        'feature_a': np.asarray(labels_a, dtype=object)[rows],
        'feature_b': np.asarray(labels_b, dtype=object)[cols],
        'method': method,
        'r': r[rows, cols],
        'z': z[rows, cols],
        'p_value': p[rows, cols],
        'q_value': np.full(len(rows), np.nan),
        'n_samples': n,
        'prevalence_a': prevalence_a[rows],
        'prevalence_b': prevalence_b[cols],
    }, columns=EDGE_COLUMNS)
    edges.attrs['two_tables'] = values_b is not None

    n_undefined = int(edges['r'].isna().sum())
    if n_undefined:
        logger.warning("%d pair(s) involve a constant feature: r and p-value are NaN and they are not tested.",
                       n_undefined)
    tested = edges['p_value'].notna()
    if correction_method is not None and tested.any():
        edges.loc[tested, 'q_value'] = multipletests(edges.loc[tested, 'p_value'], method=correction_method)[1]
    return edges


def filter_edges(edges: pd.DataFrame, alpha: float = 0.05, min_abs_r: float = None) -> pd.DataFrame:
    """
    Keep the pairs whose corrected p-value (q_value) is below alpha, and optionally whose absolute r is at least min_abs_r.

    Both positive and negative correlations are kept. If the p-values were not corrected (q_value is NaN), the raw
    p-value is used and a warning is logged.
    """
    if edges['q_value'].notna().any():
        keep = edges['q_value'] < alpha
    else:
        logger.warning("No corrected p-values (q_value): filtering on uncorrected p-values.")
        keep = edges['p_value'] < alpha
    if min_abs_r is not None:
        keep &= edges['r'].abs() >= min_abs_r
    return edges[keep].reset_index(drop=True)


def edges_to_matrix(edges: pd.DataFrame, value: str = 'r') -> pd.DataFrame:
    """Square (one table) or rectangular (two tables) matrix of one edge column, for heatmaps."""
    matrix = edges.pivot(index='feature_a', columns='feature_b', values=value)
    features_a = list(dict.fromkeys(edges['feature_a']))
    features_b = list(dict.fromkeys(edges['feature_b']))
    two_tables = edges.attrs.get('two_tables', not set(features_a) & set(features_b))
    if not two_tables:  # one table: fill both triangles and the diagonal
        features = list(dict.fromkeys(features_a + features_b))
        matrix = matrix.reindex(index=features, columns=features)
        matrix = matrix.combine_first(matrix.T)
        if value == 'r':
            for feature in features:
                matrix.loc[feature, feature] = 1.0
        return matrix.loc[features, features]
    return matrix.reindex(index=features_a, columns=features_b)


def to_gephi(edges: pd.DataFrame) -> pd.DataFrame:
    """
    Edge list in the column names Gephi and Cytoscape import directly: Source, Target, Weight (absolute r), Sign, r,
    q_value. Filter first with :func:`filter_edges`.
    """
    return pd.DataFrame({
        'Source': edges['feature_a'],
        'Target': edges['feature_b'],
        'Weight': edges['r'].abs(),
        'Sign': np.where(edges['r'] >= 0, 'positive', 'negative'),
        'r': edges['r'],
        'q_value': edges['q_value'],
    }).reset_index(drop=True)
