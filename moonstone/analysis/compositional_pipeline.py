"""
Command-line steps for compositional correlation and association: reading a feature table, selecting samples,
filtering features, and writing a methods paragraph.

The order is fixed, because the Centred Log-Ratio (CLR) depends on which features are kept:

1. Read the table. If it has taxonomy columns (kingdom ... species), sum the counts to one level.
2. Keep the samples present in the metadata, and those matching every ``--select COLUMN=VALUE``.
3. Filter features on prevalence (and optionally keep the most abundant), on the selected samples.
4. Transform and test (:func:`moonstone.analysis.correlation.feature_correlation`,
   :func:`moonstone.analysis.association.feature_association`).

Each step only drops samples or features; values are never changed. Every choice goes into ``methods.txt``.
"""
import logging

import numpy as np
import pandas as pd

from moonstone.analysis.permanova import NON_SAMPLE_COLUMNS, TAXA_LEVELS, read_taxonomy_counts, sum_to_level
from moonstone.version import __version__

logger = logging.getLogger(__name__)

DEFAULT_MIN_PREVALENCE = 0.1
DEFAULT_CORRELATION_PERMUTATIONS = 999
METHOD_NAMES = {'spearman': 'Spearman', 'pearson': 'Pearson', 'kendall': 'Kendall',
                'rho': 'proportionality (rho)', 'sparcc': 'SparCC'}
CORRECTION_NAMES = {'fdr_bh': 'the Benjamini-Hochberg false discovery rate', 'bonferroni': 'Bonferroni'}


def read_feature_table(path: str, taxa_level: str = 'species') -> tuple:
    """
    Read a table with one row per feature and one column per sample (tab- or comma-separated).

    With taxonomy columns (kingdom ... species), the counts are summed to ``taxa_level``. Otherwise the first
    column holds the feature names (for example gene names).

    Returns:
        The table (features x samples) and a list of notes for the methods paragraph.
    """
    with open(path) as handle:
        header = handle.readline()
    sep = '\t' if '\t' in header else ','
    columns = [c.strip().strip('"') for c in header.rstrip('\n').split(sep)]
    if any(level in columns for level in TAXA_LEVELS):
        table, warnings = sum_to_level(read_taxonomy_counts(path), taxa_level)
        for warning in warnings:
            logger.warning("%s: %s.", path, warning)
        return table, [f"the counts were summed to {taxa_level} level"]
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str)
    df.columns = [str(c) for c in df.columns]
    df = df[[c for c in df.columns if c not in NON_SAMPLE_COLUMNS]]
    non_numeric = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if non_numeric:
        raise ValueError(f"{path}: sample columns must be numeric, but these are not: {', '.join(non_numeric[:5])}. "
                         "The first column must hold the feature names.")
    if df.isna().any().any():
        raise ValueError(f"{path}: {int(df.isna().sum().sum())} blank cell(s). Fill them in (with 0 if the feature "
                         "was not found) before running.")
    if df.index.duplicated().any():
        raise ValueError(f"{path}: feature names are used more than once, for example "
                         f"{df.index[df.index.duplicated()].unique()[:5].tolist()}.")
    return df, []


def parse_selection(text: str) -> tuple:
    """'Type=inclusion' -> ('Type', 'inclusion')."""
    column, separator, value = text.partition('=')
    if not separator or not column.strip() or not value.strip():
        raise ValueError(f"--select needs COLUMN=VALUE, got {text!r}.")
    return column.strip(), value.strip()


def select_samples(metadata: pd.DataFrame, selections: list) -> tuple:
    """
    Keep the metadata rows matching every (column, value) pair. Values are compared as text.

    Returns:
        The selected metadata and a list of notes, one per selection, with the number of samples kept.
    """
    notes = []
    for column, value in selections:
        if column not in metadata.columns:
            raise ValueError(f"--select: {column!r} is not a metadata column.")
        as_text = metadata[column].map(lambda v: None if pd.isna(v) else str(v).strip())
        keep = as_text == value
        if not keep.any():
            found = sorted(as_text.dropna().unique().tolist())
            raise ValueError(f"--select: no sample has {column} = {value!r}. Values found: {found[:20]}.")
        notes.append(f"{column} = {value} ({int(keep.sum())} of {len(metadata)} samples)")
        metadata = metadata[keep]
    return metadata, notes


def filter_features(table: pd.DataFrame, min_prevalence: float = DEFAULT_MIN_PREVALENCE, top: int = None) -> tuple:
    """
    Keep the features present (non-zero) in at least ``min_prevalence`` of the samples, then optionally the ``top``
    ones with the highest mean relative abundance. The filter never looks at the metadata.

    Returns:
        The filtered table and a note for the methods paragraph.
    """
    if not 0 <= min_prevalence <= 1:
        raise ValueError(f"min_prevalence must be between 0 and 1, got {min_prevalence}.")
    if top is not None and top < 2:
        raise ValueError(f"top must be at least 2, got {top}.")
    total = table.shape[0]
    filtered = table[(table > 0).mean(axis=1) >= min_prevalence]
    note = f"features present in at least {min_prevalence:.0%} of samples were kept"
    if top is not None and filtered.shape[0] > top:
        sums = filtered.sum(axis=0).replace(0, np.nan)
        mean_relative = (filtered / sums).mean(axis=1)
        filtered = filtered.loc[mean_relative.sort_values(ascending=False).index[:top]]
        note += f", then the {top} with the highest mean relative abundance"
    if filtered.shape[0] < 2:
        raise ValueError(f"Only {filtered.shape[0]} feature(s) left after filtering ({note}): lower --min_prevalence.")
    return filtered, f"{note} ({filtered.shape[0]} of {total})"


def _zeros_sentence() -> str:
    return ("Zeros were replaced by multiplicative replacement (Martin-Fernandez et al., 2003), each zero becoming 65% "
            "of the smallest non-zero proportion of its sample, and values were transformed with the centred "
            "log-ratio (CLR)")


def _table_sentence(name: str, notes: list, filter_note: str) -> str:
    return f"For the {name}, " + ", then ".join(notes + [filter_note]) + ", before the transform."


def _samples_sentence(n_samples: int, selection_notes: list, group_column: str) -> str:
    sentence = f"The analysis used {n_samples} samples"
    if selection_notes:
        sentence += " selected on " + " and ".join(selection_notes)
    if group_column:
        sentence += f", one per subject (checked on {group_column})"
    return sentence + "."


def correlation_methods(method: str, n_samples: int, tables: list, selection_notes: list, permutations: int,
                        seed: int, correction_method: str, alpha: float, group_column: str = None) -> str:
    """
    Methods paragraph for a correlation run.

    Args:
        tables: one (name, notes, filter_note) tuple per table.
    """
    parts = [f"Feature correlations were computed with Moonstone {__version__}.",
             _samples_sentence(n_samples, selection_notes, group_column)]
    parts += [_table_sentence(*table) for table in tables]
    if method == 'sparcc':
        parts.append("Correlations were estimated with SparCC (Friedman and Alm, 2012), as the median of 20 Dirichlet "
                     "draws of the sample fractions, with 10 rounds of excluding pairs above 0.1.")
    else:
        parts.append(_zeros_sentence() + (", each table on its own." if len(tables) > 1 else "."))
        parts.append(f"{METHOD_NAMES[method]} coefficients were computed for every pair of features"
                     + (" across the two tables." if len(tables) > 1 else "."))
    parts.append(f"p-values were obtained from a permutation null: {permutations} permutations (seed {seed}), each "
                 "shuffling every feature across samples independently and re-transforming the table, the observed "
                 "coefficient being compared with the mean and standard deviation of the permuted ones.")
    parts.append(_correction_sentence(correction_method, alpha, 'pairs'))
    return ' '.join(parts)


def association_methods(variables: list, covariates: list, n_samples: int, n_dropped: int, table: tuple,
                        selection_notes: list, correction_method: str, alpha: float, standardize: bool = False,
                        group_column: str = None, ancombc: bool = False, references: dict = None) -> str:
    """
    Methods paragraph for an association run.

    Args:
        table: (name, notes, filter_note) of the count table.
        references: reference level of each reported categorical variable.
    """
    parts = [f"Associations between features and metadata were tested with Moonstone {__version__}.",
             _samples_sentence(n_samples, selection_notes, group_column)]
    if n_dropped:
        parts.append(f"{n_dropped} further sample(s) with missing metadata were left out (complete cases).")
    parts.append(_table_sentence(*table))
    parts.append(_zeros_sentence() + ".")
    model = ' + '.join(variables + covariates)
    parts.append(f"For each feature, an ordinary least squares model, CLR ~ {model}, was fitted"
                 + (", with continuous variables scaled to unit standard deviation." if standardize else ".")
                 + (f" {', '.join(covariates)} {'was' if len(covariates) == 1 else 'were'} included as "
                    f"{'a covariate' if len(covariates) == 1 else 'covariates'}." if covariates else ""))
    if references:
        parts.append("Categorical variables were compared with a reference level: "
                     + ", ".join(f"{variable} = {level}" for variable, level in references.items()) + ".")
    parts.append(_correction_sentence(correction_method, alpha, 'tests of every feature and variable level'))
    if ancombc:
        parts.append("As a sensitivity check, ANCOM-BC (Lin and Peddada, 2020; scikit-bio implementation) was run on "
                     "the same model after multiplicative replacement of zeros, with the same correction.")
    return ' '.join(parts)


def _correction_sentence(correction_method: str, alpha: float, what: str) -> str:
    if correction_method is None:
        return f"p-values were not corrected for multiple testing; results with p < {alpha} are reported."
    return (f"p-values were corrected for multiple testing over all {what} with {CORRECTION_NAMES[correction_method]} "
            f"method, and results with q < {alpha} are reported.")
