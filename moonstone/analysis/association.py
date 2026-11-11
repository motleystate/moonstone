"""
Association between features (taxa, genes...) and metadata variables, adjusted for covariates.

A correlation between a taxon and a clinical variable cannot account for other variables (age, treatment...). This
module fits one linear model per feature instead, on the Centred Log-Ratio (CLR) transform of the table
(:func:`moonstone.transformers.compositional.clr`)::

    clr(feature) ~ variables + covariates

with Ordinary Least Squares (OLS). This is the model MaAsLin2 (Microbiome Multivariable Association with Linear
Models) fits with ``normalization = "NONE"``, ``transform = "NONE"``, ``analysis_method = "LM"`` on a table already
CLR-transformed, and no random effect, so results can be checked against it. MaAsLin2 standardises continuous
variables by default; here they stay in their own units unless ``standardize=True``.

The coefficient of a variable is the change in CLR value (the log of the feature's abundance relative to the
geometric mean of the sample) per unit of the variable, or between a level and the reference level for a
categorical variable. Covariates adjust the model but are not reported. p-values are corrected for multiple testing
over every feature and every reported variable level together (Benjamini-Hochberg by default).

Samples are matched to the metadata by name. Samples with no metadata row, or a missing value in any variable or
covariate of the model, are left out (complete cases); the number is reported.

Analyses are cross-sectional: one sample per subject. With ``group_column`` (for example a patient ID), the function
stops if a subject has more than one sample. A random effect per subject belongs with paired analyses, later.

:func:`ancombc_check` runs ANCOM-BC (Analysis of Compositions of Microbiomes with Bias Correction; Lin and Peddada,
2020, as implemented in scikit-bio) on the same model, as a sensitivity check: associations found by both methods
are more robust than those found by one.

From Python::

    from moonstone.analysis.association import feature_association

    results = feature_association(normalized_counts_df, metadata_df, variables='disease',
                                  covariates=['age', 'antibiotics'])
    significant = results[results['q_value'] < 0.05]
"""
import logging

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from moonstone.analysis.correlation import _feature_labels, check_one_sample_per_subject
from moonstone.transformers.compositional import _check_table, clr, multiplicative_replacement

logger = logging.getLogger(__name__)

RESULT_COLUMNS = ['feature', 'variable', 'level', 'reference', 'coef', 'std_error', 'ci_low', 'ci_high', 't_value',
                  'p_value', 'q_value', 'n_samples', 'n_dropped', 'prevalence']
ANCOMBC_COLUMNS = ['feature', 'variable', 'level', 'reference', 'log_fold_change', 'std_error', 'p_value', 'q_value',
                   'n_samples', 'n_dropped']
CORRECTION_METHODS = (None, 'fdr_bh', 'bonferroni')
CODED_GROUPS_MAX_LEVELS = 5


def _as_list(names) -> list:
    if names is None:
        return []
    if isinstance(names, str):
        return [names]
    return list(names)


def _is_categorical(column: pd.Series, forced: bool) -> bool:
    return forced or not pd.api.types.is_numeric_dtype(column) or pd.api.types.is_bool_dtype(column)


def _levels(column: pd.Series) -> list:
    """Observed levels, in category order for a pandas Categorical, else sorted."""
    observed = set(column.unique())
    if isinstance(column.dtype, pd.CategoricalDtype):
        return [level for level in column.cat.categories if level in observed]
    return sorted(observed, key=str)


def _check_model_names(metadata: pd.DataFrame, variables: list, covariates: list, categorical: list,
                       reference: dict) -> list:
    if not variables:
        raise ValueError("Give at least one variable to test.")
    model_columns = variables + covariates
    if len(set(model_columns)) != len(model_columns):
        raise ValueError("A variable is listed twice (in variables, covariates or both).")
    unknown = [c for c in model_columns if c not in metadata.columns]
    if unknown:
        raise ValueError(f"Not metadata columns: {unknown}.")
    for name in categorical + list(reference):
        if name not in model_columns:
            raise ValueError(f"{name!r} is given in categorical or reference but is not in the model.")
    return model_columns


def _complete_cases(metadata: pd.DataFrame, samples: list, model_columns: list) -> pd.DataFrame:
    """Metadata of the samples with a metadata row and no missing value in the model columns, matched by name."""
    meta = metadata.set_axis([str(i) for i in metadata.index], axis=0)
    if meta.index.duplicated().any():
        raise ValueError("Sample names are used more than once in the metadata.")
    no_row = [s for s in samples if s not in meta.index]
    if no_row:
        logger.warning("%d sample(s) have no metadata row and are left out, for example %s.", len(no_row), no_row[:5])
    meta = meta.loc[[s for s in samples if s in meta.index], model_columns]
    complete = meta.notna().all(axis=1)
    if not complete.all():
        logger.warning("%d sample(s) have a missing value in %s and are left out.", int((~complete).sum()),
                       model_columns)
    return meta[complete]


def _categorical_columns(name: str, column: pd.Series, reference: dict):
    """One indicator column per level other than the reference: {column name: (values, level)}, and the reference."""
    levels = _levels(column)
    if len(levels) < 2:
        raise ValueError(f"{name!r} has a single value in the samples analysed ({levels}): nothing to compare.")
    ref = reference.get(name, levels[0])
    if ref not in levels:
        matches = [level for level in levels if str(level) == str(ref)]
        if not matches:
            raise ValueError(f"Reference {ref!r} is not a value of {name!r} in the samples analysed: {levels}.")
        ref = matches[0]
    columns = {f'{name}[{level}]': ((column == level).to_numpy(dtype=float), str(level))
               for level in levels if level != ref}
    return columns, str(ref)


def _numeric_column(name: str, column: pd.Series, reference: dict, standardize: bool) -> np.ndarray:
    if name in reference:
        raise ValueError(f"{name!r} is numeric: a reference level only applies to categorical variables (add it to "
                         "categorical if its numbers are group codes).")
    values = column.to_numpy(dtype=float)
    unique = np.unique(values)
    if len(unique) < 2:
        raise ValueError(f"{name!r} has a single value in the samples analysed: nothing to compare.")
    if 2 < len(unique) <= CODED_GROUPS_MAX_LEVELS and np.all(unique == np.round(unique)):
        logger.warning("%r has only %d distinct whole-number values and is treated as a number (one slope). If they "
                       "are group codes, add it to categorical.", name, len(unique))
    if standardize:
        values = (values - values.mean()) / values.std(ddof=1)
    return values


def _design(metadata: pd.DataFrame, samples: list, variables: list, covariates: list, categorical: list,
            reference: dict, standardize: bool):
    """
    Design matrix (intercept first) for the samples with complete metadata.

    Returns:
        X (samples x columns), the reported terms as (column, variable, level, reference) tuples, the samples kept,
        and the number of samples dropped.
    """
    model_columns = _check_model_names(metadata, variables, covariates, categorical, reference)
    meta = _complete_cases(metadata, samples, model_columns)
    columns = {'Intercept': np.ones(len(meta))}
    terms = []
    for name in model_columns:
        if _is_categorical(meta[name], name in categorical):
            indicators, ref = _categorical_columns(name, meta[name], reference)
            for key, (values, level) in indicators.items():
                columns[key] = values
                terms.append((key, name, level, ref))
        else:
            columns[name] = _numeric_column(name, meta[name], reference, standardize)
            terms.append((name, name, None, None))
    design = pd.DataFrame(columns, index=meta.index)
    reported = [term for term in terms if term[1] in variables]
    return design, reported, list(meta.index), len(samples) - len(meta)


def _check_estimable(design: pd.DataFrame):
    n, p = design.shape
    if n <= p:
        raise ValueError(f"{n} samples are too few to estimate {p} coefficients: use fewer covariates or more samples.")
    if np.linalg.matrix_rank(design.to_numpy()) < p:
        raise ValueError("Some variables or levels are collinear (one is a combination of others, for example two "
                         "columns coding the same groups), so their effects cannot be separated. Remove one.")


def _ols(design: np.ndarray, response: np.ndarray):
    """OLS fit of every column of response (samples x features) on the same design (samples x coefficients)."""
    n, p = design.shape
    xtx_inv = np.linalg.inv(design.T @ design)
    coef = xtx_inv @ design.T @ response
    residuals = response - design @ coef
    df = n - p
    sigma2 = (residuals ** 2).sum(axis=0) / df
    with np.errstate(invalid='ignore', divide='ignore'):
        se = np.sqrt(np.outer(np.diag(xtx_inv), sigma2))
        t_value = coef / se
    # a perfect fit (constant feature, or residuals at rounding level) has no usable test
    perfect = sigma2 <= 1e-20 * np.maximum((response ** 2).mean(axis=0), 1e-300)
    se[:, perfect] = np.nan
    t_value[:, perfect] = np.nan
    return coef, se, t_value, df


def _correct(results: pd.DataFrame, correction_method: str):
    tested = results['p_value'].notna()
    if correction_method is not None and tested.any():
        results.loc[tested, 'q_value'] = multipletests(results.loc[tested, 'p_value'], method=correction_method)[1]


def _prepare(counts, metadata, variables, covariates, categorical, reference, group_column, correction_method,
             standardize):
    if correction_method not in CORRECTION_METHODS:
        raise ValueError(f"correction_method must be 'fdr_bh', 'bonferroni' or None, got {correction_method!r}.")
    values = _check_table(counts)
    samples = [str(s) for s in values.columns]
    if len(set(samples)) != len(samples):
        raise ValueError("Sample names are used more than once in counts.")
    values = values.set_axis(samples, axis=1)
    design, terms, kept, n_dropped = _design(
        metadata, samples, _as_list(variables), _as_list(covariates), _as_list(categorical), dict(reference or {}),
        standardize,
    )
    if group_column is not None:
        check_one_sample_per_subject(kept, metadata, group_column)
    _check_estimable(design)
    return values, design, terms, kept, n_dropped


def feature_association(
    counts: pd.DataFrame, metadata: pd.DataFrame, variables, covariates=None, categorical=None,
    reference: dict = None, standardize: bool = False, zeros: str = 'multiplicative', delta: float = None,
    pseudocount: float = None, correction_method: str = 'fdr_bh', group_column: str = None,
) -> pd.DataFrame:
    """
    Association of every feature with metadata variables, adjusted for covariates: one OLS model per feature on CLR
    values.

    Filter the table first (for example on prevalence): the CLR depends on which features are kept, and every extra
    feature adds tests to the correction.

    Args:
        counts: table with features as rows and samples as columns (counts, normalised counts or relative abundances).
        metadata: one row per sample, sample names as index. Samples are matched to counts by name.
        variables: metadata column(s) whose association is tested and reported.
        covariates: metadata column(s) the model adjusts for; they are not reported.
        categorical: numeric column(s) to treat as groups rather than numbers (for example 1/2/3 codes). Text and
          true/false columns are always categorical.
        reference: reference level per categorical column, for example ``{'disease': 'healthy'}``. By default, the
          first category of a pandas Categorical, otherwise the first value in sorted order.
        standardize: if True, continuous variables and covariates are scaled to mean 0 and standard deviation 1, so a
          coefficient is per standard deviation (MaAsLin2's default). False by default: per unit.
        zeros, delta, pseudocount: zero handling for the CLR, see :func:`moonstone.transformers.compositional.clr`.
        correction_method: 'fdr_bh' (Benjamini-Hochberg, default), 'bonferroni', or None for no correction.
        group_column: metadata column identifying the subject; if given, stop when a subject has more than one sample.

    Returns:
        One row per feature and reported variable level: feature, variable, level and reference (None for a
        continuous variable), coef, std_error, ci_low and ci_high (95% confidence interval), t_value, p_value,
        q_value (NaN when not corrected), n_samples (used in the model), n_dropped (left out for missing metadata)
        and prevalence (share of the samples used where the feature is non-zero). Statistics are NaN for a feature
        the model fits perfectly (for example a constant one).
    """
    values, design, terms, kept, n_dropped = _prepare(
        counts, metadata, variables, covariates, categorical, reference, group_column, correction_method, standardize,
    )
    # the CLR of a sample only depends on that sample, so transforming before dropping samples changes nothing
    transformed = clr(values, zeros=zeros, delta=delta, pseudocount=pseudocount)[kept]
    coef, se, t_value, df = _ols(design.to_numpy(), transformed.to_numpy().T)
    p_value = 2 * stats.t.sf(np.abs(t_value), df)
    margin = stats.t.ppf(0.975, df) * se

    labels = _feature_labels(values.index)
    prevalence = (values[kept] > 0).mean(axis=1).to_numpy()
    position = {column: i for i, column in enumerate(design.columns)}
    blocks = []
    for column, variable, level, ref in terms:
        i = position[column]
        blocks.append(pd.DataFrame({
            'feature': labels, 'variable': variable, 'level': level, 'reference': ref,
            'coef': coef[i], 'std_error': se[i], 'ci_low': coef[i] - margin[i], 'ci_high': coef[i] + margin[i],
            't_value': t_value[i], 'p_value': p_value[i], 'q_value': np.full(len(labels), np.nan),
            'n_samples': len(kept), 'n_dropped': n_dropped, 'prevalence': prevalence,
        }, columns=RESULT_COLUMNS))
    results = pd.concat(blocks, ignore_index=True)
    n_undefined = int(results['p_value'].isna().sum())
    if n_undefined:
        logger.warning("%d result(s) have no p-value because the model fits the feature perfectly (for example a "
                       "constant feature); they are not tested.", n_undefined)
    _correct(results, correction_method)
    return results


def ancombc_check(
    counts: pd.DataFrame, metadata: pd.DataFrame, variables, covariates=None, categorical=None,
    reference: dict = None, standardize: bool = False, correction_method: str = 'fdr_bh', group_column: str = None,
    max_iter: int = 100,
) -> pd.DataFrame:
    """
    ANCOM-BC on the same model as :func:`feature_association`, as a sensitivity check.

    ANCOM-BC (scikit-bio's implementation) needs strictly positive values, so zeros are replaced first by
    multiplicative replacement (:func:`moonstone.transformers.compositional.multiplicative_replacement`). Its
    log-fold changes are on the natural log scale. They are not on the same scale as the CLR coefficients of
    :func:`feature_association`, so compare the direction and significance of the results, not the values.

    Args:
        see :func:`feature_association`. max_iter: maximum iterations of the bias estimation.

    Returns:
        One row per feature and reported variable level: feature, variable, level, reference, log_fold_change,
        std_error, p_value, q_value (corrected over every feature and reported level together), n_samples, n_dropped.
    """
    from skbio.stats.composition import ancombc

    values, design, terms, kept, n_dropped = _prepare(
        counts, metadata, variables, covariates, categorical, reference, group_column, correction_method, standardize,
    )
    labels = _feature_labels(values.index)
    table = multiplicative_replacement(values[kept]).T.set_axis(labels, axis=1)
    # plain column names, so the formula does not depend on how variables and levels are spelled
    names = {column: f'x{i}' for i, column in enumerate(design.columns[1:])}
    model_metadata = design.drop(columns='Intercept').rename(columns=names)
    result = ancombc(table, model_metadata, ' + '.join(names.values()), max_iter=max_iter, p_adjust=None).result

    blocks = []
    for column, variable, level, ref in terms:
        rows = result.xs(names[column], level='Covariate').reindex(labels)
        blocks.append(pd.DataFrame({
            'feature': labels, 'variable': variable, 'level': level, 'reference': ref,
            'log_fold_change': rows['Log(FC)'].to_numpy(), 'std_error': rows['SE'].to_numpy(),
            'p_value': rows['pvalue'].to_numpy(), 'q_value': np.full(len(labels), np.nan),
            'n_samples': len(kept), 'n_dropped': n_dropped,
        }, columns=ANCOMBC_COLUMNS))
    results = pd.concat(blocks, ignore_index=True)
    _correct(results, correction_method)
    return results
