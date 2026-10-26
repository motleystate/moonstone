"""
Log-ratio transforms for compositional count tables.

Relative abundances are compositional: the values of one sample are tied to a fixed total, so one feature going up
pushes the others down, and correlations computed on them are partly artefacts. The Centred Log-Ratio (CLR) transform
removes that constraint by taking each value relative to the geometric mean of its sample:

    clr(x_i) = log(x_i) - mean(log(x))

Tables follow Moonstone's layout: one row per feature (taxon, gene...), one column per sample. The transform is
applied to each sample (column) and returns a new table of the same shape; the input is never changed.

CLR needs strictly positive values, so zeros must be handled first, and how they are handled is always an explicit
choice:

- ``zeros="multiplicative"`` (default): multiplicative replacement (Martin-Fernandez et al., 2003). In each sample,
  zeros become a small proportion ``delta`` and the non-zero proportions are shrunk so the sample still sums to 1.
  By default ``delta`` is 65% of the smallest non-zero proportion of that sample, so it stays below every observed
  value whatever the sequencing depth or normalisation.
- ``zeros="pseudocount"``: add ``pseudocount`` to every value. There is no default value: a pseudocount of 1 is
  meaningless on normalised data whose values are below 1.

The robust CLR (:func:`rclr`) uses only the non-zero values of each sample and leaves zeros missing (NaN).

The result of CLR depends on which features are in the table, because the geometric mean is taken over them.
Filter the table first (see :mod:`moonstone.filtering`), then transform.

From Python::

    from moonstone.transformers.compositional import clr

    clr_df = clr(normalized_counts_df)
"""
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

ZERO_METHODS = ('multiplicative', 'pseudocount')
# Martin-Fernandez et al. (2003) found 65% of the detection limit gave the least distortion.
DETECTION_FRACTION = 0.65


def _check_table(df: pd.DataFrame) -> pd.DataFrame:
    """Raise on input the transforms cannot handle; return the table as floats. Never fills or recodes values."""
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"Expected a pandas DataFrame (features as rows, samples as columns), got {type(df).__name__}.")
    if df.empty:
        raise ValueError("The count table is empty.")
    non_numeric = [col for col in df.columns if not pd.api.types.is_numeric_dtype(df[col])]
    if non_numeric:
        raise ValueError(
            f"{len(non_numeric)} column(s) are not numeric, for example {non_numeric[:3]}. Keep only sample columns "
            "(put taxonomy in the index)."
        )
    values = df.astype(float)
    if values.isna().any().any():
        missing = values.columns[values.isna().any()].tolist()
        raise ValueError(
            f"The table has missing values in {len(missing)} sample(s), for example {missing[:3]}. Fix the input "
            "first; the transform does not fill them in."
        )
    if (values < 0).any().any():
        raise ValueError("The table has negative values; log-ratio transforms need counts or abundances >= 0.")
    empty_samples = values.columns[(values == 0).all()].tolist()
    if empty_samples:
        raise ValueError(
            f"{len(empty_samples)} sample(s) have only zeros, for example {empty_samples[:3]}. Their log-ratios are "
            "undefined; remove them first."
        )
    n_empty_features = int((values == 0).all(axis=1).sum())
    if n_empty_features:
        logger.warning(
            "%d feature(s) are zero in every sample. They carry no information and change the CLR of the others; "
            "filter them out before the transform.", n_empty_features
        )
    return values


def multiplicative_replacement(df: pd.DataFrame, delta: float = None) -> pd.DataFrame:
    """
    Replace zeros by multiplicative replacement and return proportions (each sample sums to 1).

    Args:
        df: table with features as rows and samples as columns, values >= 0.
        delta: proportion given to each zero. None (default) uses, for each sample, 65% of its smallest non-zero
          proportion. A float applies the same proportion to every sample.

    Returns:
        A new table of proportions, same shape, with no zeros.
    """
    values = _check_table(df)
    if delta is not None and not 0 < delta < 1:
        raise ValueError(f"delta must be a proportion between 0 and 1, got {delta}.")
    return pd.DataFrame(_replace_zeros(values.to_numpy(), delta), index=values.index, columns=values.columns)


def _replace_zeros(values: np.ndarray, delta: float = None) -> np.ndarray:
    """Multiplicative replacement on an array (features x samples) already checked by _check_table."""
    proportions = values / values.sum(axis=0)
    zeros = proportions == 0
    if delta is None:
        with np.errstate(invalid='ignore'):
            delta_per_sample = np.nanmin(np.where(zeros, np.nan, proportions), axis=0) * DETECTION_FRACTION
    else:
        delta_per_sample = np.full(proportions.shape[1], float(delta))
    zero_mass = zeros.sum(axis=0) * delta_per_sample
    if (zero_mass >= 1).any():
        raise ValueError("delta is too large: the replaced zeros would take up the whole sample. Use a smaller delta.")
    return np.where(zeros, delta_per_sample, proportions * (1 - zero_mass))


def _clr_array(values: np.ndarray, zeros: str = 'multiplicative', delta: float = None,
               pseudocount: float = None) -> np.ndarray:
    """CLR of an array (features x samples) already checked, with arguments already validated by clr()."""
    positive = _replace_zeros(values, delta) if zeros == 'multiplicative' else values + pseudocount
    logs = np.log(positive)
    return logs - logs.mean(axis=0)


def _centre_logs(positive: pd.DataFrame) -> pd.DataFrame:
    logs = np.log(positive)
    return logs - logs.mean(axis=0)


def clr(df: pd.DataFrame, zeros: str = 'multiplicative', delta: float = None, pseudocount: float = None) -> pd.DataFrame:
    """
    Centred log-ratio transform of each sample (column).

    Args:
        df: table with features as rows and samples as columns, values >= 0 (counts, normalised counts or
          relative abundances: the result is the same for any per-sample scaling).
        zeros: 'multiplicative' (default) or 'pseudocount'. See the module documentation.
        delta: only with zeros='multiplicative'; see :func:`multiplicative_replacement`.
        pseudocount: only with zeros='pseudocount'; required, must be > 0.

    Returns:
        A new table of the same shape; each column sums to 0.
    """
    if zeros not in ZERO_METHODS:
        raise ValueError(f"zeros must be one of {ZERO_METHODS}, got {zeros!r}.")
    if zeros == 'multiplicative':
        if pseudocount is not None:
            raise ValueError("pseudocount is only used with zeros='pseudocount'.")
        if delta is not None and not 0 < delta < 1:
            raise ValueError(f"delta must be a proportion between 0 and 1, got {delta}.")
    else:
        if delta is not None:
            raise ValueError("delta is only used with zeros='multiplicative'.")
        if pseudocount is None or not pseudocount > 0:
            raise ValueError("zeros='pseudocount' needs an explicit pseudocount > 0.")
    values = _check_table(df)
    return pd.DataFrame(_clr_array(values.to_numpy(), zeros, delta, pseudocount),
                        index=values.index, columns=values.columns)


def rclr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Robust centred log-ratio transform: CLR computed on the non-zero values of each sample only.

    Zeros are not replaced; they stay missing (NaN) in the result, so later steps must handle missing values.

    Args:
        df: table with features as rows and samples as columns, values >= 0.

    Returns:
        A new table of the same shape, NaN where the input was 0.
    """
    values = _check_table(df)
    return _centre_logs(values.where(values > 0))
