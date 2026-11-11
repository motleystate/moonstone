
import pandas as pd
import numpy as np

from moonstone.normalization.base import BaseNormalization


class RPKMNormalization(BaseNormalization):
    """
    Normalize feature counts using Reads Per Kilobase per Million
    mapped reads (RPKM).

    Feature lengths must be provided in base pairs (bp)
    Sequencing depth is estimated from the sum of
    counts per sample unless custom depths are provided ().

    Note:
        Filtering features before automatic normalization may underestimate
        sequencing depth and bias RPKM values 

    Reference: https://doi.org/10.1038/nmeth.1226
    """

    def __init__(
        self,
        dataframe: pd.DataFrame,
        feature_lengths: pd.Series,
        sequencing_depth: pd.Series | None = None,
    ):
        """
        :param dataframe: Count matrix with features as rows and samples as columns
        :param feature_lengths: Feature lengths in base pairs (bp)
        :param sequencing_depth: sample-specific sequencing depths (optional)
        """
        super().__init__(dataframe)
        self.feature_lengths = feature_lengths
        self.sequencing_depth = sequencing_depth

    def normalize(self):
        """
        Compute RPKM-normalized counts.
        """
        # Validate and convert counts
        counts = self.raw_df.apply(pd.to_numeric, errors="coerce")

        if counts.isna().any().any():
            raise ValueError("Counts contain missing or non-numeric values.")

        if not np.isfinite(counts.to_numpy()).all():
            raise ValueError("Counts must be finite.")

        if (counts < 0).any().any():
            raise ValueError("Counts must be non-negative.")

        # Compute or retrieve sequencing depth
        if self.sequencing_depth is None:
            sequencing_depth = counts.sum(axis=0)
        else:
            sequencing_depth = self.sequencing_depth.reindex(counts.columns)

        sequencing_depth = pd.to_numeric(
            sequencing_depth, errors="coerce"
        )

        if sequencing_depth.isna().any():
            missing = sequencing_depth.index[sequencing_depth.isna()].tolist()
            raise ValueError(f"Missing or invalid sequencing depths for: {missing}")

        if not np.isfinite(sequencing_depth).all():
            raise ValueError("Sequencing depths must be finite.")

        if (sequencing_depth <= 0).any():
            raise ValueError("Sequencing depths must be strictly positive.")

        # Validate and convert feature lengths
        feature_lengths = pd.to_numeric(
            self.feature_lengths.reindex(counts.index),
            errors="coerce",
        )

        if feature_lengths.isna().any():
            missing = feature_lengths.index[feature_lengths.isna()].tolist()
            raise ValueError(f"Missing or invalid feature lengths for: {missing}")

        if not np.isfinite(feature_lengths).all():
            raise ValueError("Feature lengths must be finite.")

        if (feature_lengths <= 0).any():
            raise ValueError("Feature lengths must be strictly positive.")

        # Compute RPKM
        rpk = counts.div(feature_lengths / 1e3, axis=0)
        rpkm = rpk.div(sequencing_depth / 1e6, axis=1)

        return rpkm

