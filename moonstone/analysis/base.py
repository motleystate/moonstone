import logging
from abc import ABC, abstractmethod

import pandas as pd

from moonstone.core.module_base import BaseDF
from moonstone.filtering.basics_filtering import NoCountsFiltering
from moonstone.transformers.mergers import MergeCountsAndMetadata

logger = logging.getLogger(__name__)


class BaseAnalysis(BaseDF, ABC):

    @abstractmethod
    def analyse(self):
        """
        Perform analysis.
        """
        pass


class MetadataBasedAnalysis(BaseAnalysis):

    def __init__(self, dataframe: pd.DataFrame, metadata: pd.DataFrame):
        super().__init__(dataframe)
        self.metadata_df = metadata
        filtering_instance = NoCountsFiltering(self.df)
        self.df = filtering_instance.filtered_df
        instance = MergeCountsAndMetadata(self.metadata_df, self.read_count_df)
        self.full_df = instance.full_df_with_features_in_columns
