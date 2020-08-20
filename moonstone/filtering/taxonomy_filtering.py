from typing import List

import pandas as pd

from moonstone.filtering.base import BaseFiltering


class TaxonomyNamesFiltering(BaseFiltering):
    """
    Filtering a Taxonomy multiindexed dataframe on index names at a chosen level.
    """

    def __init__(
        self, dataframe: pd.DataFrame, names: List[str],
        axis: int = 1, level: str = 'species', keep: bool = True
    ):
        """
        :param names: list of index names
        :param level: level of the MultiIndex to filter on
        :param keep: keep column (discard them if set to False)
        """
        super().__init__(dataframe)
        self.names = names
        self.keep = keep
        self.level = level
        self._validate_parameters()

    def _validate_parameters(self):
        pass

    def filter(self) -> pd.DataFrame:
        corresponding_rows = self.df.index.get_level_values(self.level).isin(self.names)
        if self.keep:
            return self.df.loc[corresponding_rows]
        else:
            return self.df.loc[~corresponding_rows]
