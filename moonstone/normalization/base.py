from abc import ABC, abstractmethod

from moonstone.core.module_base import BaseModule, BaseDF


class BaseNormalization(BaseModule, BaseDF, ABC):

    @abstractmethod
    def normalize(self):
        """
        Method to perform normalization
        """
        pass

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            self._normalized_df = self.normalize()
            self.df = self._normalized_df
        return self._normalized_df
