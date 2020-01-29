class BaseNormalization:

    def __init__(self, df):
        self.raw_df = df
        self.df = df

    def normalize(self):
        """
        Overload this method to perform normalization in child classes
        """
        return self.raw_df

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            self._normalized_df = self.normalize()
            self.df = self._normalized_df
        return self._normalized_df
