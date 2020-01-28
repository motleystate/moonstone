class BaseNormalization:

    def __init__(self, df):
        self.df = df

    def normalize(self):
        """
        Overload this method to perform normalization in child classes
        """
        return self.df

    @property
    def normalized_df(self):
        if getattr(self, "_normalized_df", None) is None:
            setattr(self, "_normalized_df", self.normalize())
        return self._normalized_df
