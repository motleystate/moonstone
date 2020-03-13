class BaseScaling:

    def __init__(self, raw_x):
        self.raw_x = raw_x
        self.x = raw_x
        self._scaled_x = None

    def scaler(self):
        """Overridden in child classes to perform a specific scaling."""
        return self.x

    @property
    def scale_x(self):
        if getattr(self, "_scaled_x", None) is None:
            self._scaled_x = self.scaler()
            self.x = self._scaled_x
        return self._scaled_x
