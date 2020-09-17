class BaseDownsizing:

    def __init__(self, raw_f, raw_r=None):
        self.raw_f = raw_f
        self.f = raw_f
        self._downsized_f = None
        if raw_r:
            self.raw_r = raw_r
            self.r = raw_r
            self._downsized_r = None

    def downsize_single(self):
        """Overridden in child classes to perform specified downsizing of fragment reads"""
        return self.f

    def downsize_pair(self):
        """Overridden in child classes to perform specified downsizing of paired-ends reads"""
        return self.f, self.r

    @property
    def downsized_pair(self):
        if getattr(self, "._downsized_f", None) is None:
            self._downsized_f, self_downsized_r = self.downsize_pair()
            self.f = self._downsized_f
            self.r = self._downsized_r
        return self._downsized_f, self._downsized_r

    @property
    def downsized_single(self):
        if getattr(self, "._downsized_f", None) is None:
            self._downsized_f = self.downsize_single()
            self.f = self._downsized_f
        return self._downsized_f
