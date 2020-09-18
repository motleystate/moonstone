class BaseDownsizing:

    def __init__(self, raw_file_f, raw_file_r=None):
        self.raw_file_f = raw_file_f
        self.raw_file_f = raw_file_f
        self._downsized_f = None
        if raw_file_r:
            self.raw_file_r = raw_file_r
            self.raw_file_r = raw_file_r
            self._downsized_r = None

    def downsize_single(self):
        """Overridden in child classes to perform specified downsizing of fragment reads"""
        return self.raw_file_f

    def downsize_pair(self):
        """Overridden in child classes to perform specified downsizing of paired-ends reads"""
        return self.raw_file_f, self.raw_file_r

    @property
    def downsized_pair(self):
        if getattr(self, "._downsized_f", None) is None:
            self._downsized_f, self_downsized_r = self.downsize_pair()
            self.raw_file_f = self._downsized_f
            self.raw_file_r = self._downsized_r
        return self._downsized_f, self._downsized_r

    @property
    def downsized_single(self):
        if getattr(self, "._downsized_f", None) is None:
            self._downsized_f = self.downsize_single()
            self.raw_file_f = self._downsized_f
        return self._downsized_f
