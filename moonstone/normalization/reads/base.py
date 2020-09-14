class BaseStandardization:

    def __init__(self, raw_f, raw_r=None):
        self.raw_f = raw_f
        if raw_r:
            self.raw_r = raw_r

    def downsize(self):
