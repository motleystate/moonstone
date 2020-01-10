class TransformBase:

    def __init__(self, df):
        self.raw_df = df
        self.df = df
        self.history = []

    def historize(self, action, arguments=None):
        if arguments is None:
            arguments = {}
        self.history.append([action, arguments])
