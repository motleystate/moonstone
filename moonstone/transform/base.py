import logging

logger = logging.getLogger(__name__)


class TransformBase:

    def __init__(self, df):
        self.raw_df = df
        self.df = df
        self.history = []

    def historize(self, action, col_name, arguments=None):
        if arguments is None:
            arguments = {}
        self.history.append([action, {'col_name': col_name, **arguments}])

    def rename(self, col_name, new_name):
        self.df.rename(columns={col_name: new_name}, inplace=True)
        self.historize(self.rename.__name__, col_name, {'new_name': new_name})

    def run_transform(self, col_name, method_name, method_options):
        try:
            getattr(self, method_name)(col_name, **method_options)
        except AttributeError:
            logger.warning("%s is not a valid transformation name, transformation skipped.", method_name)
