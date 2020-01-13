import logging

logger = logging.getLogger(__name__)


class TransformBase:

    def __init__(self, df):
        self.raw_df = df
        self.df = df
        self.history = []

    def historize(self, action, arguments=None):
        if arguments is None:
            arguments = {}
        self.history.append([action, arguments])

    def run_transform(self, col_name, method_name, method_options):
        try:
            getattr(self, method_name)(col_name, **method_options)
        except AttributeError:
            logger.warning("%s is not a valid transformation name, transformation skipped.", method_name)
