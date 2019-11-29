import pandas as pd
import logging
module_logger = logging.getLogger(__name__)


class BaseParser(object):

    def __init__(self, file_path):
        self.logger = module_logger
        self.logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.file_path = file_path

    @property
    def dataframe(self):
        if getattr(self, '_dataframe', None) is None:
            self.to_dataframe()
        return self._dataframe

    def to_dataframe(self):
        self._dataframe = pd.read_csv(self.file_path)
