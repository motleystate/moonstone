import logging

logger = logging.getLogger(__name__)


class BaseModule:

    def visualize(self, write_file=False):
        """
        Generate visualization for the module
        """
        logger.warning("No visualization available for this module, please overload.")
        pass

    def generate_report_data(self):
        """
        Overload this method to perform normalization in child classes
        """
        logger.warning("No report data available for this module, please overload.")
        return {
            "title": "Report title",
            "text": "Some text about the module"
        }

    @property
    def report_data(self):
        if getattr(self, "_report_data", None) is None:
            self._report_data = self.generate_report_data()
        return self._report_data
