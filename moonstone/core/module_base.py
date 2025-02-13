import logging
from abc import ABC

import pandas as pd

from moonstone.utils.plot import (
    add_default_titles_to_plotting_options,
    add_default_titles_to_plotting_options_3d,
    add_x_to_plotting_options
)

logger = logging.getLogger(__name__)


class BaseModule(ABC):
    """
    Base class for modules that needs both visualization and report_data
    """

    def _handle_plotting_options(
        self, plotting_options: dict,
        title: str, xlabel: str, ylabel: str, log_scale: bool,
        zlabel: str = None
    ):
        """
        Update or create plotting options with the title of the graph and of the different axis,
        as well as the log scale
        """
        if plotting_options is None:
            plotting_options = {}
        if log_scale:
            plotting_options = add_x_to_plotting_options(plotting_options, 'yaxes', 'type', "log")
            ylabel = f"{ylabel} (log)"
        if zlabel is not None:
            return add_default_titles_to_plotting_options_3d(
                plotting_options, title, xlabel, ylabel, zlabel
            )
        return add_default_titles_to_plotting_options(
            plotting_options, title, xlabel, ylabel
        )

    def visualize(self, write_file=False):
        """
        Generate visualization for the module
        """
        logger.warning("No visualization available for this module, please overload.")
        pass

    def generate_report_data(self) -> dict:
        """
        Overload this method to perform data reporting in child classes
        """
        logger.warning("No report data available for this module, please overload.")
        return {
            "title": "Report title",
            "data": {
                "var1": "parameters and/or",
                "var2": "results in numbers"
            }
        }

    @property
    def report_data(self) -> dict:
        if getattr(self, "_report_data", None) is None:
            self._report_data = self.generate_report_data()
        return self._report_data


class BaseDF(ABC):
    """
    Base class for modules that are based on a dataframe
    """

    def __init__(self, dataframe: pd.DataFrame):
        self.raw_df = dataframe
        self.df = dataframe
