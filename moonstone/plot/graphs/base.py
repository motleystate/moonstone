from abc import ABC, abstractmethod
from typing import Union

import pandas as pd


class BaseGraph(ABC):

    def __init__(self, data: Union[pd.Series, pd.DataFrame], plotting_options: dict = None,
                 show: bool = True, output_file: Union[bool, str] = False):
        """
        :param data: data to plot
        :param show: set to False if you don't want to show the plot
        :param output_file: name of the output file
        :param plotting_options: options of plotting that will override the default setup \n
                                 [!] Make sure the value given to an argument is of the right type \n
                                 options allowed : 'log': `bool` ; 'colorbar': `[str, List[str]]` ;
                                 'tickangle': `[int, float]`
        """
        self.data = data
        self.plotting_options = plotting_options

        if type(show) == bool:
            self.show = show
        else:
            raise ValueError('Error : show value must be a bool, %s given' % type(show).__name__)

        self.output_file = output_file
        if self.output_file is True:
            # if no name given for the output file, a generic name is generated
            if data.name is not None:
                self.output_file = data.name+"_"+self.__class__.__name__+".html"
            else:
                self.output_file = self.__class__.__name__+".html"

    @abstractmethod
    def plot_one_graph(self, title: str, xlabel: str, ylabel: str):
        """
        method that plots the graph
        needs to be defined in every child class
        """
        pass
