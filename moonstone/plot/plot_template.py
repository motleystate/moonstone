import math
from matplotlib import pyplot
import numpy as np
import pandas as pd
import pylab as plt


class Histogram:

    def __init__(self, dataframe, plotting_options, **kwargs):
        self.df = dataframe
        self.plotting_options = plotting_options
        if 'output_file' in kwargs.keys() and kwargs['output_file']:
            self.output_file = kwargs['output_file']
            if self.output_file is True:
                # no name given for the output file so a generic name will be given
                self.output_file = "histogram.html"
        else:
            self.output_file = False

    def compute_asymetric_bins(self):
        """Logish bins"""
        maximum = self.df.max()
        magnitude = int(math.log10(maximum))
        binsvalues = [-0.1]                       # to have the 0 value
        i = 0
        while i < magnitude:
            binsvalues += list(np.arange(2*10**i, 10**(i+1)+1, 10**i))
            i += 1
        # i=magnitude
        binsvalues += list(np.arange(2*10**i, maximum+10**i, 10**i))        # up to maximum
        return binsvalues

    def in_bins(self, bins_values):
        self.df = pd.cut(self.df, bins=bins_values)   # put every items in the appropriate bin

    def count(self, normalize=False):
        counts = pd.value_counts(self.df, normalize=normalize)
        self.counts = counts.reindex(self.df.cat.categories)

    def one_hist(self, title, xlabel, ylabel):
        xvalues = np.arange(1, (len(self.counts))*1.5, 1.5)
        xnames = list(self.counts.index)

        fig, ax = plt.subplots(figsize=(20, 15))

        if 'log' in self.plotting_options.keys() and self.plotting_options['log']:
            pyplot.yscale('log')
            pyplot.grid(which='minor', axis='both')
            ax.set_axisbelow(True)  # grid lines are behind the rest

        x = plt.bar(xvalues, self.counts, 1)   # noqa

        plt.xticks(xvalues, xnames)
        plt.xticks(rotation=90)
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.title(title, fontsize=30)
        plt.xlabel(xlabel, fontsize=24, labelpad=20)
        plt.ylabel(ylabel, fontsize=24, labelpad=20)
        plt.show()

        if self.output_file:
            fig.savefig(self.output_file)
