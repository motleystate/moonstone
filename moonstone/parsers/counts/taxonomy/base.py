from moonstone.parsers.base import BaseParser
from moonstone.plot import PlotTaxonomyCounts
from moonstone.utils.taxonomy import TaxonomyCountsBase


class BaseTaxonomyCountsParser(TaxonomyCountsBase, BaseParser):
    PLOT_CLASS = PlotTaxonomyCounts
