from pandas import DataFrame

from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class BaseTaxonomyCountsParser(TaxonomyCountsBase, BaseParser):
    pass
