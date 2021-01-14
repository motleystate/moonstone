from pandas import DataFrame

from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class SunbeamKraken2Parser(TaxonomyCountsBase, BaseParser):
    """
    Parse output from `Kraken2 <https://ccb.jhu.edu/software/kraken2/>`_
    merge table from `Sunbeam <https://github.com/sunbeam-labs/sunbeam/>`_ pipeline.
    """

    taxa_column = 'Consensus Lineage'
    new_otu_id_name = 'NCBI_taxonomy_ID'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def to_dataframe(self) -> DataFrame:
        df = super().to_dataframe()
        # Rename first column to NCBI_taxonomy_ID
        df.columns = [self.new_otu_id_name] + list(df.columns[1:])
        df = self.split_taxa_fill_none(df, sep="; ", merge_genus_species=True)
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
