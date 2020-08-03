from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base import TaxonomyCountsBaseParser


class SunbeamKraken2Parser(TaxonomyCountsBaseParser):
    """
    Parse output from Kraken2 merge table from Sunbeam pipeline
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
        df = df.set_index(self.taxonomical_names[:self._rank_level])
        return df
