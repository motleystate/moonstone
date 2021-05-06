from pandas import DataFrame

from moonstone.parsers.counts.taxonomy.base import BaseTaxonomyCountsParser


class SunbeamKraken2Parser(BaseTaxonomyCountsParser):
    """
    Parse output from `Kraken2 <https://ccb.jhu.edu/software/kraken2/>`_
    merge table from `Sunbeam <https://github.com/sunbeam-labs/sunbeam/>`_ pipeline.
    """

    taxa_column = 'Consensus Lineage'
    new_otu_id_name = 'NCBI_taxonomy_ID'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, parsing_options={'skiprows': 1}, **kwargs)

    def _load_data(self) -> DataFrame:
        df = super()._load_data()
        # Rename first column to NCBI_taxonomy_ID
        df.columns = [self.new_otu_id_name] + list(df.columns[1:])
        df = self.split_taxa_fill_none(df, sep="; ", merge_genus_species=True)
        df = df.set_index(self.taxonomical_names[:self.rank_level])
        return df
