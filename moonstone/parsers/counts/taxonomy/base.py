from moonstone.parsers.base import BaseParser


class TaxonomyCountsBaseParser(BaseParser):

    def _fill_none(self, taxa_df):
        """
        This function serves to obtain a data frame that fills the None values with the last valid value and
        the category where it belonged to. E.g:

        Before:
        column names     kingdom  phylum            family         genus
        value            Bacteria Bacteroidetes ... Tannerellaceae None

        After:
        column names     kingdom  phylum            family         genus
        value            Bacteria Bacteroidetes ... Tannerellaceae Tannerellaceae (family)
        """
        taxa_df_with_rank = taxa_df.apply(lambda x: x + " ({})".format(x.name))
        taxa_df_with_rank_filled_none = taxa_df_with_rank.fillna(method='ffill', axis=1)
        taxa_df_filled_none = taxa_df.combine_first(taxa_df_with_rank_filled_none)
        return taxa_df_filled_none
