import pandas as pd
from moonstone.parsers.base import BaseParser


class Metaphlan2Parser(BaseParser):
    """
    Parse output from metaphlan2 merged table
    """

    taxonomical_names = {'0': "kingdom", '1': "phylum", '2': "class",
                         '3': "order", '4': "family", '5': "genus", '6': "species"}
    taxa_column = 'ID'

    def split_taxa_columns_and_fill_none(self, df):
        """
        This function takes the block of taxa column and splits it into different ones.
        It also fill in None with latest found annotation
        """
        taxa_column = df[self.taxa_column].str.split("|", expand=True)
        taxa_in_lists = [taxa_column[i].str.split("_", expand=True) for i in range(taxa_column.shape[1])]
        taxa_df = pd.concat(taxa_in_lists, axis=1)
        return taxa_df.applymap(lambda x: x if not isinstance(x, str) else None if x.islower() else x)

    def naming_taxa_columns(self, df):
        """
        This function puts the name of the different column taxa according to the number that is stated in the
        taxa blobk by Qiime2. For example, it recognises the 0 as kingdom and add it to the column where the
        kingdoms are displayed.
        """
        column_names = [self.taxonomical_names[x] for x in self.taxonomical_names if x in df.loc[:, 1].values]
        taxa_df = df[3]
        taxa_df.columns = column_names
        return taxa_df

    def filling_missing_taxa_values(self, df):
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
        values_with_taxa_status = df.apply(lambda x: x + " ({})".format(x.name))
        filling_missing_values = values_with_taxa_status.fillna(method='ffill', axis=1)
        combining_df_and_taxa_status_in_nan = df.combine_first(filling_missing_values)
        return combining_df_and_taxa_status_in_nan

    def to_dataframe(self):
        df = super().to_dataframe()
        return df
    #         taxa_columns = self.spliting_into_taxa_columns(dataframe['OTU ID'])
    #         taxa_column_with_names = self.naming_taxa_columns(taxa_columns)
    #         complete_taxa_df = self.filling_missing_taxa_values(taxa_column_with_names)
    #         df_samples = dataframe.drop('OTU ID', axis=1)
    #         taxa_columns_and_df_samples = pd.concat([complete_taxa_df, df_samples], axis=1, sort=False)
    #         standard_taxa_df = taxa_columns_and_df_samples.set_index(list(taxa_column_with_names))
    #         setattr(self, "_dataframe", dataframe)
