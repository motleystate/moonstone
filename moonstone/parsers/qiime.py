import pandas as pd
from .base import BaseParser


class Qiime2Parser(BaseParser):
    """
    This class serves for transforming the data obtained by Qiime2 (https://qiime2.org/)
    into a dataframe that has the taxa information ordered in a multilevel index and the
    readings in the correspondent sample values. In other words, it can be used, for
    example, before the DESeq2_normalization module to organize the data properly before
    is transformed.
    Example:
    - Qiime2 gives this info: D_0__Bacteria;D_1__Bacteroidetes;...;D_4__Tannerellaceae;D_5__Macellibacteroides
    But with the module we obtain:

    column names     kingdom  phylum            family         genus
    value            Bacteria Bacteroidetes ... Tannerellaceae Macellibacteroides
    """

    taxonomical_names = {'0': "kingdom", '1': "phylum", '2': "class",
                         '3': "order", '4': "family", '5': "genus", '6': "species"}
    symbol_to_remove = '#'
    terms_to_remove = ["Ambiguous_taxa", "Unknown Family"]

    def __init__(self, file_path, skip_row_number=1):
        super().__init__(file_path)
        self.skip_row_number = skip_row_number

    def to_dataframe(self):
        """
        Importing the csv and placing the samples and block of taxa as headers of the different columns
        """
        df = pd.read_csv(self.file_path, skiprows=self.skip_row_number)
        new_column_names = [column_name.replace(self.symbol_to_remove, "") for column_name in list(df)]
        df.columns = new_column_names
        self._dataframe = df

    def spliting_into_taxa_columns(self, serie):
        """
        This function takes the block of taxa column and splits it into different ones. In addition,
        it serves to clean the data a little bit since it sets to none all ambiguities
        (like 'uncultured') we might have in the data.
        """
        taxa_column = serie.str.split(";", expand=True)
        taxa_column = taxa_column.replace(self.terms_to_remove, None)
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

    @property
    def standard_taxa_df(self):
        """
        This property can be used to obtain the formated data frame directly, without any intermidiate steps.
        """
        if getattr(self, "_standard_taxa_df", None) is None:
            taxa_columns = self.spliting_into_taxa_columns(self.dataframe['OTU ID'])
            taxa_column_with_names = self.naming_taxa_columns(taxa_columns)
            complete_taxa_df = self.filling_missing_taxa_values(taxa_column_with_names)
            df_samples = self._dataframe.drop('OTU ID', axis=1)
            taxa_columns_and_df_samples = pd.concat([complete_taxa_df, df_samples], axis=1, sort=False)
            standard_taxa_df = taxa_columns_and_df_samples.set_index(list(taxa_column_with_names))
            setattr(self, "_standard_taxa_df", standard_taxa_df.astype(int))
        return self._standard_taxa_df