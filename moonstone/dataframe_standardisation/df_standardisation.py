import pandas as pd


class Df_Standardisation:
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

    def __init__(self, filepath):
        self.filepath = filepath

    def import_df(self, filepath):
        """
        Importing the csv and placing the samples and block of taxa as headers of the different columns
        """
        df = pd.read_csv(self.filepath)
        return df.rename(columns=df.loc[0]).drop(df.index[0])

    def delete_lowercase(self, string):
        if string is not None:
            if string.islower():
                return None
            else:
                return string

    def spliting_into_taxa_columns(self, df):
        """
        This function takes the block of taxa column and splits it into different ones. In addition,
        it serves to clean the data a little bit since it sets to none all ambiguities
        (like 'uncultured') we might have in the data.
        """
        taxa_column = df['#OTU ID'].str.split(";", expand=True)
        taxa_column = taxa_column.replace("Ambiguous_taxa", "nothing")
        taxa_in_lists = [taxa_column[i].str.split("_", expand=True) for i in range(taxa_column.shape[1])]
        taxa_df = pd.concat(taxa_in_lists, axis=1)
        taxa_df = taxa_df.replace("Unknown Family", "nothing")
        return taxa_df.applymap(self.delete_lowercase)

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
        values_with_taxa_status = df.apply(lambda x: x + " (" + x.name + ")")
        filling_missing_values = values_with_taxa_status.fillna(method='ffill', axis=1)
        combining_df_and_taxa_status_in_nan = df.combine_first(filling_missing_values)
        return combining_df_and_taxa_status_in_nan

    @property
    def standard_taxa_df(self):
        """
        This property can be used to obtain the formated data frame directly, without any intermidiate steps.
        """
        if getattr(self, "_standard_taxa_df", None) is None:
            df = self.import_df(self.filepath)
            taxa_columns = self.spliting_into_taxa_columns(df)
            taxa_column_with_names = self.naming_taxa_columns(taxa_columns)
            complete_taxa_df = self.filling_missing_taxa_values(taxa_column_with_names)
            df_samples = df.drop('#OTU ID', axis=1)
            taxa_columns_and_df_samples = pd.concat([complete_taxa_df, df_samples], axis=1, sort=False)
            standard_taxa_df = taxa_columns_and_df_samples.set_index(list(taxa_column_with_names))
            setattr(self, "_standard_taxa_df", standard_taxa_df.astype(int))
        return self._standard_taxa_df
