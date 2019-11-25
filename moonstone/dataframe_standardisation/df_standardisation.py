import pandas as pd


class Df_Standardisation:

    taxonomical_names = {'0': "kingdom", '1': "phylum", '2': "class",
                         '3': "order", '4': "family", '5': "genus", '6': "species"}

    def __init__(self, df):
        self.df = df

    def import_df(self, filepath):
        df = pd.read_csv(filepath)
        return df.rename(columns=df.loc[0]).drop(df.index[0])

    def delete_lowercase(self, string):
        if string is not None:
            if string.islower():
                return None
            else:
                return string

    def spliting_taxa_columns(self, df):
        taxa_column = df['#OTU ID'].str.split(";", expand=True)
        taxa_column = taxa_column.replace("Ambiguous_taxa", "nothing")
        taxa_in_lists = [taxa_column[i].str.split("_", expand=True) for i in range(taxa_column.shape[1])]
        taxa_df = pd.concat(taxa_in_lists, axis=1)
        taxa_df = taxa_df.replace("Unknown Family", "nothing")
        return taxa_df.applymap(self.delete_lowercase)

    def get_taxa_column_name(self, df):
        if getattr(self, "_get_taxa_column_name", None) is None:
            column_names = [taxonomical_names[x] for x in self.taxonomical_names if x in df.loc[:, 1].values]
            setattr(self, "_tax_get_taxa_column_name", column_names)
        return self._tax_get_taxa_column_name

    def naming_taxa_columns(self, df):
        taxa_df = df[3]
        taxa_df.columns = self.getting_taxa_column_name
        return taxa_df

    def filling_missing_taxa_values(self, df):
        values_with_taxa_status = df.apply(lambda x: x + " (" + x.name + ")")
        filling_missing_values = values_with_taxa_status.fillna(method='ffill', axis=1)
        combining_df_and_taxa_status_in_nan = df.combine_first(filling_missing_values)
        return combining_df_and_taxa_status_in_nan

    @property
    def standard_taxa_df(self, filepath):
        if getattr(self, "_standard_taxa_df", None) is None:
            df = self.import_df(filepath)
            taxa_column = self.spliting_taxa_columns(df)
            taxa_column_with_names = self.naming_taxa_columns(taxa_column)
            taxa_df = self.filling_missing_taxa_values(taxa_column_with_names)
            df_samples = df.drop('#OTU ID', axis=1)
            taxa_columns_and_df_samples = pd.concat([taxa_df, df_samples], axis=1, sort=False)
            standard_taxa_df = taxa_columns_and_df_samples.set_index(self.get_taxa_column_name)
            setattr(self, "_standard_taxa_df", standard_taxa_df)
        return self._standard_taxa_df
