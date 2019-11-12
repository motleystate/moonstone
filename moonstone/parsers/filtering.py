class Filtering(object):
    def __init__(self, df, mean_filter):
        self.df = df
        self.mean_filter = mean_filter

    def by_mean(self):
        print("Filtering option selected...filtering data by mean reads.")
        df_filtered = self.df[self.df.columns[self.df.mean() >= self.mean_filter]]
        num_vars = df_filtered.shape[1]
        print(f'Retained {num_vars} taxa with at least {self.mean_filter} mean reads.')
        return df_filtered
