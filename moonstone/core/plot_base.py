class PlotTemplateDataFrame:

    def __init__(self, dataframe, **kwarg):
        self.df = dataframe
        if 'output' in kwargs.keys() and kwargs['output']:
            self.output = kwargs['output']

    def in_bins(df, bins_values, normalize=False):
        if isinstance(df, pd.Series):
            out = pd.cut(df, bins = bins_values)   #put every seq in the appropriate bin
            counts = pd.value_counts(out, normalize=normalize)
            counts = counts.reindex(out.cat.categories)
        else:
            print("ERROR : Wrong given argument type : need a Series")     # what error to give ?
            return "ERROR"
    return counts

    def one_hist(y, title, xlabel, ylabel, **kwargs):
        xvalues = np.arange(1,(len(y))*1.5,1.5)
        xnames=list(y.index)
        
        fig, ax = plt.subplots(figsize=(20,15))
        
        if 'log' in kwargs.keys() and kwargs['log']:
            pyplot.yscale('log')
            pyplot.grid(which='minor', axis='both')
            ax.set_axisbelow(True)  # grid lines are behind the rest
    
        x = plt.bar(xvalues,y,1)
        plt.xticks(xvalues, xnames)
        plt.xticks(rotation=90)
        ax.tick_params(axis='both', which='major', labelsize=20)
        plt.title(title, fontsize=30)
        plt.xlabel(xlabel, fontsize=24, labelpad=20)
        plt.ylabel(ylabel, fontsize=24, labelpad=20)
        plt.show()

        if self.output:
            fig.savefig(self.output)