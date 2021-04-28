from moonstone.filtering import TaxonomyMeanFiltering
from moonstone.plot.graphs import BarGraph
from moonstone.parsers.base import BaseParser
from moonstone.utils.taxonomy import TaxonomyCountsBase


class BaseTaxonomyCountsParser(TaxonomyCountsBase, BaseParser):
    def plot_most_abundant_taxa(
        self,
        mean_taxa: float = None,
        taxa_number: int = 20,
        taxa_level: str = "species",
        **kwargs,
    ):
        data_df = self.dataframe
        if mean_taxa is not None:
            data_df = TaxonomyMeanFiltering(data_df, mean_taxa).filtered_df
        # Select top `taxa_number`
        top_taxa = (
            data_df.groupby(taxa_level)
            .sum()
            .sum(axis=1)
            .sort_values(ascending=False)[:taxa_number]
            .index
        )
        number_of_taxa = len(
            top_taxa
        )  # Can be different than threshold if not enough taxa
        # Filter for top species
        abundances = data_df.groupby(taxa_level).sum().loc[top_taxa]
        percentage_presence = (abundances != 0).sum(axis=1) / abundances.shape[1] * 100
        # Make graph
        graph = BarGraph(percentage_presence.iloc[::-1])
        # Plotting options
        plotting_options = {
            "layout": {
                "title": f"{number_of_taxa} most abundant {taxa_level} - Total sum of abundances",
                "xaxis_title": "Percentage Sample",
                "yaxis_title": "Species",
            }
        }
        if mean_taxa is not None:
            plotting_options["layout"]["title"] = "{} (mean among samples > {})".format(
                plotting_options["layout"]["title"], mean_taxa
            )
        graph.plot_one_graph(
            orientation="h", plotting_options=plotting_options, **kwargs
        )
