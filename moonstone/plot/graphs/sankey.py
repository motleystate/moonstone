from typing import Union

import plotly.graph_objects as go

from moonstone.plot.graphs.base import BaseGraph


class SankeyDiagram(BaseGraph):
    def plot_one_graph(
        self,
        plotting_options: dict = None,
        show: bool = True,
        output_file: Union[bool, str] = False,
    ) -> go.Figure:
        # @TODO: Figure out how to order the values in the columns
        # + add option to personalize colors

        # in case different columns have the same values
        sankey_df = self.data.apply(lambda col: col.astype(str) + " (col=" + col.name + ")", axis=0)
        columns = sankey_df.columns

        # Getting unique labels
        tmp = sankey_df.stack(dropna=False).reset_index(drop=True)
        vc_labels = tmp.value_counts()
        labels = tmp.unique()

        # Mapping labels to indices
        label_map = {label: i for i, label in enumerate(labels)}

        # Creating source-target pairs
        source = []
        target = []
        values = []
        for i in range(1, len(columns)):
            vc = sankey_df[[columns[i-1], columns[i]]].value_counts().reset_index()
            for c in [columns[i-1], columns[i]]:
                vc[c] = vc[c].astype("category")
                vc[c] = vc[c].cat.set_categories(labels, ordered=True)
            vc = vc.sort_values(by=[columns[i-1], columns[i]])
            source += list(vc[columns[i-1]].map(label_map))
            target += list(vc[columns[i]].map(label_map))
            values += list(vc["count"])

        node_label = [k.split(" (col=")[0]+"<br>"+str(vc_labels.loc[k]) for k in label_map.keys()]
        x = None
        y = None

        fig = go.Figure(go.Sankey(
                arrangement='snap',
                node=dict(
                        pad=15, thickness=20,
                        line=dict(color='black', width=0.5),
                        label=node_label,
                        x = x, y = y,
                ),
                link=dict(
                        source=source, target=target, value=values,
                        color="rgba(203, 204, 204, 0.2)",
                )
        ))

        if plotting_options is not None:
            fig = self._handle_plotting_options_plotly(fig, plotting_options)

        self._handle_output_plotly(fig, show, output_file)

        return fig


