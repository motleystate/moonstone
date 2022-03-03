from typing import Union

import plotly.graph_objects as go


def check_list_type(items, expected_type: Union[tuple, type]) -> bool:
    """
    check if the type of the different items correspond to the expectedtype(s)
    """
    return all(isinstance(s, expected_type) for s in items)


def add_x_to_plotting_options(
    plotting_options: dict, option_cat: str, x: str, defaultvalue
) -> dict:
    """
    don't override given plotting_options, meaning it only add the default value
    if value not already defined in plotting_options
    """
    if plotting_options is None:
        plotting_options = {}
    if option_cat not in plotting_options.keys():
        plotting_options[option_cat] = {x: defaultvalue}
        return plotting_options
    elif x not in plotting_options[option_cat].keys():
        plotting_options[option_cat][x] = defaultvalue
        return plotting_options
    else:  # if x already specified => nothing to change there
        return plotting_options


def add_default_titles_to_plotting_options(
    plotting_options: dict, title: str, xlabel: str, ylabel: str
) -> dict:
    plotting_options = add_x_to_plotting_options(
        plotting_options, "layout", "title_text", title
    )
    plotting_options = add_x_to_plotting_options(
        plotting_options, "layout", "title_x", 0.5
    )
    plotting_options = add_x_to_plotting_options(
        plotting_options, "xaxes", "title_text", xlabel
    )
    plotting_options = add_x_to_plotting_options(
        plotting_options, "yaxes", "title_text", ylabel
    )
    return plotting_options


def add_groups_annotations(fig, x_coor: list, groups: list) -> go.Figure:
    """
    Args
        fig: The figure to which you want to add labels annotation above the data.
        x_coor: List of x coordinates triplet [(x start of background, x of text annotation, x end of background)].
        groups: List of groups' name, should be in the same order as x_coor.
    """
    i = 0
    color_bg = ["#FFFFFF", "#a7bcdb"]
    while i < len(groups):
        # adding background color
        fig.add_shape(
            type="rect",
            x0=x_coor[i][0],
            y0=100,
            x1=x_coor[i][2],
            y1=104,
            line=dict(
                width=0,
            ),
            fillcolor=color_bg[i % 2],
        )
        # adding text annotation (group name)
        fig.add_annotation(
            x=x_coor[i][1],
            y=102,
            xref="x",
            yref="y",
            text=groups[i],
            showarrow=False,
            font=dict(family="Arial", size=14),
        )
        if i < (len(groups) - 1):
            # adding line separating groups
            fig.add_shape(
                type="line",
                x0=x_coor[i][2],
                y0=100,
                x1=x_coor[i][2],
                y1=0,
                line=dict(width=1, dash="solid", color="white"),
            )
        i += 1
    return fig
