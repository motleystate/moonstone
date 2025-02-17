from typing import Union

import plotly.graph_objects as go

from moonstone.utils.dict_operations import merge_dict


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
    plotting_options = plotting_options.copy()   # so that it doesn't affect the original plotting_options
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


def add_default_titles_to_plotting_options_3d(
    plotting_options: dict, title: str, xlabel: str, ylabel: str, zlabel: str
) -> dict:
    plotting_options = add_x_to_plotting_options(
        plotting_options, "layout", "title_text", title
    )
    scene_dic = {"xaxis_title": xlabel, "yaxis_title": ylabel, "zaxis_title": zlabel}
    if "scene" in plotting_options["layout"].keys():
        scene_dic = merge_dict(plotting_options["layout"]["scene"], scene_dic)

    plotting_options["layout"]["scene"] = scene_dic

    return plotting_options


def add_groups_annotations(
    fig, x_coor: list, y_coor: tuple, groups: list,
    color_bg=["#FFFFFF", "#a7bcdb"],
) -> go.Figure:
    """
    Add labels annotations, in the form of rectangles for the label with text annotations and a line separating the
    groups going from the bottom of the rectangle to the bottom of the data.

    Args
        fig: The figure to which you want to add labels annotation above the data.
        x_coor: List of x coordinates triplet [(x start of background, x of text annotation, x end of background)].
        y_coor: Tuple of 3 y coordinates : (y bottom of line, y bottom of rectangle, y top of rectangle)
        groups: List of groups' name, should be in the same order as x_coor.
    """
    i = 0
    y_med = (y_coor[1] + y_coor[2])/2
    n_cbg = len(color_bg)
    while i < len(groups):
        # adding background color
        fig.add_shape(
            type="rect",
            x0=x_coor[i][0],
            y0=y_coor[1],
            x1=x_coor[i][2],
            y1=y_coor[2],
            line=dict(
                width=0,
            ),
            fillcolor=color_bg[i % n_cbg],
        )
        # adding text annotation (group name)
        fig.add_annotation(
            x=x_coor[i][1],
            y=y_med,
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
                y0=y_coor[1],
                x1=x_coor[i][2],
                y1=y_coor[0],
                line=dict(width=1, dash="solid", color="white"),
            )
        i += 1
    return fig
