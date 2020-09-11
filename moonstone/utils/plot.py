from typing import Union


def check_list_type(items, expected_type: Union[tuple, type]) -> bool:
    """
    check if the type of the different items correspond to the expectedtype(s)
    """
    return all(isinstance(s, expected_type) for s in items)


def add_x_to_plotting_options(plotting_options: dict, option_cat: str, x: str, defaultvalue):
    """
    don't override given plotting_options, meaning it only add the default value
    if value not already defined in plotting_options
    """
    if option_cat not in plotting_options.keys():
        plotting_options[option_cat] = {x: defaultvalue}
        return plotting_options
    elif x not in plotting_options[option_cat].keys():
        plotting_options[option_cat][x] = defaultvalue
        return plotting_options
    else:                               # if x already specified => nothing to change there
        return plotting_options


def add_default_titles_to_plotting_options(plotting_options: dict, title: str, xlabel: str, ylabel: str):
    plotting_options = add_x_to_plotting_options(plotting_options, 'layout', 'title_text', title)
    plotting_options = add_x_to_plotting_options(plotting_options, 'layout', 'title_x', 0.5)
    plotting_options = add_x_to_plotting_options(plotting_options, 'xaxes', 'title_text', xlabel)
    plotting_options = add_x_to_plotting_options(plotting_options, 'yaxes', 'title_text', ylabel)
    return plotting_options
