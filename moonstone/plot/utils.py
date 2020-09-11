from typing import Union


def check_list_type(items, expected_type: Union[tuple, type]) -> bool:
    """
    check if the type of the different items correspond to the expectedtype(s)
    """
    return all(isinstance(s, expected_type) for s in items)


def add_x_to_plotting_options(plotting_options: dict, x: str, defaultvalue):
    """
    don't override given plotting_options, meaning it only add the default value
    if value not already defined in plotting_options
    ~~~ MAYBE DELETE THIS METHOD ~~~
    """
    if x not in plotting_options.keys():    # if x not already specified (not given or not of the right type)
        plotting_options[x] = defaultvalue
        return plotting_options
    else:
        return plotting_options
