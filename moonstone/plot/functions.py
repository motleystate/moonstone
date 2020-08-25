import sys
import typing
from typing import List
import warnings


def check_list_of(listx, expectedtype):    # expectedtype can be tuple of type
    """
    check the type of what's in listx
    """
    return all(isinstance(s, expectedtype) for s in listx)


def check_type(to_check, expectedtype: type):
    """
    check if type of 'to_check' is right
    """
    # type(List[str]) = typing._GenericAlias in python3.7 and typing.GenericMeta in python3.6
    if sys.version_info[1] >= 7:
        typeList = typing._GenericAlias
        # typeListName = expectedtype._name       # for now only typing type is List so not necessary
    else:
        typeList = typing.GenericMeta
        # typeListName = str(p._gorg).split('.')[1]
    if type(expectedtype) == typeList:      # if List[<type>] (<type> = str, int etc.)
        if type(to_check) == list:
            if check_list_of(to_check, expectedtype.__args__[0]) is True:
                return True
        return 'List['+expectedtype.__args__[0].__name__+']'
    else:
        if type(to_check) == expectedtype:
            return True
        else:
            return expectedtype.__name__


def check_type_s(to_check, expectedtype_s: type):
    """
    check if type of 'to_check' is right.
    Dispatcher if several types are allowed (then expectedtype_s is a list) or not
    """
    if type(expectedtype_s) == list:    # if more than one type allowed
        type_in_str_for_warning = []
        for i in expectedtype_s:
            s = check_type(to_check, i)
            if s is True:
                return True
            else:
                type_in_str_for_warning += [s]
        return " or ".join(type_in_str_for_warning)   # return what will be in the error raised
    else:                               # if only one type
        return check_type(to_check, expectedtype_s)                        # return what will be in the error raised


def check_types_in_plotting_options(plotting_options: dict):
    """
    check types of plotting options given by the user.
    Create a new dictionary cleaned_plotting_options without the options that are of the wrong type
    """
    expectedtype = {'log': bool, 'colorbar': [str, List[str]], 'tickangle': [int, float]}
    cleaned_plotting_options = {}
    for i in plotting_options.keys():
        if i in expectedtype.keys():
            s = check_type_s(plotting_options[i], expectedtype[i])
            if s is True:
                cleaned_plotting_options[i] = plotting_options[i]
            else:
                if type(plotting_options[i]) == list and 'List' in s:
                    # warning when the problem is not the overall type but the type of elements in list
                    warnings.warn(('Warning : %s value in plotting_options must be a %s. \
Please check the type of elements in the list given. Value overridden') % (i, s))
                else:
                    warnings.warn('Warning : %s value in plotting_options must be a %s, %s given. Value overridden' %
                                  (i, s, type(plotting_options[i]).__name__))
    return cleaned_plotting_options


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
