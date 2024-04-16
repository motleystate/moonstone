from typing import Union, Tuple


"""Operations on dictionnarie(s)."""
from collections.abc import MutableMapping


def merge_dict(first_dict: dict, second_dict: dict):
    """
    Merging two dictionaries. If key in both, first_dict's key is kept.
    Also deals with dictionaries inside the dictionaries.
    """
    merged_dict = second_dict.copy()
    for key in first_dict.keys():
        if merged_dict.get(key, None) is None:
            merged_dict[key] = first_dict[key]
        else:
            if type(merged_dict[key]) is dict:
                merged_dict[key] = merge_dict(first_dict[key], second_dict[key])
            else:
                merged_dict[key] = first_dict[key]
    return merged_dict


def filter_dict(
    dic: dict, authorized_keys: list, raisingwarning_keys: list = None
) -> Union[dict, Tuple[dict, list]]:
    """
    Args:
        dic: Dictionary to filter.
        authorized_keys: List of keys to keep in dictionary.
        raisingwarning_keys: List of keys that if they aren't on the authorized_keys list but are in dic should
          raise a warning.
    """
    if raisingwarning_keys is None:
        return dict((k, dic[k]) for k in authorized_keys if k in dic)
    else:
        newdic = {}
        raisewarning = []
        for k in dic.keys():
            if k in authorized_keys:
                newdic[k] = dic[k]
            elif k in raisingwarning_keys:
                raisewarning += [k]
        return newdic, raisewarning


def _flatten_dict_gen(d, parent_key, sep, flatten_list):
    """
    src: https://www.freecodecamp.org/news/how-to-flatten-a-dictionary-in-python-in-4-different-ways/
    """
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            yield from flatten_dict(v, new_key, sep=sep).items()
        if flatten_list and isinstance(v, list):
            yield from flatten_dict({str(x): v[x] for x in range(len(v))}, new_key, sep=sep).items()
        else:
            yield new_key, v


def flatten_dict(d: MutableMapping, parent_key: str = '', sep: str = '.', flatten_list: bool = False):
    """
    Flatten an entangled dictionary.

    Example:
    entangled dictionary = {'a': {'b': 1, 'c': 2}}
    flattened dictionary = {
        'a': {'b': 1, 'c': 2},
        'a.b': 1,
        'a.c': 2
    }

    src: https://www.freecodecamp.org/news/how-to-flatten-a-dictionary-in-python-in-4-different-ways/

    Args:
        - d: dictionary to flatten
        - parent_key: name of the keys under which d was registered. It is used as prefix for the keys in d.
        - sep: character or string of character used to separate parent_key/prefix and keys of d.
        - flatten_list: should list be flatten. Key would be parent_key followed by position in list.
          Example: {'a': ['b', 'c']} -> {'a.0': 'b', 'a.1': 'c'} 
    """
    return dict(_flatten_dict_gen(d, parent_key, sep, flatten_list))


def search_for_key_in_entangled_dictionaries(dic, key_of_interest):
    """
    To search for a particular key inside a entangled dictionary without knowing the exact location

    Args:
        - dic: entangled dictionary to search into.
        - key_of_interest: key to search in to search for.
    """
    d = flatten_dict(dic, flatten_list=False)
    d_values = {}
    for k in d.keys():
        if key_of_interest in k.split("."):
            d_values[k] = d[k]
    return d_values


def super_pop(dic: dict, list_of_keys: list, default):
    """
    To remove item at a specified index inside a dictionary.
    Useful if you don't know if the specified index is really there but you don't want to get an error.

    Args:
        - dic: entangled dictionary to search into.
        - list_of_keys: specified index as list of the successive keys inside the entangled dictionary.
        - default: value returned if item is not at the specified index
    """
    for k in list_of_keys:
        if k in dic.keys():
            dic = dic[k]
        else:
            return default
    return dic
