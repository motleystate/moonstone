from typing import Union, Tuple


"""Operations on dictionnarie(s)."""


def merge_dict(first_dict: dict, second_dict: dict):
    """
    Merging two dictionaries. If key in both, first_dict's key is kept.
    Also deals with dictionaries inside the dictionaries.
    """
    merged_dict = second_dict
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
