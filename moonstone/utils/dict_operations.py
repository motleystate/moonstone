"""Operations on dictionnarie(s)."""


def merge_dict(first_dict: dict, second_dict: dict):
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
