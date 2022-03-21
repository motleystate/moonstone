from unittest import TestCase

from moonstone.utils.dict_operations import merge_dict, filter_dict


class TestMergeDict(TestCase):
    def test_simple_dicts(self):
        first_dict = {"a": 3, "b": 4}
        second_dict = {"a": 5, "c": 7}
        expected_dict = {"a": 3, "b": 4, "c": 7}
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)

    def test_with_subdict(self):
        first_dict = {"a": 3, "sub": {"b": 4, "c": 6}}
        second_dict = {"a": 4, "sub": {"b": 7, "d": 6}}
        expected_dict = {"a": 3, "sub": {"b": 4, "c": 6, "d": 6}}
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)

    def test_with_subdict_variant_one(self):
        first_dict = {"a": 3, "sub": {"b": 4, "c": 6}}
        second_dict = {"a": 4}
        expected_dict = {"a": 3, "sub": {"b": 4, "c": 6}}
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)

    def test_with_subdict_variant_two(self):
        first_dict = {"a": 3}
        second_dict = {"a": 4, "sub": {"b": 4, "c": 6}}
        expected_dict = {"a": 3, "sub": {"b": 4, "c": 6}}
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)

    def test_with_multiple_subdict(self):
        first_dict = {"a": 3, "sub": {"b": 4, "subsub": {"z": 3}}}
        second_dict = {"a": 4, "sub": {"b": 7, "d": 6, "subsub": {"z": 5}}}
        expected_dict = {"a": 3, "sub": {"b": 4, "d": 6, "subsub": {"z": 3}}}
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)


class TestFilterDict(TestCase):

    def test_with_raisingwarning_keys(self):
        old_dic = {"ananas": 3, "banana": 4, "coconut": 7, "durian": 6}
        authorized_keys = ["ananas", "coconut"]
        raisingwarning_keys = ["durian"]
        new_dic, raisewarning = filter_dict(
            old_dic,
            authorized_keys,
            raisingwarning_keys
        )
        self.assertDictEqual(new_dic, {"ananas": 3, "coconut": 7})
        self.assertListEqual(raisewarning, ["durian"])

    def test_without_raisingwarning_keys(self):
        old_dic = {"ananas": 3, "banana": 4, "coconut": 7, "durian": 6}
        authorized_keys = ["ananas", "coconut"]
        new_dic = filter_dict(
            old_dic,
            authorized_keys
        )
        self.assertDictEqual(new_dic, {"ananas": 3, "coconut": 7})
