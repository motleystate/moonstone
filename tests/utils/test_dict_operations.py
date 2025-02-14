from unittest import TestCase

from moonstone.utils.dict_operations import (
    filter_dict,
    merge_dict,
    flatten_dict,
    search_for_key_in_entangled_dictionaries,
    super_pop
)


class TestMergeDict(TestCase):

    def test_simple_dicts(self):
        first_dict = {"a": 3, "b": 4}
        second_dict = {"a": 5, "c": 7}
        expected_second_dict = {"a": 5, "c": 7}
        expected_dict = {"a": 3, "b": 4, "c": 7}

        # Testing the merging
        self.assertDictEqual(merge_dict(first_dict, second_dict), expected_dict)

        # ... but also testing that second_dict has not been changed
        self.assertDictEqual(second_dict, expected_second_dict)

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


class TestFlattenDict(TestCase):

    def test_flatten_dict(self):
        tested_dict = {'a': {'b': 1, 'c': 2}, 'd': 3, 'e': {'f': 4, 'g': {'h': 5}, 'i': 6}}
        expected_dict = {
            'a.b': 1,
            'a.c': 2,
            'a': {'b': 1, 'c': 2},
            'd': 3,
            'e.f': 4,
            'e.g.h': 5,
            'e.g': {'h': 5},
            'e.i': 6,
            'e': {'f': 4, 'g': {'h': 5}, 'i': 6}
        }
        self.assertDictEqual(flatten_dict(tested_dict), expected_dict)
        # Checking tested_dict is unchanged
        self.assertDictEqual(
            tested_dict,
            {'a': {'b': 1, 'c': 2}, 'd': 3, 'e': {'f': 4, 'g': {'h': 5}, 'i': 6}}
            )

    def test_flatten_dict_empty(self):
        tested_dict = {}
        expected_dict = {}
        self.assertDictEqual(flatten_dict(tested_dict), expected_dict)

    def test_flatten_dict_without_flattening_list(self):
        tested_dict = {'a': ['b', 'c'], 'd': 3, 'e': {'f': 4, 'g': {}, 'i': 6}}
        expected_dict = {
            'a': ['b', 'c'],
            'd': 3,
            'e.f': 4,
            'e.g': {},
            'e.i': 6,
            'e': {'f': 4, 'g': {}, 'i': 6}
        }
        self.assertDictEqual(flatten_dict(tested_dict), expected_dict)

    def test_flatten_dict_with_flattening_list(self):
        tested_dict = {'a': ['b', 'c'], 'd': 3, 'e': {'f': 4, 'g': {}, 'i': 6}}
        expected_dict = {
            'a.0': 'b',
            'a.1': 'c',
            'd': 3,
            'e.f': 4,
            'e.g': {},
            'e.i': 6,
            'e': {'f': 4, 'g': {}, 'i': 6}
            }
        self.assertDictEqual(flatten_dict(tested_dict, flatten_list=True), expected_dict)
        # Checking tested_dict is unchanged
        self.assertDictEqual(
            tested_dict,
            {'a': ['b', 'c'], 'd': 3, 'e': {'f': 4, 'g': {}, 'i': 6}}
            )

    def test_search_for_key_in_entangled_dictionaries(self):
        tested_dict = {'a': ['b', 'c'], 'd': 3, 'e': {'f': 4, 'g': {'h': 5}, 'i': {}}, 'h': 6}

        # Simple case: first level key of entangled dictionary that appears only once
        self.assertDictEqual(
            search_for_key_in_entangled_dictionaries(tested_dict, 'd'),
            {'d': 3}
            )

        # Case: key that appears in a later level of the entangled dictionary
        self.assertDictEqual(
            search_for_key_in_entangled_dictionaries(tested_dict, 'f'),
            {'e.f': 4}
            )

        # Case: key that appears more than once in the entangled dictionary
        self.assertDictEqual(
            search_for_key_in_entangled_dictionaries(tested_dict, 'h'),
            {'e.g.h': 5, 'h': 6}
            )

        # Case: key that DOESN'T appears in the entangled dictionary
        self.assertDictEqual(
            search_for_key_in_entangled_dictionaries(tested_dict, 'b'),  # 'b' is not a key
            {}
            )

    def test_super_pop(self):
        tested_dict = {'a': ['b', 'c'], 'd': 3, 'e': {'f': 4, 'g': {'h': 5}, 'i': {}}, 'h': 6}

        # Case: first level key of entangled dictionary
        # Not a real practical case because it would be more efficient to use the python build-in method pop
        self.assertEqual(super_pop(tested_dict, ["d"], None), 3)

        # Case: key that appears in a later level of the entangled dictionary
        self.assertEqual(super_pop(tested_dict, ["e", "f"], None), 4)

        # Case: key that doesn't appear at the last level
        self.assertEqual(super_pop(tested_dict, ["e", "k"], None), None)

        # Case: key that doesn't appear at the first level
        self.assertEqual(super_pop(tested_dict, ["k", "k"], None), None)
