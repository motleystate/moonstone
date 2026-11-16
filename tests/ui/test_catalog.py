import argparse
import inspect
import subprocess
import sys
from unittest import TestCase

from moonstone import examples
from moonstone.main import build_parser
from moonstone.ui.catalog import TOOLS, api_objects


class TestCatalog(TestCase):

    def test_every_form_option_exists_in_the_parser(self):
        for tool in TOOLS.values():
            dests = {action.dest for action in tool.build_parser()._actions}
            self.assertLessEqual(set(tool.options) | set(tool.advanced) | set(tool.required), dests, tool.key)

    def test_every_command_line_option_is_on_a_form(self):
        """A flag added to the command line must be added to the catalog."""
        shown = {dest for tool in TOOLS.values() if tool.module == 'moonstone.main'
                 for dest in tool.options + tool.advanced}
        public = {action.dest for action in build_parser()._actions
                  if action.option_strings and action.help != argparse.SUPPRESS}
        self.assertEqual(public - shown - {'help', 'skip_prompt', 'correlation'}, set())

    def test_every_option_of_the_other_commands_is_on_a_form(self):
        for key in ('reads', 'reads_table', 'normalize', 'read_threshold'):
            tool = TOOLS[key]
            public = {action.dest for action in tool.build_parser()._actions
                      if action.option_strings and action.help != argparse.SUPPRESS} - {'help'}
            fixed = {action.dest for action in tool.build_parser()._actions
                     if set(action.option_strings) & set(tool.fixed)}
            placed = set(tool.options) | set(tool.advanced) | set(tool.required) | fixed
            if tool.outdir_flag:
                placed.add(tool.outdir_flag.lstrip('-').replace('-', '_'))
            self.assertEqual(public - placed, set(), key)

    def test_groups_and_order(self):
        self.assertEqual(list(dict.fromkeys(tool.group for tool in TOOLS.values())), ['From reads', 'Prepare', 'Analyse'])
        self.assertEqual([tool.key for tool in TOOLS.values() if tool.group == 'From reads'], ['reads', 'reads_table'])
        self.assertEqual(TOOLS['reads'].input_kind, 'reads_dir')
        self.assertEqual(TOOLS['reads'].choice_labels['remove_host'], {'human': 'Human (9606)', 'mouse': 'Mouse (10090)'})

    def test_examples_exist(self):
        for tool in TOOLS.values():
            if not tool.example:
                continue
            with self.subTest(tool=tool.key):
                self.assertIn(tool.example['countfile'], examples.NAMES)
                self.assertGreater(len(tool.read_samples(examples.path(tool.example['countfile']))), 0)

    def test_every_documented_object_has_a_docstring(self):
        missing = [name for group in api_objects().values() for name, obj in group if not inspect.getdoc(obj)]
        self.assertEqual(missing, [])

    def test_documented_objects_are_listed_once(self):
        names = [name for group in api_objects().values() for name, obj in group]
        self.assertEqual(len(names), len(set(names)))
        self.assertIn('moonstone.parsers.SunbeamKraken2Parser', names)

    def test_catalog_imports_without_streamlit(self):
        code = "import sys; sys.modules['streamlit'] = None; import moonstone.ui.catalog"
        self.assertEqual(subprocess.run([sys.executable, '-c', code]).returncode, 0)
