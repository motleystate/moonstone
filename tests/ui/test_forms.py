from unittest import TestCase

from moonstone.ui.catalog import TOOLS
from moonstone.ui.forms import check_command, form_actions, preview, to_argv


class TestForms(TestCase):

    def test_form_actions_follow_the_catalog_order(self):
        shown, advanced = form_actions(TOOLS['permanova'])
        self.assertEqual(tuple(a.dest for a in shown), TOOLS['permanova'].options)
        self.assertEqual(tuple(a.dest for a in advanced), TOOLS['permanova'].advanced)

    def test_example_commands_are_valid(self):
        for tool in TOOLS.values():
            argv = to_argv(tool, 'counts.tsv', 'meta.csv', 'out', tool.example['values'])
            self.assertIsNone(check_command(tool, argv), tool.key)

    def test_permanova_command(self):
        argv = to_argv(TOOLS['permanova'], 'c.tsv', 'm.csv', 'out',
                       {'permanova': ['disease'], 'covariates': ['age'], 'taxa_level': 'species', 'seed': 0})
        self.assertEqual(argv, ['c.tsv', 'm.csv', 'out', '-os', '--permanova', 'disease', '--covariates', 'age'])

    def test_permanova_with_no_variable_still_sends_the_flag(self):
        self.assertEqual(to_argv(TOOLS['permanova'], 'c', 'm', 'o', {})[-1], '--permanova')

    def test_flags_numbers_and_lists(self):
        argv = to_argv(TOOLS['classify'], 'c', 'm', 'o',
                       {'classify': 'disease', 'model': ['svm', 'random_forest'], 'roc': True, 'folds': 3,
                        'explore_kernels': False, 'filtering': None})
        self.assertEqual(argv, ['c', 'm', 'o', '-os', '--classify', 'disease', '--model', 'svm', 'random_forest',
                                '--roc', '--folds', '3'])

    def test_read_threshold_command(self):
        self.assertEqual(to_argv(TOOLS['read_threshold'], 'c.csv', None, 'out', {}),
                         ['c.csv', '--save-plots', '--save-filtered', '--output-dir', 'out'])

    def test_missing_required_value(self):
        argv = to_argv(TOOLS['classify'], 'c', 'm', 'o', {})
        self.assertEqual(check_command(TOOLS['classify'], argv), '--classify is required.')

    def test_parser_message_is_returned(self):
        argv = to_argv(TOOLS['correlation'], 'c', 'm', 'o', {'alpha': 2.0})
        self.assertEqual(check_command(TOOLS['correlation'], argv), '--alpha must be between 0 and 1.')

    def test_preview_quotes_paths(self):
        text = preview(TOOLS['permanova'], ['Mes données/c.tsv', 'm.csv', 'out', '-os', '--permanova'])
        self.assertEqual(text, "moonstone 'Mes données/c.tsv' m.csv out -os --permanova")
        self.assertTrue(preview(TOOLS['read_threshold'], ['c.csv']).startswith(
            'python -m moonstone.filtering.read_threshold_analysis '))
