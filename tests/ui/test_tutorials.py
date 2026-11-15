from unittest import TestCase

from moonstone.ui.catalog import TOOLS
from moonstone.ui.tutorials import TUTORIALS, Step, Tutorial, run_tutorial

# The classes of docs/user_guide/available_classes.rst
LISTED = {
    'MetadataParser', 'YAMLBasedMetadataParser', 'GeneCountsParser', 'Picrust2PathwaysParser', 'SunbeamKraken2Parser',
    'Metaphlan2Parser', 'Metaphlan3Parser', 'MetaphlanReadStatsParser', 'Qiime2Parser', 'NoCountsFiltering',
    'NamesFiltering', 'NaNPercentageFiltering', 'NumberOfDifferentValuesFiltering', 'TaxonomyNamesFiltering',
    'TaxonomyMeanFiltering', 'MeanFiltering', 'GeometricMeanNormalization', 'TotalCountsNormalization',
    'RandomSelection', 'TaxonomyRandomSelection', 'PlotCountsStats', 'PlotTaxonomyCounts', 'PlotMetadataStats',
    'Classification', 'DifferentialAnalysis', 'ShannonIndex', 'SimpsonInverseIndex', 'Chao1Index',
    'FaithsPhylogeneticDiversity', 'BrayCurtis', 'WeightedUniFrac', 'UnweightedUniFrac',
}
# Described in a tutorial but not run: no example file in their format (or no phylogenetic tree) is shipped.
NOT_RUN = {
    'YAMLBasedMetadataParser', 'Picrust2PathwaysParser', 'Metaphlan2Parser', 'Metaphlan3Parser',
    'MetaphlanReadStatsParser', 'Qiime2Parser', 'FaithsPhylogeneticDiversity', 'WeightedUniFrac', 'UnweightedUniFrac',
}


class TestTutorials(TestCase):

    def test_every_tutorial_runs_on_the_examples(self):
        for tutorial in TUTORIALS.values():
            with self.subTest(tutorial=tutorial.key):
                results = run_tutorial(tutorial)
                self.assertEqual(len(results), len(tutorial.steps))
                for step, (value, printed) in zip(tutorial.steps, results):
                    if step.code:
                        self.assertTrue(value is not None or printed, f'{tutorial.key}: a code step shows nothing')

    def test_every_tool_has_a_tutorial(self):
        self.assertEqual({tutorial.tool for tutorial in TUTORIALS.values() if tutorial.tool}, set(TOOLS))

    def test_every_listed_class_is_used_in_a_tutorial(self):
        code = '\n'.join(step.code for tutorial in TUTORIALS.values() for step in tutorial.steps)
        text = '\n'.join(step.text for tutorial in TUTORIALS.values() for step in tutorial.steps)
        for name in LISTED - NOT_RUN:
            self.assertIn(name, code)
        for name in NOT_RUN:
            self.assertIn(name, text)

    def test_steps_share_one_namespace_and_show_the_last_expression(self):
        tutorial = Tutorial('x', 'X', 'g', 's', (Step('a', 'a = 2\nprint("hello")'), Step('b', 'a * 21'), Step('c')))
        self.assertEqual(run_tutorial(tutorial), [(None, 'hello\n'), (42, ''), (None, '')])

    def test_a_failing_step_is_named(self):
        broken = Tutorial('x', 'X', 'g', 's', (Step('ok', '1 + 1'), Step('bad', 'undefined_name')))
        with self.assertRaisesRegex(RuntimeError, "'x', step 2 .*NameError"):
            run_tutorial(broken)
