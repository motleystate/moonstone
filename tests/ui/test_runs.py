import os
import shutil
import tempfile
from unittest import TestCase

from moonstone import examples
from moonstone.ui.catalog import TOOLS
from moonstone.ui.forms import to_argv
from moonstone.ui.runs import list_runs, new_run_dir, result_files, run

FEWER_PERMUTATIONS = {'permanova': {'permutations': 99}, 'correlation': {'permutations': 100}}


class TestRuns(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def test_each_example_runs(self):
        for tool in TOOLS.values():
            if not tool.example:
                continue
            with self.subTest(tool=tool.key):
                run_dir = new_run_dir(self.tmp, tool.key)
                meta = examples.path('metadata.csv') if tool.needs_metadata else None
                argv = to_argv(tool, examples.path(tool.example['countfile']), meta, str(run_dir),
                               {**tool.example['values'], **FEWER_PERMUTATIONS.get(tool.key, {})})
                lines = list(run(tool, argv, run_dir))
                self.assertEqual((run_dir / 'exit_code.txt').read_text().strip(), '0', ''.join(lines))
                for pattern in tool.outputs:
                    self.assertTrue(list(run_dir.glob(pattern)), pattern)
                    self.assertTrue(any(path.match(pattern) for path in result_files(run_dir, tool.outputs)), pattern)
                self.assertEqual(list_runs(self.tmp)[0].exit_code, 0)

    def test_two_runs_in_the_same_second_get_two_directories(self):
        first, second = new_run_dir(self.tmp, 'permanova'), new_run_dir(self.tmp, 'permanova')
        self.assertNotEqual(first, second)
        self.assertTrue(first.is_dir() and second.is_dir())

    def test_failed_run_is_kept_and_listed(self):
        run_dir = new_run_dir(self.tmp, 'permanova')
        argv = to_argv(TOOLS['permanova'], '/nope.tsv', '/nope.csv', str(run_dir), {})
        output = ''.join(run(TOOLS['permanova'], argv, run_dir))
        latest = list_runs(self.tmp)[0]
        self.assertEqual((latest.tool_key, latest.exit_code), ('permanova', 1))
        self.assertIn('--permanova', latest.command)
        self.assertIn('nope.tsv', output)
        self.assertIn('nope.tsv', (run_dir / 'run.log').read_text())
        self.assertNotIn('command.txt', [path.name for path in result_files(run_dir)])

    def test_path_with_a_space_runs(self):
        folder = os.path.join(self.tmp, 'Mes données')
        os.mkdir(folder)
        for name in ['species_counts.tsv', 'metadata.csv']:
            shutil.copy(examples.path(name), folder)
        run_dir = new_run_dir(os.path.join(folder, 'runs'), 'correlation')
        argv = to_argv(TOOLS['correlation'], os.path.join(folder, 'species_counts.tsv'),
                       os.path.join(folder, 'metadata.csv'), str(run_dir),
                       {'select': ['Type=inclusion'], 'permutations': 100})
        lines = list(run(TOOLS['correlation'], argv, run_dir))
        self.assertEqual((run_dir / 'exit_code.txt').read_text().strip(), '0', ''.join(lines))

    def test_unfinished_run(self):
        run_dir = new_run_dir(self.tmp, 'associate')
        (run_dir / 'command.txt').write_text('moonstone ...\n')
        self.assertIsNone(list_runs(self.tmp)[0].exit_code)

    def test_no_runs(self):
        self.assertEqual(list_runs(os.path.join(self.tmp, 'nowhere')), [])
