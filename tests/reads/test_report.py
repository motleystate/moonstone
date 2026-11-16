import importlib.util
import os
import subprocess
import sys
import tempfile
from unittest import TestCase

from moonstone import __version__
from moonstone.reads import WORKFLOW_DIR

SCRIPT = os.path.join(WORKFLOW_DIR, 'scripts', 'make_report.py')


def write_report(path, unclassified, root, human=0, species=()):
    """A Kraken2-style report: percent, clade reads, direct reads, rank, taxid, name."""
    rows = [(unclassified, 0, 'U', 0, 'unclassified'), (root, 0, 'R', 1, 'root')]
    if human:
        rows.append((human, human, 'S', 9606, 'Homo sapiens'))
    rows += [(n, n, 'S', 1000 + i, 'sp%d' % i) for i, n in enumerate(species)]
    with open(path, 'w') as handle:
        for row in rows:
            handle.write('0.00\t%d\t%d\t%s\t%d\t%s\n' % row)


class TestRunReport(TestCase):
    """The run report, on small hand-written Kraken2 reports (ported from Moonbeam's tests/test_report.py)."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.host, self.final, self.bracken = (os.path.join(self.tmp, name) for name in ('host', 'final', 'brk'))
        for folder in (self.host, self.final, self.bracken):
            os.mkdir(folder)
        # good: 1,000,000 reads, 25% host, 60% of the remaining classified
        write_report(self.host + '/good.taxa.report', 0, 1000000, human=250000)
        write_report(self.final + '/good.taxa.report', 300000, 450000)
        # bad: few reads, almost nothing classified
        write_report(self.host + '/bad.taxa.report', 0, 5000, human=10)
        write_report(self.final + '/bad.taxa.report', 4900, 90)
        for sample, species in (('good', (10, 20, 30)), ('bad', (5,))):
            write_report(self.bracken + '/%s.brackenCorrected_taxa.report' % sample, 0, sum(species), species=species)

    def make(self, *arguments):
        out = os.path.join(self.tmp, 'report.md')
        subprocess.run([sys.executable, SCRIPT, *arguments, '--out', out], check=True)
        with open(out) as handle:
            return handle.read()

    def test_rows_and_warnings(self):
        table = os.path.join(self.tmp, 'all.tsv')
        with open(table, 'w') as handle:
            handle.write('# Constructed from biom file\n#OTU ID\tgood.taxa\n')      # 'bad' is missing
        text = self.make('--samples', 'good', 'bad', '--final-dir', self.final, '--host-dir', self.host,
                         '--bracken-dir', self.bracken, '--tables', table, '--settings', '{"paired": true}')
        self.assertIn('| good | 1000000 | 25.0 | 750000 | 60.0 | 40.0 | 3 |', text)
        self.assertIn('bad: only 4990 reads', text)
        self.assertIn('bad: 1.8% classified', text)
        self.assertIn('bad: not found in all.tsv', text)
        self.assertNotIn('good:', text.split('## Warnings')[1].split('## Output')[0])

    def test_several_host_taxids_are_summed(self):
        write_report(self.host + '/two.taxa.report', 0, 1000, human=100)
        with open(self.host + '/two.taxa.report', 'a') as handle:
            handle.write('0.00\t200\t200\tS\t10090\tMus musculus\n')
        write_report(self.final + '/two.taxa.report', 0, 700)
        text = self.make('--samples', 'two', '--final-dir', self.final, '--host-dir', self.host,
                         '--settings', '{"host_taxid": [9606, 10090]}')
        self.assertIn('| two | 1000 | 30.0 | 700 |', text)

    def test_without_host_mapping_or_bracken(self):
        text = self.make('--samples', 'good', '--final-dir', self.final)
        self.assertIn('| good | NA | NA | 750000 | 60.0 | 40.0 |', text)
        self.assertNotIn('Species (Bracken)', text)

    def test_report_records_the_moonstone_version(self):
        text = self.make('--samples', 'good', '--final-dir', self.final)
        self.assertIn('- Moonstone: ' + __version__, text)
        self.assertNotIn('commit', text)

    def test_a_missing_tool_does_not_break_the_report(self):
        spec = importlib.util.spec_from_file_location('make_report', SCRIPT)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(module.tool_version(['a-tool-that-does-not-exist', '--version']), 'not found')
        module.tool_version = lambda command: 'Bracken v3.1'
        self.assertEqual(module.bracken_version(), 'Bracken v3.1')
        module.tool_version = lambda command: 'bracken: illegal option -- v'      # versions before 3.x have no -v
        self.assertTrue(module.bracken_version().startswith('unknown'))
