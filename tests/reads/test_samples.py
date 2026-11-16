import os
import tempfile
from unittest import TestCase

from moonstone.reads.samples import Sample, discover, extension, link, read_sheet, write_sheet


class TestSamples(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def touch(self, *names):
        for name in names:
            path = os.path.join(self.tmp, name)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, 'w').close()

    def codes(self, findings):
        return [finding.code for finding in findings]

    def test_three_naming_conventions(self):
        self.touch('S1_1.fastq.gz', 'S1_2.fastq.gz', 'S2_R1.fastq.gz', 'S2_R2.fastq.gz',
                   'S3_R1_001.fastq.gz', 'S3_R2_001.fastq.gz')
        samples, findings = discover(self.tmp)
        self.assertEqual(findings, [])
        self.assertEqual([s.name for s in samples], ['S1', 'S2', 'S3'])
        self.assertEqual(os.path.basename(samples[2].read2), 'S3_R2_001.fastq.gz')

    def test_underscore_in_sample_name(self):
        self.touch('A_1_1.fastq.gz', 'A_1_2.fastq.gz')
        samples, findings = discover(self.tmp)
        self.assertEqual([(s.name, os.path.basename(s.read1)) for s in samples], [('A_1', 'A_1_1.fastq.gz')])
        self.assertEqual(findings, [])

    def test_sub_folders_and_match(self):
        self.touch('run1/x_1.fq.gz', 'run1/x_2.fq.gz', 'run2/y_1.fq.gz', 'run2/y_2.fq.gz', 'run1/notes.txt')
        self.assertEqual([s.name for s in discover(self.tmp)[0]], ['x', 'y'])
        self.assertEqual([s.name for s in discover(self.tmp, match='^x')[0]], ['x'])

    def test_single_end(self):
        self.touch('S1.fastq.gz', 'S2.fastq.gz')
        samples, findings = discover(self.tmp, paired=False)
        self.assertEqual([(s.name, s.read2) for s in samples], [('S1', None), ('S2', None)])
        self.assertEqual(findings, [])

    def test_no_read_file(self):
        self.touch('notes.txt')
        samples, findings = discover(self.tmp)
        self.assertEqual((samples, self.codes(findings)), ([], ['no_reads']))
        self.assertEqual(self.codes(discover('/nonexistent')[1]), ['not_a_directory'])

    def test_mate_without_pair_is_named(self):
        self.touch('S1_1.fastq.gz', 'S1_2.fastq.gz', 'S2_1.fastq.gz')
        findings = discover(self.tmp)[1]
        self.assertEqual(self.codes(findings), ['mate_missing'])
        self.assertIn('S2', findings[0].message)

    def test_duplicate_sample_names_both_paths(self):
        self.touch('run1/S1_1.fastq.gz', 'run1/S1_2.fastq.gz', 'run2/S1_1.fastq.gz', 'run2/S1_2.fastq.gz')
        findings = discover(self.tmp)[1]
        self.assertEqual(set(self.codes(findings)), {'duplicate_sample'})
        self.assertIn('run1', findings[0].message)
        self.assertIn('run2', findings[0].message)

    def test_mixed_paired_and_single(self):
        self.touch('S1_1.fastq.gz', 'S1_2.fastq.gz', 'lonely.fastq.gz')
        findings = discover(self.tmp)[1]
        self.assertEqual(self.codes(findings), ['not_paired'])
        self.assertIn('lonely.fastq.gz', findings[0].message)

    def test_mixed_extensions(self):
        self.touch('S1_1.fastq.gz', 'S1_2.fastq.gz', 'S2_1.fq', 'S2_2.fq')
        self.assertEqual(self.codes(discover(self.tmp)[1]), ['mixed_extensions'])

    def test_unsafe_sample_name(self):
        self.touch('my sample_1.fastq.gz', 'my sample_2.fastq.gz')
        self.assertEqual(self.codes(discover(self.tmp)[1]), ['sample_name'])

    def test_link_uses_the_names_the_workflow_expects(self):
        self.touch('run1/S3_R1_001.fq.gz', 'run1/S3_R2_001.fq.gz')
        samples = discover(self.tmp)[0]
        target = os.path.join(self.tmp, 'links')
        link(samples, target)
        self.assertEqual(sorted(os.listdir(target)), ['S3_1.fq.gz', 'S3_2.fq.gz'])
        self.assertEqual(os.path.realpath(os.path.join(target, 'S3_2.fq.gz')), os.path.realpath(samples[0].read2))
        self.assertEqual(extension(samples[0].read1), 'fq.gz')

    def test_link_never_overwrites(self):
        self.touch('S1_1.fastq.gz', 'S1_2.fastq.gz')
        samples = discover(self.tmp)[0]
        target = os.path.join(self.tmp, 'links')
        link(samples, target)
        with self.assertRaises(FileExistsError):
            link(samples, target)

    def test_sheet_round_trip(self):
        samples = [Sample('S1', '/a/S1_1.fastq.gz', '/a/S1_2.fastq.gz'), Sample('S2', '/a/S2.fastq.gz', None)]
        sheet = os.path.join(self.tmp, 'samples.tsv')
        write_sheet(samples, sheet)
        self.assertEqual(open(sheet).readline(), 'sample\tread1\tread2\n')
        self.assertEqual(read_sheet(sheet), samples)
