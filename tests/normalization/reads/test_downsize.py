import gzip
import os
import random
import tempfile
from unittest import TestCase

import pandas as pd

from moonstone.normalization.reads.downsize_dir import DownsizeDir, pair_up, plot_reads
from moonstone.normalization.reads.plot_placeholder import get_reads, make_annotations, plotting_info
from moonstone.normalization.reads.read_downsize import DownsizePair


def fastq_records(n_reads, mate, prefix):
    """Illumina 1.8+ style records: '@<read id> <mate>:N:0:<index>'."""
    return ''.join(f'@{prefix}:{i} {mate}:N:0:ACGT\nACGTACGTAC\n+\nFFFFFFFFFF\n' for i in range(n_reads))


def write_fastq(path, n_reads, mate, prefix, gz=False):
    text = fastq_records(n_reads, mate, prefix)
    if gz:
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    else:
        with open(path, 'w') as handle:
            handle.write(text)


def read_ids(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as handle:
        lines = handle.read().splitlines()
    return [line.split(' ')[0] for line in lines[0::4]]


class TestDownsizePair(TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp() + '/'
        write_fastq(self.dir + 'S1_R1.fastq', 1000, 1, 'M:1:FC:1:1101:100')
        write_fastq(self.dir + 'S1_R2.fastq', 1000, 2, 'M:1:FC:1:1101:100')
        self.read_info = ['@M:1:FC:1:1101:100:0', '1', 1000, 'Uncompressed/FASTQ']

    def downsize(self, out_prefix, n=50, seed=62375):
        pair = DownsizePair('S1_R1.fastq', 'S1_R2.fastq', read_info=self.read_info, in_dir=self.dir,
                            out_dir=self.dir + out_prefix, n=n, seed=seed)
        pair.downsize_pair()
        return (read_ids(self.dir + out_prefix + 'downsized.S1_R1.fastq'),
                read_ids(self.dir + out_prefix + 'downsized.S1_R2.fastq'))

    def test_writes_requested_number_of_reads(self):
        forward, reverse = self.downsize('a.')
        self.assertEqual(len(forward), 50)
        self.assertEqual(len(reverse), 50)

    def test_forward_and_reverse_keep_the_same_reads_in_order(self):
        forward, reverse = self.downsize('a.')
        self.assertEqual(forward, reverse)

    def test_reads_are_whole_records_from_the_input(self):
        forward, _ = self.downsize('a.')
        self.assertTrue(set(forward) <= set(read_ids(self.dir + 'S1_R1.fastq')))
        with open(self.dir + 'a.downsized.S1_R1.fastq') as handle:
            lines = handle.read().splitlines()
        self.assertTrue(all(line == '+' for line in lines[2::4]))
        self.assertTrue(all(line == 'ACGTACGTAC' for line in lines[1::4]))

    def test_same_seed_same_reads(self):
        self.assertEqual(self.downsize('a.', seed=1), self.downsize('b.', seed=1))

    def test_different_seed_different_reads(self):
        self.assertNotEqual(self.downsize('a.', seed=1), self.downsize('b.', seed=2))

    def test_selects_the_drawn_reads(self):
        forward, _ = self.downsize('a.', n=50, seed=7)
        expected = sorted(random.Random(7).sample(range(1000), 50))
        self.assertEqual(forward, [f'@M:1:FC:1:1101:100:{i}' for i in expected])

    def test_downsize_to_the_full_depth_keeps_every_read_once(self):
        # Rarefying to the depth of the smallest sample: that sample keeps all its reads, the first one included.
        forward, reverse = self.downsize('a.', n=1000)
        self.assertEqual(forward, read_ids(self.dir + 'S1_R1.fastq'))
        self.assertEqual(reverse, read_ids(self.dir + 'S1_R2.fastq'))

    def test_more_reads_than_available_raises(self):
        with self.assertRaisesRegex(ValueError, 'to 1001 reads: it has only 1000'):
            self.downsize('a.', n=1001)

    def test_without_read_info_and_input_in_another_directory(self):
        out_dir = tempfile.mkdtemp() + '/'
        DownsizePair('S1_R1.fastq', 'S1_R2.fastq', in_dir=self.dir, out_dir=out_dir, n=50).downsize_pair()
        forward = read_ids(out_dir + 'downsized.S1_R1.fastq')
        self.assertEqual(len(forward), 50)
        self.assertEqual(forward, read_ids(out_dir + 'downsized.S1_R2.fastq'))

    def test_gzip_pair_stays_gzipped(self):
        write_fastq(self.dir + 'S2_R1.fastq.gz', 1000, 1, 'M:1:FC:1:1101:200', gz=True)
        write_fastq(self.dir + 'S2_R2.fastq.gz', 1000, 2, 'M:1:FC:1:1101:200', gz=True)
        DownsizePair('S2_R1.fastq.gz', 'S2_R2.fastq.gz', in_dir=self.dir, out_dir=self.dir, n=50,
                     read_info=['@M:1:FC:1:1101:200:0', '1', 1000, 'application/gzip']).downsize_pair()
        forward = read_ids(self.dir + 'downsized.S2_R1.fastq.gz')
        self.assertEqual(len(forward), 50)
        self.assertEqual(forward, read_ids(self.dir + 'downsized.S2_R2.fastq.gz'))


class TestPairUp(TestCase):

    def test_pairs_by_header_and_mate(self):
        info = {
            'x_2.fq': ['@A:1', '2', 10, 'Uncompressed/FASTQ'],
            'y_R1.fastq': ['@B:7', '1', 10, 'Uncompressed/FASTQ'],
            'x_1.fq': ['@A:1', '1', 10, 'Uncompressed/FASTQ'],
            'y_R2.fastq': ['@B:7', '2', 10, 'Uncompressed/FASTQ'],
        }
        pairs = pair_up(info)
        self.assertEqual(sorted(zip(pairs[0::2], pairs[1::2])), [('x_1.fq', 'x_2.fq'), ('y_R1.fastq', 'y_R2.fastq')])

    def test_unpaired_file_is_left_out(self):
        info = {
            'x_1.fq': ['@A:1', '1', 10, 'Uncompressed/FASTQ'],
            'x_2.fq': ['@A:1', '2', 10, 'Uncompressed/FASTQ'],
            'lonely_R1.fq': ['@C:3', '1', 10, 'Uncompressed/FASTQ'],
        }
        with self.assertLogs('moonstone.normalization.reads.downsize_dir', level='WARNING') as logs:
            self.assertEqual(pair_up(info), ['x_1.fq', 'x_2.fq'])
        self.assertIn("No pair found for 1 file(s), not downsized: ['lonely_R1.fq']", logs.output[0])

    def test_header_that_starts_another_header_is_not_a_match(self):
        info = {
            'E_R1.fq': ['@M:1:FC:1:1101:10737:1266', '1', 10, 'Uncompressed/FASTQ'],
            'E_R2.fq': ['@M:1:FC:1:1101:10737:1266', '2', 10, 'Uncompressed/FASTQ'],
            'F_R1.fq': ['@M:1:FC:1:1101:10737:12660', '1', 10, 'Uncompressed/FASTQ'],
            'F_R2.fq': ['@M:1:FC:1:1101:10737:12660', '2', 10, 'Uncompressed/FASTQ'],
        }
        self.assertEqual(pair_up(info), ['E_R1.fq', 'E_R2.fq', 'F_R1.fq', 'F_R2.fq'])


class TestDownsizeDir(TestCase):

    def setUp(self):
        self.dir = tempfile.mkdtemp() + '/'
        # Pairs named three different ways, found by their read headers, not their names.
        write_fastq(self.dir + 'A_R1.fastq', 200, 1, 'M:1:FC:1:1101:101')
        write_fastq(self.dir + 'A_R2.fastq', 200, 2, 'M:1:FC:1:1101:101')
        write_fastq(self.dir + 'sampleB_1.fq', 300, 1, 'M:1:FC:1:1101:202')
        write_fastq(self.dir + 'sampleB_2.fq', 300, 2, 'M:1:FC:1:1101:202')
        write_fastq(self.dir + 'C.forward.fastq.gz', 250, 1, 'M:1:FC:1:1101:303', gz=True)
        write_fastq(self.dir + 'C.reverse.fastq.gz', 250, 2, 'M:1:FC:1:1101:303', gz=True)

    def test_read_info(self):
        downsizer = DownsizeDir(n=20, in_dir=self.dir)
        self.assertEqual(downsizer.read_info('sampleB_2.fq'),
                         {'sampleB_2.fq': ['@M:1:FC:1:1101:202:0', '2', 300, 'Uncompressed/FASTQ']})
        self.assertEqual(downsizer.read_info('C.forward.fastq.gz'),
                         {'C.forward.fastq.gz': ['@M:1:FC:1:1101:303:0', '1', 250, 'application/gzip']})

    def test_default_output_directory(self):
        downsizer = DownsizeDir(n=20, in_dir=self.dir)
        self.assertEqual(downsizer.out_dir, self.dir + 'downsized/')
        self.assertTrue(os.path.isdir(downsizer.out_dir))

    def test_downsizes_every_pair_to_the_same_depth(self):
        downsizer = DownsizeDir(n=20, in_dir=self.dir)
        downsizer.down_dir_pair()
        outputs = sorted(os.listdir(downsizer.out_dir))
        self.assertEqual(outputs, ['downsized.A_R1.fastq', 'downsized.A_R2.fastq', 'downsized.C.forward.fastq.gz',
                                   'downsized.C.reverse.fastq.gz', 'downsized.sampleB_1.fq', 'downsized.sampleB_2.fq'])
        for forward, reverse in zip(outputs[0::2], outputs[1::2]):
            with self.subTest(pair=forward):
                forward_ids = read_ids(downsizer.out_dir + forward)
                self.assertEqual(len(forward_ids), 20)
                self.assertEqual(forward_ids, read_ids(downsizer.out_dir + reverse))

    def test_files_that_are_not_fastq_are_skipped(self):
        with open(self.dir + 'README.txt', 'w') as handle:
            handle.write('notes about the run\n')
        with open(self.dir + 'logo.png', 'wb') as handle:
            handle.write(b'\x89PNG\r\n\x1a\n' + bytes(32))
        downsizer = DownsizeDir(n=20, in_dir=self.dir)
        with self.assertLogs('moonstone.normalization.reads.downsize_dir', level='WARNING') as logs:
            self.assertEqual(downsizer.read_info('README.txt'), {})
            self.assertEqual(downsizer.read_info('logo.png'), {})
        self.assertIn('Skipping README.txt', logs.output[0])
        self.assertIn('Skipping logo.png: image/png', logs.output[1])
        downsizer.down_dir_pair()
        self.assertEqual(len(os.listdir(downsizer.out_dir)), 6)

    def test_pair_with_too_few_reads_left_out(self):
        downsizer = DownsizeDir(n=220, in_dir=self.dir)
        with self.assertLogs('moonstone.normalization.reads.downsize_dir', level='WARNING') as logs:
            downsizer.down_dir_pair()
        self.assertIn('Leaving out A_R1.fastq and its pair: 200 reads, fewer than 220', '\n'.join(logs.output))
        self.assertEqual(sorted(os.listdir(downsizer.out_dir)),
                         ['downsized.C.forward.fastq.gz', 'downsized.C.reverse.fastq.gz',
                          'downsized.sampleB_1.fq', 'downsized.sampleB_2.fq'])

    def test_plot_reads_table(self):
        info = {'a': ['@A:1', '1', 100, 'x'], 'b': ['@B:1', '1', 300, 'x']}
        expected = pd.DataFrame({'reads': [100, 300]}, index=['a', 'b'])
        pd.testing.assert_frame_equal(plot_reads(info), expected)


class TestPlotPlaceholder(TestCase):

    def test_plotting_info_and_get_reads(self):
        info = {'a': ['@A:1', '1', 100, 'x'], 'b': ['@B:1', '1', 300, 'x']}
        df_info = plotting_info(info)
        pd.testing.assert_frame_equal(df_info, pd.DataFrame({'reads': [100, 300]}, index=['a', 'b']))
        files, reads = get_reads(df_info)
        self.assertEqual(files, ['a', 'b'])
        self.assertEqual(list(reads), [100, 300])

    def test_make_annotations(self):
        df_info = pd.DataFrame({'reads': [100, 300, 200]}, index=['a', 'b', 'c'])
        annotation = make_annotations(df_info)
        for text in ['Number of Samples = 3<br>', 'Total Reads = 600<br>', 'Mean Reads = 200<br>',
                     'Min Reads = 100<br>', 'Max Reads = 300<br>', '50% of samples have at least 200 reads<br>',
                     '90% of samples have at least 120 reads<br>', '25% of samples have at least 250 reads<br>']:
            with self.subTest(text=text):
                self.assertIn(text, annotation)
