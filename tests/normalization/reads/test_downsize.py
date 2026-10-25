import gzip
import os
import tempfile
from unittest import TestCase

import pandas as pd

from moonstone.normalization.reads.downsize_dir import DownsizeDir, pair_up, plot_reads
from moonstone.normalization.reads.plot_placeholder import get_reads, plotting_info
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
        self.assertEqual(pair_up(info), ['x_1.fq', 'x_2.fq'])


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
