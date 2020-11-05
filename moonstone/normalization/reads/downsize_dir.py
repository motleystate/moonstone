import logging
import os
import gzip
import filetype
from moonstone.normalization.reads.read_downsize import DownsizePair

logger = logging.getLogger(__name__)


def pair_up(seq_files_info):
    paired_list = []
    query = None

    for key in seq_files_info:
        if key not in paired_list:  # Any one pairs should generate ONLY one loop.
            query = seq_files_info[key][0]  # set the query to the header of the sequence file
            file_counter = 0  # Confirm that we find two files
            fwd_reads_file = None  # reset the file results
            rev_reads_file = None

            for key_a in seq_files_info:  # Loop through again to find the matching header
                if query in seq_files_info[key_a][0]:  # If the headers match it is either forward or reverse reads.
                    if seq_files_info[key_a][1] == '1':
                        fwd_reads_file = key_a
                        file_counter += 1
                    if seq_files_info[key_a][1] == '2':
                        rev_reads_file = key_a
                        file_counter += 1
                    if fwd_reads_file and rev_reads_file:
                        logger.info('\nForward file = %s\nReverse file = %s' % (fwd_reads_file, rev_reads_file))
                        paired_list.append(fwd_reads_file)
                        paired_list.append(rev_reads_file)

    logger.info(f'List of Paired Reads Files:\n{paired_list}')
    return paired_list


class DownsizeDir:
    """Used to downsize all reads in a given directory to the same number of reads. Reads are downsized by
    random selection of raw reads generating a subset from which alpha diversity can be calculated.
    Note that removal of data, while useful for diversity assessment, is no longer considered good practice.
    https://doi.org/10.1371/journal.pcbi.1003531
    """
    def __init__(self, n=1000, seed=62375, in_dir='./', out_dir=''):
        self.in_dir = in_dir
        self.downsize_to = n
        self.seed = seed

        if out_dir:
            self.out_dir = out_dir
        else:
            self.out_dir = in_dir + 'downsized/'

        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')

    def detect_seq_reads(self):
        """The provided directory might contain files that are not sequence reads.
        This module attempts to identify ONLY files with sequence data."""
        logger.info(f'Detecting sequence files in {self.in_dir}')
        seq_files = [f for f in os.listdir(self.in_dir) if os.path.isfile(self.in_dir + f)]
        logger.info(f'List of Sequencing Files Found:\n{seq_files}')
        return seq_files

    def read_info(self, files):
        """Gather information on the number of reads for each of the sequence reads in the given directory.
        Number of reads can be plotted or reported. Files names and headers are used to match pairs.
        Both compressed and gzipped files are accepted."""

        seq_files_info = {}
        for fh in files:
            detect_type = None
            detect_type = filetype.guess(self.in_dir + fh)
            if not detect_type:
                logger.info('Assuming uncompressed fastq file for %s' % fh)
                file = open(self.in_dir + fh, 'r')
                read_num = sum(1 for _ in open(self.in_dir + fh)) // 4
                header = file.readline().split(' ')[0]
                file.seek(0, 0)
                pair = file.readline().split(' ')[1][0]
                file.close()
                seq_files_info[fh] = [header, pair, read_num, 'Uncompressed/FASTQ']

            if detect_type:
                if detect_type.mime == 'application/gzip':
                    logger.info('Detected gzipped file for %s' % fh)
                    file = gzip.open(self.in_dir + fh, 'r')
                    read_num = sum(1 for _ in gzip.open(self.in_dir + fh)) // 4
                    header = file.readline().decode().split(' ')[0]
                    file.seek(0, 0)
                    pair = file.readline().decode().split(' ')[1][0]
                    file.close()
                    seq_files_info[fh] = [header, pair, read_num, detect_type.mime]

        return seq_files_info

    def down_dir_pair(self):
        files_to_downsize = self.detect_seq_reads()
        file_info_dict = self.read_info(files_to_downsize)
        list_to_downsize = pair_up(file_info_dict)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        else:
            logger.info('Looks like %s exists.' % self.out_dir)
        for k in range(len(list_to_downsize)//2):  # number of files divided by 2: one instance per pair
            instance = DownsizePair(raw_file_f=list_to_downsize[k * 2],
                                    raw_file_r=list_to_downsize[k * 2 + 1],
                                    in_dir=self.in_dir, out_dir=self.out_dir,
                                    n=self.downsize_to)
            if file_info_dict[list_to_downsize[k*2]][3] == 'Uncompressed/FASTQ':
                logger.info(f'{list_to_downsize[k*2]} run uncompressed')
                instance.downsize_pair()

            if file_info_dict[list_to_downsize[k*2]][3] == 'application/gzip':
                logger.info(f'{list_to_downsize[k*2]} run gzip')
                instance.downsize_pair_gzip()

        logger.info('Done!')