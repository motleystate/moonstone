import logging
import os
import gzip
import filetype
from moonstone.normalization.reads.read_downsize import DownsizePair

logger = logging.getLogger(__name__)


def pair_up(file_info):
    paired_list = []
    query = None

    for i in file_info:
        if i[0] not in paired_list:  # Pairs should not generate two independent loops.
            query = i[1][0]  # set the query to the header of the sequence file
            file_counter = 0  # Confirm that we find two files
            fwd_reads = None  # reset the file results
            rev_reads = None

            for q in file_info:  # Loop through again to find the matching header
                if query in q[1][0]:  # If the headers match it is either forward or reverse reads.
                    if q[1][1] == '1':
                        fwd_reads = q[0]
                        file_counter += 1
                    if q[1][1] == '2':
                        rev_reads = q[0]
                        file_counter += 1
                    if fwd_reads and rev_reads:
                        logger.info('Forward file = %s\nReverse file = %s' % (fwd_reads, rev_reads))
                        #  print(file_counter)
                        paired_list.append(fwd_reads)
                        paired_list.append(rev_reads)

    return paired_list


class DownsizeDir:
    """Used to downsize all reads in a given directory to the same number of reads. Reads are downsized by
    random selection of raw reads generating a subset from which alpha diversity can be calculated.
    Note that removal of data, while useful for diversity assessment, is no longer considered good practice.
    https://doi.org/10.1371/journal.pcbi.1003531
    """
    def __init__(self, path, n=1000, out_dir='downsized/'):
        self.path = path
        self.downsize_to = n
        self.out_dir = out_dir
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')

    def detect_seq_reads(self):
        """The provided directory might contain files that are not sequence reads.
        This module attempts to identify ONLY files with sequence data."""
        logger.info(f'Detecting sequence files in {self.path}')
        files = [f for f in os.listdir(self.path) if os.path.isfile(self.path + f)]
        return files

    def read_info(self, files):
        """Gather information on the number of reads for each of the sequence reads in the given directory.
        Number of reads can be plotted or reported. Files names and headers are used to match pairs.
        Both compressed and gzipped files are accepted."""

        file_info = []
        for fh in files:
            detect_type = None
            detect_type = filetype.guess(self.path + fh)
            if not detect_type:
                logger.info('Assuming uncompressed fastq file for %s' % fh)
                file = open(self.path + fh, 'r')
                read_num = sum(1 for _ in open(self.path + fh)) // 4
                header = file.readline().split(' ')[0]
                file.seek(0, 0)
                pair = file.readline().split(' ')[1][0]
                file.close()
                file_info.append([self.path + fh, [header, pair, read_num]])

            if detect_type:
                if detect_type.mime == 'application/gzip':
                    logger.info('Detected gzipped file for %s' % fh)
                    file = gzip.open(self.path + fh, 'r')
                    read_num = sum(1 for _ in gzip.open(self.path + fh)) // 4
                    header = file.readline().decode().split(' ')[0]
                    file.seek(0, 0)
                    pair = file.readline().decode().split(' ')[1][0]
                    file.close()
                    file_info.append([self.path + fh, [header, pair, read_num]])

        return file_info

    def down_dir_pair(self):
        files = self.detect_seq_reads()
        information = self.read_info(files)
        to_downsize = pair_up(information)
        for k in range(len(to_downsize)//2):  # number of files divided by 2: one instance per pair
            instance = DownsizePair(raw_file_f=to_downsize[k*2], raw_file_r=to_downsize[k*2+1])
            instance.downsize_pair()
