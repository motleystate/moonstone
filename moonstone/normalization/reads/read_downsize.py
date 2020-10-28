import logging
import random
import gzip

from moonstone.normalization.reads.base import BaseDownsizing

logger = logging.getLogger(__name__)


class DownsizePair(BaseDownsizing):
    """Normalization for the purposes of assessing diversity. Reads are downsized by random selection of raw reads
    generating a subset from which alpha diversity can be calculated.
    Note that removal of data, while useful for diversity assessment, is no longer considered good practice.
    https://doi.org/10.1371/journal.pcbi.1003531
    """

    def __init__(self, raw_file_f, raw_file_r, starting_reads=None, in_dir='./', out_dir='./', n=1000, seed=62375):
        """Paired reads assumes forward and reverse FASTQ files.
        n is the number of reads that will be randomly picked, with a default of 1000.
        A random seed is preset to 62375 to allow for reproducibility"""

        super().__init__(raw_file_f, raw_file_r)
        self.downsize_to = n
        self.seed = seed
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.starting_reads = starting_reads

        if self.raw_file_f == self.raw_file_r:
            logger.error(f"Files {self.raw_file_f} and {self.raw_file_r} are the same! Expected Forward and Reverse!")

    def count_starting_reads(self, file):


    def downsize_pair(self):
        """Selects a pseudo-random list of reads from the sequence file and returns the downsized file in the
        same format. The seed for generating the list of reads to select is set during instantiation.
        """

        if not self.starting_reads:
            records: int = sum(1 for _ in open(self.in_dir + self.raw_file_f)) // 4
        else:
            records = self.starting_reads

        logger.info('Found %i reads' % records)
        random.seed(self.seed)
        rand_reads: list = sorted([random.randint(0, records - 1) for _ in range(self.downsize_to)])

        forward_reads = open(self.in_dir + self.raw_file_f, 'r')
        reverse_reads = open(self.in_dir + self.raw_file_r, 'r')
        downsized_forward = open(self.out_dir+"downsized."+self.raw_file_f, "w+")
        downsized_reverse = open(self.out_dir+"downsized."+self.raw_file_r, "w+")

        rec_no = -1
        for rr in rand_reads:
            # Read records in the file (groups of 4 for fastq)
            # Until the record being read is no longer less that one of the ordered random (oxymoron?) records
            while rec_no < rr:
                rec_no += 1
                for i in range(4):
                    forward_reads.readline()
                    reverse_reads.readline()
            # Once rec_no == rr (we find a matching record), we write that record, forward and reverse.
            for i in range(4):
                downsized_forward.write(forward_reads.readline())
                downsized_reverse.write(reverse_reads.readline())
            rec_no += 1

        # Close raw read files.
        forward_reads.close()
        reverse_reads.close()
        # Reset file objects of downsized reads so we can count them. Count them.
        downsized_forward.seek(0)
        downsized_reverse.seek(0)
        downsized_forward_count = sum(1 for _ in downsized_forward) / 4
        downsized_reverse_count = sum(1 for _ in downsized_reverse) / 4
        # Close the downsized files.
        downsized_forward.close()
        downsized_reverse.close()

        logger.info('Wrote %i reads to to %s.\nWrote %i reads to %s' % (downsized_forward_count, downsized_forward.name,
                                                                        downsized_reverse_count, downsized_reverse.name)
                    )

    def downsize_pair_gzip(self):
        """Same as 'downsize_pair' module, but made for gzip compressed files. This module returns files using the
        same compression"""

        records: int = sum(1 for _ in gzip.open(self.in_dir + self.raw_file_f)) // 4
        logger.info('Found %i reads' % records)
        random.seed(self.seed)
        rand_reads: list = sorted([random.randint(0, records - 1) for _ in range(self.downsize_to)])

        forward_reads = gzip.open(self.in_dir + self.raw_file_f, 'rb')
        reverse_reads = gzip.open(self.in_dir + self.raw_file_r, 'rb')
        downsized_forward = gzip.open(self.out_dir+"downsized."+self.raw_file_f, "wb+")
        downsized_reverse = gzip.open(self.out_dir+"downsized."+self.raw_file_r, "wb+")

        rec_no = -1
        reads_written = 0
        for rr in rand_reads:
            while rec_no < rr:
                rec_no += 1
                for i in range(4):
                    forward_reads.readline()
                    reverse_reads.readline()
            for i in range(4):
                downsized_forward.write(forward_reads.readline())
                downsized_reverse.write(reverse_reads.readline())
                reads_written += 1
            rec_no += 1

        forward_reads.close()
        reverse_reads.close()
        downsized_forward.close()
        downsized_reverse.close()
        with gzip.open(self.out_dir+"downsized."+self.raw_file_f, "rb") as f:
            downsized_forward_count = sum(1 for _ in f) / 4
        with gzip.open(self.out_dir+"downsized."+self.raw_file_r, "rb") as r:
            downsized_reverse_count = sum(1 for _ in r) / 4

        logger.info('\nWrote %i reads to to %s.\nWrote %i reads to %s' % (downsized_forward_count,
                                                                          downsized_forward.name,
                                                                          downsized_reverse_count,
                                                                          downsized_reverse.name)
                    )

    def downsize_single(self):
        pass
