import logging
import random

from moonstone.normalization.reads.base import BaseDownsizing

logger = logging.getLogger(__name__)


class DownsizePair(BaseDownsizing):
    """Normalization for the purposes of assessing diversity. Reads are downsized by random selection of raw reads
    generating a subset from which alpha diversity can be calculated.
    Note that removal of data, while useful for diversity assessment, is no longer considered good practice.
    https://doi.org/10.1371/journal.pcbi.1003531
    """

    def __init__(self, raw_f, raw_r, n=1000):
        """Paired reads assumes forward and reverse FASTQ files.
        n is the number of reads that will be randomly picked, with a default of 1000."""

        super().__init__(raw_f, raw_r)
        self.raw_f = raw_f
        self.raw_r = raw_r
        self.downsize_to = n

    def downsize_pair(self):
        if self.raw_f == self.raw_r:
            print(f"\nWarning files {self.raw_f} and {self.raw_r} are the same! Expected Forward and Reverse!\n")

        records: int = sum(1 for _ in open(self.raw_f)) // 4
        rand_reads: list = sorted([random.randint(0, records - 1) for _ in range(self.downsize_to)])

        # forward_reads, reverse_reads = open(self.raw_f), open(self.raw_r)
        downsized_forward, downsized_reverse = \
            open(self.raw_f + ".downsized", "w"), open(self.raw_r + ".downsized", "w")
        rec_no = - 1
        for rr in rand_reads:
            # Read records in the file (groups of 4 for fastq)
            # Until the record being read is no longer less that one of the ordered random (oxymoron?) records
            while rec_no < rr:
                rec_no += 1
                for i in range(4):
                    self.raw_f.readline()
                for i in range(4):
                    self.raw_r.readline()
            # Once rec_no == rr (we find a matching record), we write that record, forward and reverse.
            for i in range(4):
                downsized_forward.write(self.raw_f.readline())
                downsized_reverse.write(self.raw_r.readline())
            rec_no += 1

        print("wrote to %s, %s" % (downsized_forward.name, downsized_reverse.name))
