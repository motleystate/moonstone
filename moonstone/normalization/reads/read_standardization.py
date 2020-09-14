import logging
import random

from moonstone.normalization.reads.base import BaseStandardization

logger = logging.getLogger(__name__)


class DownSizePair(BaseStandardization):

    def __init__(self, raw_f, raw_r, n=1000):
        self.raw_f = raw_f
        self.raw_r = raw_r
        self.downsize_to = n

    def random_paired_reads(self):
        if self.raw_f == self.raw_r:
            print(f"\nWarning files {self.raw_f} and {self.raw_r} are the same! Expected Forward and Reverse!\n")

        records: int = sum(1 for _ in open(self.raw_f)) // 4
        rand_reads: list = sorted([random.randint(0, records - 1) for _ in range(self.downsize_to)])

        self.raw_f, self.raw_r = open(self.raw_f), open(self.raw_r)
        suba, subb = open(self.raw_f + ".subset", "w"), open(self.raw_r + ".subset", "w")
        rec_no = - 1
        for rr in rand_reads:

            # Read records in the file (groups of 4 for fastq)
            # Until the record being read is no longer less that one of the ordered random (oxymoron?) records
            while rec_no < rr:
                rec_no += 1
                for i in range(4): self.raw_f.readline()
                for i in range(4): self.raw_r.readline()
            # Once rec_no == rr (we find a matching record), we write that record, forward and reverse.
            for i in range(4):
                suba.write(self.raw_f.readline())
                subb.write(self.raw_r.readline())
            rec_no += 1

        print("wrote to %s, %s" % (suba.name, subb.name))

