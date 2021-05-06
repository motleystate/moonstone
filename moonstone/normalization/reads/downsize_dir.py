import logging
import os
import gzip
import filetype
import multiprocessing as mp
import pandas as pd
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


def plot_reads(file_info_dict):
    logger.info('Generating plot of number of reads')
    # generate a dataframe from the file information dictionary
    # to include filename as the index
    files: list = []
    reads: list = []
    for key in file_info_dict:
        files.append(key)
        reads.append(file_info_dict[key][2])

    df = pd.DataFrame(index=files, data=reads, columns=['reads'])
    return df


class DownsizeDir:
    """Used to downsize all reads in a given directory to the same number of reads. Reads are downsized by
    random selection of raw reads generating a subset from which alpha diversity can be calculated.
    Note that removal of data, while useful for diversity assessment, is no longer considered good practice.
    https://doi.org/10.1371/journal.pcbi.1003531
    """
    def __init__(self, n=1000, processes=1, seed=62375, in_dir='./', out_dir=''):
        logger.info(f'Starting instance of {__class__.__name__} in {__name__}.')
        self.in_dir = in_dir
        self.downsize_to = n
        self.seed = seed

        if out_dir:
            self.out_dir = out_dir
        else:
            self.out_dir = in_dir + 'downsized/'
            logger.info('No output directory specified.\nCreating default: %s ' % self.out_dir)
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        else:
            logger.info('Looks like %s exists.' % self.out_dir)

        if processes > mp.cpu_count():
            logger.warning('Number of requested processes [%i] is greater that the number of system CPUs [%i]' %
                           (processes, mp.cpu_count()))
            self.processes = mp.cpu_count()
            logger.info('Number of processes set to maximum number of detected CPUs [%i].' % self.processes)
        else:
            self.processes = processes
            logger.info('Number of processes set to %i ' % self.processes)

    def detect_seq_reads(self):
        """The provided directory might contain files that are not sequence reads.
        This module attempts to identify ONLY files with sequence data."""
        logger.info(f'Detecting sequence files in {self.in_dir}')
        seq_files = [f for f in os.listdir(self.in_dir) if os.path.isfile(self.in_dir + f)]
        logger.info(f'List of Sequencing Files Found:\n{seq_files}')
        return seq_files

    def read_info(self, seq_file):
        """Gather information on the number of reads for each of the sequence reads in the given directory.
        Number of reads can be plotted or reported. Files names and headers are used to match pairs.
        Both compressed and gzipped files are accepted.

        Function returns a dictionary where the filename is the key and the value is a list of information:
        {file: [header, F/R, Number of reads, format]}
        # e.g. {'forward.fastq': ['@A00709:44:HYG57DSXX:2:1101:10737:1266', '1', 100257, 'Uncompressed/FASTQ']"""

        seq_files_info = {}
        detect_type = filetype.guess(self.in_dir + seq_file)
        if detect_type:
            if detect_type.mime == 'application/gzip':
                logger.info('Detected gzipped file for %s' % seq_file)
                file = gzip.open(self.in_dir + seq_file, 'r')
                read_num = sum(1 for _ in gzip.open(self.in_dir + seq_file)) // 4
                header = file.readline().decode().split(' ')[0]
                file.seek(0, 0)
                pair = file.readline().decode().split(' ')[1][0]
                file.close()
                seq_files_info[seq_file] = [header, pair, read_num, detect_type.mime]
                return seq_files_info

        if not detect_type:
            logger.info('Assuming uncompressed fastq file for %s' % seq_file)
            file = open(self.in_dir + seq_file, 'r')
            read_num = sum(1 for _ in open(self.in_dir + seq_file)) // 4
            header = file.readline().split(' ')[0]
            file.seek(0, 0)
            pair = file.readline().split(' ')[1][0]
            file.close()
            seq_files_info[seq_file] = [header, pair, read_num, 'Uncompressed/FASTQ']
            return seq_files_info

    def down_dir_pair(self):
        files_to_downsize = self.detect_seq_reads()
        logging.info('Found %i files.' % len(files_to_downsize))

        '''This is a quick but efficient multiprocessing implementation to handle retrieving information from files
        in the target directory. The Pool is created, the number of workers = the class 'processes' attribute.
        Results are expected as a dictionary, so the resulting 'list of dictionaries' is converted with handy
        dict comprehension.
        '''
        with mp.Pool(processes=self.processes) as pool:
            results = pool.map(self.read_info, files_to_downsize, chunksize=1)
        file_info_dict = {k: v for result in results for k, v in result.items()}

        list_to_downsize = pair_up(file_info_dict)

        worker_parameters = []
        for k in range(len(list_to_downsize)//2):  # number of files divided by 2: one instance per pair
            worker_parameters.append({'raw_file_f': list_to_downsize[k * 2],
                                      'raw_file_r': list_to_downsize[k * 2 + 1],
                                      'read_info': file_info_dict[list_to_downsize[k * 2]],
                                      'in_dir': self.in_dir, 'out_dir': self.out_dir, 'n': self.downsize_to
                                      })

        with mp.Pool(processes=self.processes) as pool:
            check = pool.map(self.instantiate, worker_parameters, chunksize=1)  # noqa

        logger.info('Done!')

    def instantiate(self, wp):
        logger.info('Instantiating with parameters: %s' % wp)
        instance = DownsizePair(**wp)
        instance.downsize_pair()
