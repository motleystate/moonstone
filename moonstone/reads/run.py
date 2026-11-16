"""``moonstone-reads run``: from a folder of sequencing reads to Bracken-corrected counts."""
import argparse
import os
import subprocess
import sys
import tempfile

from moonstone.reads import config, launch, samples
from moonstone.validation import Finding

RUN_FILES = ('config.yaml', 'samples.tsv', 'reads', 'output')


def build_parser() -> argparse.ArgumentParser:
    """Parser of ``moonstone-reads run``. Database defaults come from the site file, if there is one."""
    site = config.site()
    hosts = ', '.join(f'{name} ({taxid})' for name, taxid in config.HOSTS.items())
    parser = argparse.ArgumentParser(
        prog='moonstone-reads run',
        description='Trim adapters (Trimmomatic), remove host reads, classify with Kraken2 and correct with Bracken. '
                    'Writes the Bracken-corrected raw counts of all samples, and a run report.')
    parser.add_argument('reads_dir', help='Folder with the read files (sub-folders are searched). Paired files end '
                                          'with _1/_2, _R1/_R2 or _R1_001/_R2_001 before .fastq.gz or .fq.gz')
    parser.add_argument('outdir', help='Run directory, created by this command; it must not already hold a run')
    parser.add_argument('--db', metavar='DIR', default=site.get('db'),
                        help='Kraken2 database used to classify the reads (with its Bracken k-mer file)')
    parser.add_argument('--host-db', metavar='DIR', default=site.get('host_db'),
                        help='Kraken2 database of the host genomes; it must contain the hosts chosen with --remove-host')
    parser.add_argument('--remove-host', nargs='*', choices=list(config.HOSTS), default=list(config.HOSTS),
                        help=f'Hosts whose reads are removed before classification: {hosts}. Give no name to skip '
                             'host removal. Only reads classified exactly as the host are removed '
                             '(default: human mouse)')
    parser.add_argument('--single-end', action='store_true', help='The reads are single-end (default: paired-end)')
    parser.add_argument('--match', metavar='REGEX',
                        help='Use only the read files whose name matches this regular expression')
    parser.add_argument('--keep-trimmed-reads', action='store_true',
                        help='Keep the trimmed reads. This allows mapping with another database without trimming '
                             'again (default: deleted once no step needs them)')
    parser.add_argument('--keep-kraken-output', action='store_true',
                        help='Keep the per-read Kraken2 files (large). The per-sample reports are always kept '
                             '(default: deleted once no step needs them)')
    parser.add_argument('--slurm', action='store_true',
                        help='Submit the run to Slurm and return at once; follow it with: moonstone-reads status OUTDIR')
    parser.add_argument('--no-trimming', action='store_true', help='Skip adapter trimming')
    parser.add_argument('--no-bracken', action='store_true', help='Skip Bracken: Kraken2 counts only')
    parser.add_argument('--read-len', type=int, default=150,
                        help='Read length for Bracken; the database needs the matching k-mer file (default: 150)')
    parser.add_argument('--bracken-threshold', type=int, default=5,
                        help='Bracken: minimum number of reads for a species to be re-estimated (default: 5)')
    parser.add_argument('--min-reads', type=int, default=100000,
                        help='Report: warn when fewer reads reach the classification (default: 100000)')
    parser.add_argument('--min-classified-pct', type=float, default=20.0,
                        help='Report: warn when less than this percentage of reads is classified (default: 20)')
    parser.add_argument('--cores', type=int, default=4, help='Without --slurm: number of cores to use (default: 4)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Write the run directory and list the steps, without running or submitting them')
    return parser


def check_inputs(reads_dir: str, options: dict) -> tuple:
    """
    Everything `moonstone-reads run` checks before it writes anything: the read files (found, paired, named),
    then the databases, adapters and programs. Nothing is written.

    :param options: the options of the command by name (``db``, ``host_db``, ``remove_host``, ``single_end``,
        ``match`` ...); missing ones take the command's defaults.
    :return: (samples found, findings); the run can start when there is no finding.
    """
    os.environ['PATH'] = launch.environment()['PATH']       # find the tools of this Python's environment
    paired = not options.get('single_end')
    found, findings = samples.discover(reads_dir, paired=paired, match=options.get('match') or None)
    if findings:
        return found, findings
    if not options.get('db'):
        return found, [Finding('error', 'no_database', '--db is required: give the Kraken2 database used to classify '
                                                       'the reads (there is no default in the site file).')]
    options = dict(options, adapter_dir=options.get('adapter_dir') or config.site().get('adapter_dir'))
    with tempfile.TemporaryDirectory() as scratch:
        samples.link(found, os.path.join(scratch, 'reads'))
        trial = config.build(os.path.join(scratch, 'reads'), scratch, paired, samples.extension(found[0].read1), options)
        return found, config.validate(trial)


def _fail(messages) -> int:
    for message in messages:
        print(f'Error: {message}', file=sys.stderr)
    return 1


def main(argv=None) -> int:
    """Check the inputs, write the run directory, then run the workflow or submit it to Slurm."""
    args = build_parser().parse_args(argv)
    site = config.site()
    outdir = os.path.abspath(args.outdir)
    if any(os.path.exists(os.path.join(outdir, name)) for name in RUN_FILES):
        return _fail([f'{outdir} already holds a run; choose another directory. Nothing is overwritten.'])

    # everything is checked before anything is written
    options = dict(vars(args), adapter_dir=site.get('adapter_dir'))
    found, findings = check_inputs(args.reads_dir, options)
    if findings:
        return _fail(finding.message for finding in findings)
    extension = samples.extension(found[0].read1)

    os.makedirs(outdir, exist_ok=True)
    samples.link(found, os.path.join(outdir, 'reads'))
    samples.write_sheet(found, os.path.join(outdir, 'samples.tsv'))
    config.write(config.build(os.path.join(outdir, 'reads'), os.path.join(outdir, 'output'), not args.single_end,
                              extension, options), os.path.join(outdir, 'config.yaml'))
    print(f'{len(found)} samples ({"single" if args.single_end else "paired"}-end). Run directory: {outdir}')

    if args.slurm and not args.dry_run:
        command = launch.snakemake_command(outdir, profile=site.get('profile') or launch.SLURM_PROFILE)
        try:
            job_id = launch.submit(outdir, command, submit=site.get('submit', 'sbatch'),
                                   head_options=site.get('head_options', '-c 1 --mem 4G'))
        except RuntimeError as error:
            return _fail([f'The run was not submitted: {error}'])
        print(f'Submitted to Slurm as job {job_id}. Follow it with: moonstone-reads status {outdir}')
        return 0
    command = launch.snakemake_command(outdir, cores=args.cores, dry_run=args.dry_run)
    return subprocess.call(command, env=launch.environment())


if __name__ == '__main__':
    sys.exit(main())
