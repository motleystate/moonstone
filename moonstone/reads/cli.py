"""The ``moonstone-reads`` command: run, follow and cancel the reads workflow."""
import subprocess
import sys

from moonstone.reads import launch, run, table

USAGE = '''usage: moonstone-reads COMMAND ...

From sequencing reads to count tables.

commands:
  run READS_DIR OUTDIR ...    trim, remove host reads, classify (Kraken2) and correct (Bracken)
  status OUTDIR               where a run stands
  cancel OUTDIR               cancel a run submitted to Slurm
  table TABLE OUTPUT          a run's table in Moonstone's layout, counts unchanged

Options of a command: moonstone-reads COMMAND --help
Normalization is a separate step: python -m moonstone.normalization.table --help
'''


def _status(arguments) -> int:
    if len(arguments) != 1:
        print('usage: moonstone-reads status OUTDIR', file=sys.stderr)
        return 2
    status = launch.status(arguments[0])
    print(f'{status.state}' + (f' (Slurm job {status.job_id})' if status.job_id else ''))
    if status.total:
        print(f'{status.done} of {status.total} steps done')
    if status.state in ('failed', 'unknown') and status.log_tail:
        print('\nEnd of the Snakemake log:\n' + status.log_tail)
    return 0


def _cancel(arguments) -> int:
    if len(arguments) != 1:
        print('usage: moonstone-reads cancel OUTDIR', file=sys.stderr)
        return 2
    try:
        launch.cancel(arguments[0])
    except ValueError as error:
        print(f'Error: {error}', file=sys.stderr)
        return 1
    except (OSError, subprocess.CalledProcessError) as error:
        print(f'Error: scancel could not cancel the job: {error}', file=sys.stderr)
        return 1
    print('Cancellation requested.')
    return 0


COMMANDS = {'run': run.main, 'status': _status, 'cancel': _cancel, 'table': table.main}


def main(argv=None) -> int:
    arguments = sys.argv[1:] if argv is None else list(argv)
    if not arguments or arguments[0] in ('-h', '--help'):
        print(USAGE)
        return 0
    if arguments[0] not in COMMANDS:
        print(f'moonstone-reads: unknown command {arguments[0]!r}\n\n{USAGE}', file=sys.stderr)
        return 2
    return COMMANDS[arguments[0]](arguments[1:])


if __name__ == '__main__':
    sys.exit(main())
