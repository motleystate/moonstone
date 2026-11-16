"""Run the reads workflow: build the Snakemake command, submit it to Slurm, follow it, cancel it."""
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Optional

from moonstone.reads import WORKFLOW_DIR

SNAKEFILE = os.path.join(WORKFLOW_DIR, 'Snakefile')
SLURM_PROFILE = os.path.join(WORKFLOW_DIR, 'profiles', 'slurm')
JOB_ID_FILE = 'head_job_id.txt'
REPORT = os.path.join('output', 'pipeline_report.md')
SUBMITTED = re.compile(r'Submitted (?:batch )?job (\d+)')
PROGRESS = re.compile(r'(\d+) of (\d+) steps \(\d+%\) done')
# Slurm job states, as printed by squeue and sacct
RUNNING_STATES = {'RUNNING', 'COMPLETING', 'CONFIGURING', 'SUSPENDED', 'RESIZING'}
PENDING_STATES = {'PENDING', 'REQUEUED', 'RESV_DEL_HOLD', 'REQUEUE_HOLD', 'REQUEUE_FED'}
LOG_TAIL_LINES = 30


@dataclass(frozen=True)
class Status:
    """Where a run stands. ``state`` is one of: not started, pending, running, finished, failed, cancelled, unknown."""
    state: str
    job_id: Optional[str]
    done: int           # workflow steps finished
    total: int
    log_tail: str


def environment() -> dict:
    """
    The environment for Snakemake and Slurm commands: the folder of this Python is added to PATH (so the tools of
    its environment are found even when it is not activated), and memory settings inherited from an enclosing
    Slurm job are removed, because they conflict with the memory requests of the jobs submitted from it.
    """
    env = {key: value for key, value in os.environ.items()
           if key not in ('SLURM_MEM_PER_CPU', 'SLURM_MEM_PER_GPU', 'SLURM_MEM_PER_NODE')}
    here = os.path.dirname(os.path.abspath(sys.executable))
    paths = env.get('PATH', '').split(os.pathsep)
    if here not in paths:
        env['PATH'] = os.pathsep.join(paths + [here])
    return env


def snakemake_command(run_dir: str, cores: Optional[int] = None, profile: Optional[str] = None,
                      dry_run: bool = False) -> list:
    """The Snakemake command of a run directory: with ``profile`` for Slurm, else on ``cores`` local cores."""
    command = [sys.executable, '-m', 'snakemake', '-s', SNAKEFILE, '--configfile',
               os.path.join(run_dir, 'config.yaml'), '--directory', run_dir]
    command += ['--profile', profile] if profile else ['-c', str(cores or 1)]
    return command + (['-n'] if dry_run else [])


def submit(run_dir: str, command: list, submit: str = 'sbatch', head_options: str = '-c 1 --mem 4G') -> str:
    """
    Submit ``command`` as a Slurm job (the "head job", which then submits one job per step) and return its ID.

    The ID is written to head_job_id.txt as soon as the submission command prints it, even if that command then
    fails on something else (a notification, for example): the job exists and must not be submitted again.

    :param submit: the submission command, for example ``sbatch`` or ``sj-submit --project reads --``.
    :param head_options: its options for the head job, for example ``-p atm -c 1 --mem 4G``.
    :raise RuntimeError: no job was submitted; the message is the output of the submission command.
    """
    arguments = shlex.split(submit) + shlex.split(head_options) + ['--wrap', shlex.join(command)]
    try:
        result = subprocess.run(arguments, cwd=run_dir, env=environment(), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
    except OSError as error:
        raise RuntimeError(f'Could not run {arguments[0]}: {error}') from error
    found = SUBMITTED.search(result.stdout)
    if not found:
        raise RuntimeError(result.stdout.strip() or f'{arguments[0]} printed no job ID (exit code {result.returncode}).')
    with open(os.path.join(run_dir, JOB_ID_FILE), 'w') as handle:
        handle.write(found.group(1) + '\n')
    return found.group(1)


def _job_id(run_dir: str) -> Optional[str]:
    try:
        with open(os.path.join(run_dir, JOB_ID_FILE)) as handle:
            return handle.read().strip() or None
    except OSError:
        return None


def _slurm(*arguments) -> str:
    """The first word printed by a Slurm query, or '' (also when the command is absent or fails)."""
    if not shutil.which(arguments[0]):
        return ''
    result = subprocess.run(arguments, capture_output=True, text=True, env=environment())
    return (result.stdout.split() or [''])[0] if result.returncode == 0 else ''


def _slurm_state(job_id: str) -> str:
    queued = _slurm('squeue', '-h', '-j', job_id, '-o', '%T')
    if queued in PENDING_STATES:
        return 'pending'
    if queued in RUNNING_STATES:
        return 'running'
    ended = _slurm('sacct', '-n', '-X', '-j', job_id, '-o', 'State')
    if ended.startswith('CANCELLED'):
        return 'cancelled'
    if ended in PENDING_STATES:
        return 'pending'
    if ended in RUNNING_STATES:
        return 'running'
    # COMPLETED without the report: Snakemake ended without finishing the workflow
    return 'failed' if ended else 'unknown'


def status(run_dir: str) -> Status:
    """
    The state of a run, read from its directory and from Slurm. It needs nothing from the process that started
    the run, so it gives the same answer after the interface or the terminal was closed.
    """
    done, total, tail = 0, 0, ''
    logs = sorted(glob.glob(os.path.join(run_dir, '.snakemake', 'log', '*.log')))
    if logs:
        with open(logs[-1], errors='replace') as handle:
            lines = handle.readlines()
        tail = ''.join(lines[-LOG_TAIL_LINES:])
        progress = PROGRESS.findall(''.join(lines))
        if progress:
            done, total = map(int, progress[-1])
    job_id = _job_id(run_dir)
    if os.path.exists(os.path.join(run_dir, REPORT)):
        state = 'finished'
    elif job_id is None:
        state = 'unknown' if logs else 'not started'
    else:
        state = _slurm_state(job_id)
    return Status(state, job_id, done, total, tail)


def cancel(run_dir: str) -> None:
    """Cancel the head job of a run; Snakemake then cancels the jobs it submitted."""
    job_id = _job_id(run_dir)
    if job_id is None:
        raise ValueError(f'{run_dir}: no Slurm job was submitted for this run.')
    subprocess.run(['scancel', job_id], check=True, env=environment())
