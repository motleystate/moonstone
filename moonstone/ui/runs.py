"""Run a tool's command into a new directory, and list the runs already made."""
import os
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterator, Optional

from moonstone import __version__
from moonstone.ui.catalog import Tool
from moonstone.ui.forms import preview

BOOKKEEPING = ('command.txt', 'run.log', 'exit_code.txt')


@dataclass(frozen=True)
class Run:
    """A run directory: the tool, the command, and the exit code (None when the run did not finish)."""
    path: Path
    tool_key: str
    command: str
    exit_code: Optional[int]


def new_run_dir(root: str, tool_key: str) -> Path:
    """Create and return <root>/<date>-<time>-<tool_key>; a directory is never reused."""
    base = Path(root) / f'{datetime.now():%Y%m%d-%H%M%S}-{tool_key}'
    candidate, attempt = base, 1
    while True:
        try:
            candidate.mkdir(parents=True)
            return candidate
        except FileExistsError:
            attempt += 1
            candidate = base.with_name(f'{base.name}-{attempt}')


def run(tool: Tool, argv: list, run_dir: Path) -> Iterator[str]:
    """
    Run the command and yield its output line by line. The command is recorded in command.txt before it
    starts, the output in run.log, and the exit code in exit_code.txt when the process ends.
    """
    started = datetime.now().isoformat(timespec='seconds')
    (run_dir / 'command.txt').write_text(f'{preview(tool, argv)}\nmoonstone {__version__}\n{started}\n')
    # unbuffered so that lines arrive as they are printed; Agg so that no plot window is opened
    environment = dict(os.environ, PYTHONUNBUFFERED='1', MPLBACKEND='Agg')
    with open(run_dir / 'run.log', 'w') as log:
        process = subprocess.Popen([sys.executable, '-m', tool.module, *argv], stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True, env=environment)
        for line in process.stdout:
            log.write(line)
            log.flush()
            yield line
        exit_code = process.wait()
    (run_dir / 'exit_code.txt').write_text(f'{exit_code}\n')


def list_runs(root: str) -> list:
    """The runs under ``root``, newest first."""
    runs = []
    for command_file in sorted(Path(root).glob('*/command.txt'), reverse=True):
        path = command_file.parent
        exit_file = path / 'exit_code.txt'
        runs.append(Run(path, path.name.split('-')[2] if path.name.count('-') >= 2 else '',
                        command_file.read_text().splitlines()[0],
                        int(exit_file.read_text()) if exit_file.exists() else None))
    return runs


def result_files(run_dir: Path) -> list:
    """The files a run wrote, without the interface's own bookkeeping files."""
    return sorted(path for path in Path(run_dir).iterdir() if path.is_file() and path.name not in BOOKKEEPING)
