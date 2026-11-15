"""From form values to a command line, checked by the tool's own parser."""
import contextlib
import io
import re
import shlex
from typing import Optional

from moonstone.ui.catalog import Tool


def _actions(tool: Tool) -> dict:
    return {action.dest: action for action in tool.build_parser()._actions}


def form_actions(tool: Tool) -> tuple:
    """The parser actions of the tool's form: (shown, under "More options"), in catalog order."""
    by_dest = _actions(tool)
    return tuple(by_dest[d] for d in tool.options), tuple(by_dest[d] for d in tool.advanced)


def flag(action) -> str:
    """The long form of an option, e.g. --group_column."""
    return max(action.option_strings, key=len)


def _empty(value) -> bool:
    return value is None or value is False or (isinstance(value, (str, list, tuple)) and len(value) == 0)


def to_argv(tool: Tool, countfile: str, metadata: Optional[str], outdir: str, values: dict) -> list:
    """
    The arguments of the command: file arguments first, then the tool's fixed flags, then one entry per form
    option. A value that is empty or equal to the parser default sends nothing, except for a required option,
    whose flag is always sent.
    """
    argv = [countfile]
    if tool.needs_metadata:
        argv.append(metadata or '')
    if tool.outdir_flag is None:
        argv.append(outdir)
    argv += tool.fixed
    by_dest = _actions(tool)
    for dest in tool.options + tool.advanced:
        action, value = by_dest[dest], values.get(dest)
        unset = _empty(value) or value == action.default
        if unset and dest not in tool.required:
            continue
        argv.append(flag(action))
        if _empty(value) or value is True:
            continue
        argv += [str(v) for v in value] if isinstance(value, (list, tuple)) else [str(value)]
    if tool.outdir_flag is not None:
        argv += [tool.outdir_flag, outdir]
    return argv


def check_command(tool: Tool, argv: list) -> Optional[str]:
    """The parser's own error message for this command, or None when the command is valid."""
    by_dest = _actions(tool)
    for dest in tool.required:
        action = by_dest[dest]
        position = argv.index(flag(action)) + 1
        if action.nargs in (None, '+') and (position == len(argv) or argv[position].startswith('-')):
            return f'{flag(action)} is required.'
    stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(stderr), contextlib.redirect_stdout(io.StringIO()):
            tool.parse(list(argv))
    except SystemExit as stop:
        if stop.code:
            lines = stderr.getvalue().strip().splitlines()
            return re.sub(r'^.*?: error: ', '', lines[-1]) if lines else 'The command is not valid.'
    return None


def preview(tool: Tool, argv: list) -> str:
    """The command as it would be typed in a terminal."""
    program = 'moonstone' if tool.module == 'moonstone.main' else f'python -m {tool.module}'
    return f'{program} {shlex.join(argv)}'
