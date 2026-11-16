"""Streamlit rendering of the interface's pages. Everything shown about a tool comes from the catalog and its parser."""
import argparse
import dataclasses
import fnmatch
import importlib.util
import inspect
import json
import os
import pprint
from datetime import datetime
from pathlib import Path
import subprocess

import pandas as pd
import streamlit as st

from moonstone import examples
from moonstone.reads import launch as reads_launch
from moonstone.ui import checks, forms, runs
from moonstone.ui.catalog import METADATA_RULES, TOOLS, Tool, api_objects, resolve
from moonstone.ui.tutorials import TUTORIALS, Tutorial, run_tutorial
from moonstone.validation import Finding, check_file_format, convert_locale, convert_to_utf8, has_errors
from moonstone.validation.format_check import EXCEL_EXTENSIONS

# Plain-language names of the options; the flag itself is shown next to the name, its help text comes from the parser.
LABELS = {
    'classify': 'Variable to predict', 'model': 'Models', 'kernel': 'SVM kernel', 'roc': 'Plot ROC curves',
    'group_column': 'Subject of each sample', 'explore_kernels': 'Compare SVM kernels',
    'feature_selection': 'Select features', 'folds': 'Cross-validation folds', 'seed': 'Random seed',
    'filtering': 'Minimum mean reads', 'permanova': 'Variables to test (none: all)', 'covariates': 'Covariates',
    'strata': 'Blocks to shuffle within', 'taxa_level': 'Taxonomic level', 'select': 'Keep only the samples where',
    'categorical': 'Numeric variables to treat as groups', 'no_read_threshold': 'Skip the read-threshold filter',
    'permutations': 'Permutations', 'threads': 'Parallel tests (0: all cores)', 'correlation_method': 'Method',
    'second_countfile': 'Second count table', 'min_prevalence': 'Minimum prevalence', 'top': 'Keep the top N features',
    'correction': 'Multiple-testing correction', 'alpha': 'Significance threshold', 'associate': 'Variables to test',
    'reference': 'Reference levels', 'standardize': 'Standardize continuous variables', 'ancombc': 'Also run ANCOM-BC',
    'db': 'Kraken2 database', 'host_db': 'Host Kraken2 database', 'remove_host': 'Hosts to remove (none: keep all reads)',
    'single_end': 'Single-end reads', 'match': 'Only the files matching', 'keep_trimmed_reads': 'Keep the trimmed reads',
    'keep_kraken_output': 'Keep the per-read Kraken2 files', 'slurm': 'Submit to Slurm', 'no_trimming': 'Skip trimming',
    'no_bracken': 'Skip Bracken', 'read_len': 'Read length (Bracken)', 'bracken_threshold': 'Bracken threshold',
    'min_reads': 'Warn below this number of reads', 'min_classified_pct': 'Warn below this % classified',
    'cores': 'Cores (without Slurm)', 'dry_run': 'Only list the steps', 'method': 'Method',
}
# Filled by the app with its pages, so that a result can be sent to the page of the next step.
PAGES = {}
# After a run of the tool (key): the output (relative to the run directory) that the next tool takes, and its button.
NEXT_STEP = {'reads': ('reads_table', 'output/mapping/all_samples.bracken.tsv', 'Make a Moonstone table'),
             'reads_table': ('normalize', 'counts_raw.csv', 'Normalize')}
TABLE_COLUMNS = 50
PREVIEW_ROWS = 20
TABLE_ROWS = 1000
MAX_LEVELS = 20    # a column with more distinct values is not offered as COLUMN=VALUE choices
ALERTS = {'error': st.error, 'warning': st.warning, 'info': st.info}
ALERT_ICONS = {'error': ':material/error:', 'warning': ':material/warning:', 'info': ':material/info:'}


def runs_root() -> str:
    """The directory under which run directories are created."""
    return st.session_state.setdefault('runs_root', os.path.abspath('moonstone_runs'))


def keep_form_values():
    """Keep what was typed in a form when another page is opened (Streamlit drops the state of hidden widgets)."""
    for key in list(st.session_state):
        if key.split('.')[0] in (*TOOLS, 'check') and not key.endswith(('.example', '.run', '.view', '.convert')):
            st.session_state[key] = st.session_state[key]


def _stamp(path):
    try:
        status = os.stat(path)
        return status.st_mtime_ns, status.st_size
    except (OSError, TypeError):
        return None


# The three functions below read the input files: cached on the file's modification time and size.
# ``options``: the form values as JSON, for the tools whose checks depend on them (a folder of reads). The results
# expire after two minutes, because a folder can change without its own modification time changing.
@st.cache_data(max_entries=64, ttl=120, show_spinner='Checking the files...')
def _findings(tool_key, countfile, metadata, second_countfile, stamps, options=''):
    return checks.check_inputs(TOOLS[tool_key], countfile, metadata, second_countfile,
                               json.loads(options) if options else None)


@st.cache_data(max_entries=64, show_spinner=False)
def _columns(path, stamp):
    return checks.metadata_columns(path)


@st.cache_data(max_entries=64, ttl=120, show_spinner=False)
def _summary(tool_key, path, stamp, options=''):
    return checks.summary(TOOLS[tool_key], path, json.loads(options) if options else None)


def show_findings(findings: list):
    """One alert per finding, with where it was found."""
    for finding in findings:
        where = [f'column {finding.column!r}' if finding.column is not None else '',
                 f'row {finding.row}' if finding.row is not None else '',
                 f'value {finding.example!r}' if finding.example is not None else '']
        where = ', '.join(part for part in where if part)
        ALERTS[finding.level](finding.message + (f' ({where})' if where else ''), icon=ALERT_ICONS[finding.level])


# ---------------------------------------------------------------------------------------------------------------------
# Home
# ---------------------------------------------------------------------------------------------------------------------
def _card(pages, key, title, text, caption):
    with st.container(border=True, height='stretch'):
        st.markdown(f'#### {title}')
        st.write(text)
        st.caption(caption)
        if st.button('Open', key=f'home.{key}', icon=':material/arrow_forward:') and key in pages:
            st.switch_page(pages[key])


def render_home(pages: dict):
    """One card per tool, by group. ``pages`` maps a card key to its st.Page."""
    st.title('Moonstone')
    st.markdown(
        'Analysis of metagenomics counts, from sequencing reads or a count table to statistics and plots. '
        'Every page builds a `moonstone` command, checks your files before it runs, and keeps each run in its '
        'own directory with the command that produced it.')
    starts = [('tutorials', ':material/school: Tutorials', 'Worked examples for every tool and for the Python '
               'functions, run on the example files.', 'Start here after a long break.'),
              ('runs', ':material/history: Previous runs', 'Reopen a run: its command, log and results.',
               'Every run keeps the command that produced it.')]
    check = ('check', ':material/fact_check: Check my files', 'Is my file in the format Moonstone expects? '
             'Repair files saved by Excel in a French locale.', 'Needs: any count table or metadata file.')
    cards = {'Start': starts, 'From reads': [], 'Prepare': [check], 'Analyse': []}
    for tool in TOOLS.values():
        needs = 'a count table and a metadata file' if tool.needs_metadata else 'a count table'
        if tool.input_kind == 'reads_dir':
            needs = 'sequencing reads and Kraken2 databases'
        cards[tool.group].append((tool.key, f'{tool.icon} {tool.title}', tool.question,
                                  f'Needs: {needs}. {tool.input_rules[0]}'))
    for group, group_cards in cards.items():
        group_cards = [card for card in group_cards if card[0] in pages]
        if group_cards:
            st.subheader(group)
            columns = st.columns(3)
            for position, card in enumerate(group_cards):
                with columns[position % 3]:
                    _card(pages, *card)


# ---------------------------------------------------------------------------------------------------------------------
# Tool page
# ---------------------------------------------------------------------------------------------------------------------
def _label(action) -> str:
    text = LABELS.get(action.dest, action.dest.replace('_', ' ').capitalize())
    return f'{text} `{forms.flag(action)}`'


def _keep_valid(key, options, multi):
    """Drop a remembered choice that is no longer offered (for example after another metadata file is loaded)."""
    value = st.session_state.get(key)
    if multi:
        st.session_state[key] = [item for item in (value or []) if item in options]
    elif value not in options:
        st.session_state[key] = None


def _widget(tool: Tool, action, columns: dict):
    """The widget of one parser option; returns its value in the form forms.to_argv expects."""
    key = f'{tool.key}.{action.dest}'
    label, help_text = _label(action), action.help
    multi = action.nargs in ('+', '*') or isinstance(action, argparse._AppendAction)
    if isinstance(action, argparse._StoreTrueAction):
        return st.checkbox(label, key=key, help=help_text)
    if action.choices:
        options = list(action.choices)
        st.session_state.setdefault(key, list(action.default or []) if multi else action.default)
        _keep_valid(key, options, multi)
        shown = tool.choice_labels.get(action.dest, {})
        return (st.multiselect if multi else st.selectbox)(label, options, key=key, help=help_text,
                                                           format_func=lambda choice: shown.get(choice, choice))
    metavar = action.metavar or ''
    if metavar in ('variable', 'column') or '=' in metavar:
        if not columns:    # no readable metadata: free text, several values separated by spaces
            text = st.text_input(label, key=f'{key}.text', help=help_text,
                                 placeholder=metavar + (' ...' if multi else ''))
            return text.split() if multi else (text or None)
        options = list(columns)
        if '=' in metavar:
            options = [f'{column}={value}' for column, values in columns.items() if len(values) <= MAX_LEVELS
                       for value in values]
        _keep_valid(key, options, multi)
        if multi:
            return st.multiselect(label, options, key=key, help=help_text, placeholder='Choose from the metadata')
        return st.selectbox(label, options, index=None, key=key, help=help_text, placeholder='None')
    if action.type in (int, float):
        if action.default is None:
            return st.number_input(label, value=None, step=1 if action.type is int else None, key=key,
                                   help=help_text, placeholder='Default')
        st.session_state.setdefault(key, action.default)
        return st.number_input(label, step=1 if action.type is int else 0.01, key=key, help=help_text)
    if action.default:      # for example a database path from the site file
        st.session_state.setdefault(key, action.default)
    return st.text_input(label, key=key, help=help_text) or None


def _load_example(tool: Tool):
    st.session_state[f'{tool.key}.countfile'] = examples.path(tool.example['countfile'])
    if tool.example['metadata']:
        st.session_state[f'{tool.key}.metadata'] = examples.path(tool.example['metadata'])
    for dest, value in tool.example['values'].items():
        st.session_state[f'{tool.key}.{dest}'] = value


def _path_input(tool: Tool, name: str, label: str, help_text: str):
    path = st.text_input(label, key=f'{tool.key}.{name}', help=help_text, placeholder='/path/to/file').strip()
    return path or None


def _run_view(tool: Tool):
    key = tool.key
    missing = [how for module, how in tool.requires if importlib.util.find_spec(module) is None]
    for how in missing:
        st.warning(f'This tool cannot run here: it needs {how}.', icon=':material/warning:')
    if tool.example:
        st.button('Use the example files', key=f'{key}.example', on_click=_load_example, args=(tool,),
                  icon=':material/science:', help='Fill this form with example files and options that work.')
    countfile = _path_input(tool, 'countfile', tool.input_label, ' '.join(tool.input_rules))
    summary = st.empty()
    metadata = _path_input(tool, 'metadata', 'Metadata', ' '.join(METADATA_RULES)) if tool.needs_metadata else None
    columns = _columns(metadata, _stamp(metadata)) if metadata else {}

    shown, advanced = forms.form_actions(tool)
    values = {}
    for actions, container in [(shown, lambda: st.container(border=True)), (advanced, lambda: st.expander('More options'))]:
        if actions:       # a tool without options gets no empty box
            grid = container().columns(2)
            for position, action in enumerate(actions):
                with grid[position % 2]:
                    values[action.dest] = _widget(tool, action, columns)

    # a folder of reads is checked with the options of the form (databases, hosts, pairing)
    options = json.dumps(values, sort_keys=True, default=str) if tool.input_kind == 'reads_dir' else ''
    if countfile and _summary(key, countfile, _stamp(countfile), options):
        summary.caption(_summary(key, countfile, _stamp(countfile), options))
    ready = bool(countfile) and (bool(metadata) or not tool.needs_metadata) and not missing
    second = values.get('second_countfile')
    findings = []
    if countfile and (metadata or not tool.needs_metadata):
        stamps = (_stamp(countfile), _stamp(metadata), _stamp(second))
        findings = list(_findings(key, countfile, metadata, second, stamps, options))
    planned = Path(runs_root()) / f'{datetime.now():%Y%m%d-%H%M%S}-{key}'
    argv = forms.to_argv(tool, countfile or tool.input_label.upper().replace(' ', '_'), metadata or 'METADATA',
                         str(planned), values)
    command_error = forms.check_command(tool, argv)
    if command_error:
        findings.append(Finding('error', 'command', command_error))
    show_findings(findings)
    if ready and not findings:
        st.success('The files and options are valid.', icon=':material/check_circle:')

    st.caption('The command this form runs (a new run directory is created for each run):')
    st.code(forms.preview(tool, argv), language='bash', wrap_lines=True)
    if st.button('Run', key=f'{key}.run', type='primary', icon=':material/play_arrow:',
                 disabled=not ready or has_errors(findings)):
        run_dir = runs.new_run_dir(runs_root(), key)
        argv = forms.to_argv(tool, countfile, metadata, str(run_dir), values)
        with st.status(f'Running {tool.title}...', expanded=True) as status:
            output, lines = st.empty(), []
            for line in runs.run(tool, argv, run_dir):
                lines.append(line)
                output.code(''.join(lines[-25:]), language=None)
            failed = (run_dir / 'exit_code.txt').read_text().strip() != '0'
            submitted = (run_dir / reads_launch.JOB_ID_FILE).exists()
            done = 'submitted to Slurm' if submitted else 'finished'
            status.update(label=f'{tool.title} failed' if failed else f'{tool.title} {done}',
                          state='error' if failed else 'complete', expanded=False)
        st.session_state[f'{key}.last_run'] = str(run_dir)
    if st.session_state.get(f'{key}.last_run'):
        render_results(Path(st.session_state[f'{key}.last_run']))


def _preview_table(name: str):
    path = examples.path(name)
    table = pd.read_csv(path, sep='\t' if name.endswith('.tsv') else ',', nrows=PREVIEW_ROWS, dtype=str)
    st.dataframe(table, hide_index=True)
    with open(path, 'rb') as handle:
        st.download_button(f'Download {name}', handle.read(), file_name=name, key=f'example.{name}',
                           icon=':material/download:', on_click='ignore')


def _inputs_view(tool: Tool):
    st.subheader(tool.input_label)
    st.markdown('\n'.join(f'- {rule}' for rule in tool.input_rules))
    if tool.example:
        st.caption(f'Example: `{tool.example["countfile"]}` (simulated values, first {PREVIEW_ROWS} rows)')
        _preview_table(tool.example['countfile'])
    if tool.needs_metadata:
        st.subheader('Metadata')
        st.markdown('\n'.join(f'- {rule}' for rule in METADATA_RULES))
        st.caption(f'Example: `{tool.example["metadata"]}` (simulated values, first {PREVIEW_ROWS} rows)')
        _preview_table(tool.example['metadata'])
    st.subheader('What a run writes')
    st.table(pd.DataFrame({'File': list(tool.outputs), 'Content': list(tool.outputs.values())}), hide_index=True)


def _documentation_view(tool: Tool):
    shown, advanced = forms.form_actions(tool)
    if shown or advanced:
        st.subheader('Options')
        rows = [{'Option': f'{forms.flag(action)} {action.metavar or ""}'.strip(),
                 'Default': '' if action.default in (None, False, []) else str(action.default),
                 'What it does': action.help} for action in shown + advanced]
        st.table(pd.DataFrame(rows), hide_index=True)
    st.subheader('Reference')
    for dotted in tool.docs:
        with st.expander(dotted):
            show_docstring(resolve(dotted))
    st.link_button('User guide', tool.guide_url, icon=':material/menu_book:')


def show_docstring(obj):
    """An object's signature and docstring."""
    if not inspect.ismodule(obj):
        try:
            st.code(f'{obj.__name__}{inspect.signature(obj)}', language='python', wrap_lines=True)
        except (TypeError, ValueError):
            pass
    # shortcut: docstrings are shown as plain text, in three mixed styles; render with docutils if that reads badly.
    st.code(inspect.getdoc(obj) or 'No documentation.', language=None, wrap_lines=True)


def _tutorial_view(tool: Tool):
    for tutorial in TUTORIALS.values():
        if tutorial.tool == tool.key:
            render_tutorial(tutorial)


TOOL_VIEWS = {'Run': _run_view, 'Tutorial': _tutorial_view, 'Inputs and examples': _inputs_view,
              'Documentation': _documentation_view}


def render_tool(tool: Tool):
    """The page of one catalog tool."""
    st.title(f'{tool.icon} {tool.title}')
    st.markdown(f'**{tool.question}**')
    # on_change='rerun': only the open tab is computed
    tabs = st.tabs(list(TOOL_VIEWS), key=f'{tool.key}.view', on_change='rerun')
    for tab, view in zip(tabs, TOOL_VIEWS.values()):
        if tab.open:
            with tab:
                view(tool)


# ---------------------------------------------------------------------------------------------------------------------
# Results of a run
# ---------------------------------------------------------------------------------------------------------------------
def _show_table(path: Path, separator: str):
    """A preview of a table: at most TABLE_ROWS rows and TABLE_COLUMNS columns are read and sent to the browser."""
    with open(path, 'rb') as handle:
        first = handle.readline()
        lines_after_first = sum(1 for _ in handle)
    skip = 1 if first.startswith(b'# ') else 0       # the tables of the reads workflow start with a comment line
    n_columns = len(pd.read_csv(path, sep=separator, skiprows=skip, nrows=0).columns)
    table = pd.read_csv(path, sep=separator, skiprows=skip, nrows=TABLE_ROWS,
                        usecols=range(min(n_columns, TABLE_COLUMNS)))
    st.dataframe(table, hide_index=True)
    n_rows, shown = lines_after_first - skip, []
    if n_rows > TABLE_ROWS:
        shown.append(f'first {TABLE_ROWS:,} of {n_rows:,} rows')
    if n_columns > TABLE_COLUMNS:
        shown.append(f'first {TABLE_COLUMNS} of {n_columns:,} columns')
    if shown:
        st.caption(f'Shown: {" and ".join(shown)}; download the file for all of it.')


def _show_file(path: Path):
    suffix = path.suffix.lower()
    if suffix in ('.csv', '.tsv'):
        _show_table(path, '\t' if suffix == '.tsv' else ',')
    elif suffix in ('.png', '.jpg', '.jpeg', '.svg'):
        st.image(str(path))
    elif suffix == '.pdf':
        try:
            st.pdf(str(path), height=600, key=f'pdf.{path}')
        except Exception:    # the PDF viewer is an optional Streamlit component
            st.caption('No PDF viewer installed (pip install "streamlit[pdf]"); download the file to see it.')
    elif suffix == '.md':
        st.markdown(path.read_text(errors='replace'))
    elif suffix in ('.txt', '.log', '.yaml', '.yml'):
        st.code(path.read_text(errors='replace'), language='yaml' if suffix in ('.yaml', '.yml') else None,
                wrap_lines=True)
    else:
        st.caption('No preview for this type of file.')
    st.download_button(f'Download {path.name}', path.read_bytes(), file_name=path.name, key=f'download.{path}',
                       icon=':material/download:', on_click='ignore')


def _slurm_status(run_dir: Path) -> str:
    """The state of a run submitted to Slurm, refreshed while it is pending or running. Returns the state."""
    active = reads_launch.status(str(run_dir)).state in ('pending', 'running')

    def block():
        status = reads_launch.status(str(run_dir))      # read from the run directory and Slurm: survives a restart
        if active and status.state not in ('pending', 'running'):
            st.rerun()                                  # the run just ended: show its results
        steps = f' · {status.done} of {status.total} steps done' if status.total else ''
        st.markdown(f'**Status: {status.state}** (Slurm job {status.job_id}){steps}')
        if status.total:
            st.progress(status.done / status.total)
        if status.state in ('pending', 'running'):
            st.caption('Checked again every 30 seconds. The run continues if you close this page or the interface.')
            with st.popover('Cancel the run', icon=':material/cancel:'):
                st.write(f'Cancel Slurm job {status.job_id} and the jobs it started?')
                if st.button('Yes, cancel the run', key=f'cancel.{run_dir.name}'):
                    try:
                        reads_launch.cancel(str(run_dir))
                        st.success('Cancellation requested.', icon=':material/check_circle:')
                    except (ValueError, OSError, subprocess.CalledProcessError) as error:
                        st.error(f'The run could not be cancelled: {error}', icon=':material/error:')
        elif status.state != 'finished' and status.log_tail:
            with st.expander('End of the Snakemake log', expanded=status.state == 'failed'):
                st.code(status.log_tail, language=None, wrap_lines=True)

    st.fragment(block, run_every=30 if active else None)()
    return reads_launch.status(str(run_dir)).state


def _next_step(tool: Tool, run_dir: Path):
    """Offer to send the result of a run to the tool that takes it next."""
    if tool is None or tool.key not in NEXT_STEP:
        return
    target, relative, label = NEXT_STEP[tool.key]
    result = run_dir / relative
    if not result.exists():
        return
    if st.button(label, key=f'next.{run_dir.name}', icon=':material/arrow_forward:',
                 help=f'Open "{TOOLS[target].title}" with this result filled in.'):
        st.session_state[f'{target}.countfile'] = str(result)
        if target in PAGES:
            st.switch_page(PAGES[target])
        st.info(f'`{result}` is filled in on the page "{TOOLS[target].title}".', icon=':material/info:')


def render_results(run_dir: Path):
    """The files of a run directory, the tool's main outputs first, then the log and the command."""
    run_dir = Path(run_dir)
    tool = TOOLS.get(run_dir.name.split('-')[2]) if run_dir.name.count('-') >= 2 else None
    exit_file, log_file, command_file = (run_dir / name for name in ('exit_code.txt', 'run.log', 'command.txt'))
    st.subheader('Results')
    st.caption(f'Run directory: `{run_dir}`')
    log = log_file.read_text(errors='replace') if log_file.exists() else ''
    failed = exit_file.exists() and exit_file.read_text().strip() != '0'
    if failed:
        st.error('The run failed. The end of its log:\n\n```\n' + '\n'.join(log.splitlines()[-8:]) + '\n```',
                 icon=':material/error:')
    elif command_file.exists() and log_file.exists() and not exit_file.exists():
        st.warning('This run did not finish.', icon=':material/warning:')
    if not failed and (run_dir / reads_launch.JOB_ID_FILE).exists():
        _slurm_status(run_dir)
    outputs = tool.outputs if tool else {}
    files = runs.result_files(run_dir, outputs)
    described = {path: text for path in files for pattern, text in outputs.items()
                 if fnmatch.fnmatch(path.relative_to(run_dir).as_posix(), pattern)}
    for pattern in outputs:        # in the order of the catalog: the main result first
        for path in files:
            if path in described and fnmatch.fnmatch(path.relative_to(run_dir).as_posix(), pattern):
                with st.expander(f'**{path.relative_to(run_dir)}**: {described[path]}', expanded=True):
                    _show_file(path)
    _next_step(tool, run_dir)
    if tool is not None and tool.key in ('reads_table', 'normalize'):
        for path in files:
            if path in described and path.suffix == '.csv':
                st.caption('Path of the table, to paste in an analysis page:')
                st.code(str(path), language=None)
    for path in files:
        if path not in described:
            with st.expander(path.name):
                _show_file(path)
    if log:
        with st.expander('Log of the run'):
            st.code(log, language=None, wrap_lines=True)
    if command_file.exists():
        st.caption('Command, Moonstone version and start time:')
        st.code(command_file.read_text(), language='bash', wrap_lines=True)


# ---------------------------------------------------------------------------------------------------------------------
# Check my files
# ---------------------------------------------------------------------------------------------------------------------
LOCALE_CODES = {'semicolon_separator', 'decimal_comma', 'french_thousands', 'ambiguous_comma', 'day_first_dates',
                'month_first_dates', 'ambiguous_dates', 'mixed_date_order'}
SEPARATORS = {';': 'semicolon', ',': 'comma', '\t': 'tab'}
DATE_ORDERS = {None: 'leave dates unchanged', 'dmy': 'day first (dd/mm/yyyy)', 'mdy': 'month first (mm/dd/yyyy)'}


def _conversion(path: str, codes: set):
    """Offer the conversion that fits the findings; the result is always a new file."""
    stem, extension = os.path.splitext(path)
    st.subheader('Repair a copy of the file')
    if 'not_utf8' in codes:
        st.write('The file is not UTF-8. Its encoding cannot be detected reliably: choose the one it was saved in.')
        encoding = st.selectbox('Encoding of the file', ['cp1252', 'latin-1', 'mac_roman'], key='check.encoding',
                                help='cp1252 is what Excel on Windows writes for Western European languages.')

        def convert(output):
            return convert_to_utf8(path, output, encoding)
    else:
        st.write('The file uses non-English conventions. Declare them; numbers and dates are rewritten with a dot '
                 'as decimal mark and ISO dates, and the columns are separated by commas.')
        grid = st.columns(3)
        separator = grid[0].selectbox('Separator in the file', list(SEPARATORS), format_func=SEPARATORS.get,
                                      index=0 if 'semicolon_separator' in codes else 1, key='check.sep')
        decimal = grid[1].selectbox('Decimal mark in the file', [',', '.'], key='check.decimal')
        guessed = 'dmy' if 'day_first_dates' in codes else 'mdy' if 'month_first_dates' in codes else None
        date_order = grid[2].selectbox('Dates in the file', list(DATE_ORDERS), format_func=DATE_ORDERS.get,
                                       index=list(DATE_ORDERS).index(guessed), key='check.dates')

        def convert(output):
            return convert_locale(path, output, sep=separator, decimal=decimal, date_order=date_order)
    output = st.text_input('Write the converted copy to', value=f'{stem}.moonstone{extension}',
                           help='A new file. The original is never changed, and an existing file is never overwritten.')
    if st.button('Convert', key='check.convert', type='primary', icon=':material/build:'):
        try:
            result = convert(output)
        except (ValueError, OSError) as error:
            st.error(str(error), icon=':material/error:')
            return
        st.success(f'Wrote `{result.output}` ({len(result.changes)} value(s) changed). Use that file from now on.',
                   icon=':material/check_circle:')
        show_findings(result.notes + result.findings)
        if result.changes:
            st.caption('Every value that was changed:')
            st.dataframe(result.changes_frame().head(TABLE_ROWS), hide_index=True)


def render_check():
    """Check one file against the format Moonstone expects, and against the tool that will read it."""
    st.title(':material/fact_check: Check my files')
    st.markdown('Moonstone expects UTF-8 text, comma or tab separators, a dot as decimal mark and ISO dates. '
                'Files saved by Excel in a French locale often differ on every point. Nothing is changed in '
                'your file: a repair, when one is possible, is written to a new file.')
    path = st.text_input('File to check', key='check.file', placeholder='/path/to/file.csv').strip()
    grid = st.columns(2)
    tool_key = grid[0].selectbox('To be used as the count table of', list(TOOLS), index=None, key='check.tool',
                                 format_func=lambda key: TOOLS[key].title, placeholder='Any tool: check the format only',
                                 help='Also checks that this tool can read the table.')
    metadata = grid[1].text_input('Metadata file to match the samples with', key='check.metadata',
                                  placeholder='optional', disabled=tool_key is None).strip()
    if not path:
        return
    if tool_key:
        tool = TOOLS[tool_key]
        tool = dataclasses.replace(tool, needs_metadata=tool.needs_metadata and bool(metadata))
        findings = checks.check_inputs(tool, path, metadata or None)
    elif not os.path.isfile(path):
        findings = [Finding('error', 'file_missing', f'{path}: file not found.')]
    else:
        findings = check_file_format(path)
    show_findings(findings)
    if not findings:
        st.success('Nothing to report: the file is in the expected format.', icon=':material/check_circle:')
    elif not has_errors(findings):
        st.success('No error: the file can be used.', icon=':material/check_circle:')
    codes = {finding.code for finding in findings}
    if os.path.isfile(path) and ('not_utf8' in codes or codes & LOCALE_CODES):
        _conversion(path, codes)
    if os.path.isfile(path) and not has_errors(findings):
        spreadsheet = path.lower().rsplit('.', 1)[-1] in EXCEL_EXTENSIONS
        read = pd.read_excel if spreadsheet else lambda source, **options: pd.read_csv(
            source, sep=None, engine='python', **options)
        with st.expander(f'First {PREVIEW_ROWS} rows'):
            st.dataframe(read(path, nrows=PREVIEW_ROWS, dtype=str), hide_index=True)


# ---------------------------------------------------------------------------------------------------------------------
# Previous runs
# ---------------------------------------------------------------------------------------------------------------------
def render_runs():
    """The runs under the run directory, and the results of the one picked."""
    st.title(':material/history: Previous runs')
    listed = runs.list_runs(runs_root())
    if not listed:
        st.info(f'No run yet in `{runs_root()}`. Change the run directory in the sidebar to look elsewhere.',
                icon=':material/info:')
        return
    status = {None: 'did not finish', 0: 'finished'}

    def state(run):       # a run submitted to Slurm: its state there; otherwise how its command ended
        if run.exit_code == 0 and (run.path / reads_launch.JOB_ID_FILE).exists():
            return reads_launch.status(str(run.path)).state
        return status.get(run.exit_code, 'failed')

    st.dataframe(pd.DataFrame([{
        'Run': run.path.name, 'Tool': TOOLS[run.tool_key].title if run.tool_key in TOOLS else run.tool_key,
        'Status': state(run), 'Command': run.command} for run in listed]), hide_index=True)
    picked = st.selectbox('Open a run', [run.path.name for run in listed], key='runs.pick')
    render_results(Path(runs_root()) / picked)


# ---------------------------------------------------------------------------------------------------------------------
# Python API
# ---------------------------------------------------------------------------------------------------------------------
def _first_line(obj) -> str:
    return (inspect.getdoc(obj) or '').strip().split('\n')[0]


def render_api():
    """Every public class and function of the Python API, with its signature and docstring."""
    st.title(':material/menu_book: Python API')
    st.markdown('The classes and functions behind the tools, for use in a script or a notebook. '
                'The Tutorials page shows them at work on the example files.')
    query = st.text_input('Search', key='api.search', placeholder='A name, or a word of the documentation',
                          icon=':material/search:').strip().lower()
    for group, objects in api_objects().items():
        matches = [(name, obj) for name, obj in objects
                   if query in name.lower() or query in (inspect.getdoc(obj) or '').lower()]
        if not matches:
            continue
        st.subheader(group)
        for name, obj in matches:
            with st.expander(f'**{obj.__name__}** {_first_line(obj)}'):
                st.code(f'from {name.rpartition(".")[0]} import {obj.__name__}', language='python')
                show_docstring(obj)
                methods = [(method_name, method) for method_name, method in inspect.getmembers(obj, inspect.isfunction)
                           if not method_name.startswith('_')] if inspect.isclass(obj) else []
                if methods:
                    st.markdown('**Methods**\n\n' + '\n'.join(
                        f'- `{method_name}`: {_first_line(method)}' for method_name, method in methods))


# ---------------------------------------------------------------------------------------------------------------------
# Tutorials
# ---------------------------------------------------------------------------------------------------------------------
# cache_resource: results hold figures and tables that are shown, never modified
@st.cache_resource(max_entries=32, show_spinner='Running the tutorial on the example files...')
def _tutorial_results(key: str) -> list:
    return run_tutorial(TUTORIALS[key])


def _show_value(value, key: str):
    """The result of a tutorial step: a table, a figure, an image file, or text."""
    if value is None:
        return
    if isinstance(value, pd.Series):
        value = value.to_frame()
    module = type(value).__module__
    if isinstance(value, pd.DataFrame):
        st.dataframe(value)
    elif module.startswith('plotly'):
        st.plotly_chart(value, key=key)
    elif module.startswith('matplotlib'):
        st.pyplot(value)
    elif isinstance(value, str) and value.lower().endswith(('.png', '.jpg')) and os.path.isfile(value):
        st.image(value)
    else:
        st.code(pprint.pformat(value, width=110), language=None, wrap_lines=True)


def render_tutorial(tutorial: Tutorial):
    """A tutorial: each step's text, its code, and the result of running that code on the example files."""
    st.header(tutorial.title)
    st.caption(tutorial.summary)
    try:
        results = _tutorial_results(tutorial.key)
    except RuntimeError as error:
        st.error(f'{error}. The code is shown without its results.', icon=':material/error:')
        results = [(None, '')] * len(tutorial.steps)
    for number, (step, (value, printed)) in enumerate(zip(tutorial.steps, results), start=1):
        st.markdown(step.text)
        if step.code:
            st.code(step.code, language='python', wrap_lines=True)    # nothing hidden off the right edge
            if printed:
                st.code(printed, language=None, wrap_lines=True)
            _show_value(value, key=f'tutorial.{tutorial.key}.{number}')
        st.space('small')


def _open_tutorial(group: str):
    st.session_state['tutorials.pick'] = st.session_state[f'tutorials.menu.{group}']


def render_tutorials():
    """The list of tutorials by group, and the one picked."""
    st.title(':material/school: Tutorials')
    st.markdown('Worked examples for every tool and for the Python functions behind them. Each piece of code is '
                'run on the example files when you open the tutorial, so what you see is what the code does today. '
                'Copy any block into a script or a notebook.')
    picked = st.session_state.setdefault('tutorials.pick', next(iter(TUTORIALS)))
    menu, body = st.columns([1, 4], gap='medium')
    with menu:
        groups = dict.fromkeys(tutorial.group for tutorial in TUTORIALS.values())
        for group in groups:
            keys = [key for key, tutorial in TUTORIALS.items() if tutorial.group == group]
            # one list per group; only the group of the open tutorial shows a selection
            st.session_state[f'tutorials.menu.{group}'] = picked if picked in keys else None
            st.radio(group, keys, format_func=lambda key: TUTORIALS[key].title, key=f'tutorials.menu.{group}',
                     on_change=_open_tutorial, args=(group,))
    with body:
        render_tutorial(TUTORIALS[picked])
