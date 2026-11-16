"""Check a tool's input files before a run, with the readers the run itself will use."""
import dataclasses
import os
import tempfile
from typing import Optional

from moonstone.parsers import handleMetadata
from moonstone.reads import run as reads_run
from moonstone.ui.catalog import Tool
from moonstone.validation import Finding, check_file_format, has_errors, match_samples


def _read_metadata(path: str):
    """The metadata as the command line reads it (its copy of the file goes to a temporary directory)."""
    with tempfile.TemporaryDirectory() as scratch:
        return handleMetadata.Inputs(path).openmeta(scratch)


def _named(path: str, findings: list) -> list:
    name = os.path.basename(path)
    return [dataclasses.replace(f, message=f'{name}: {f.message}') for f in findings]


def _check_table(path: str, sep: Optional[str], read_samples) -> tuple:
    """Findings for one file, and its sample IDs (None when the file cannot be read)."""
    if not os.path.isfile(path):
        return [Finding('error', 'file_missing', f'{path}: file not found.')], None
    findings = check_file_format(path, sep=sep)
    if has_errors(findings):
        return _named(path, findings), None
    try:
        samples = [str(sample) for sample in read_samples(path)]
    except Exception as error:  # whatever stops the tool's reader stops the run: report it, do not crash the page
        message = str(error) if os.path.basename(path) in str(error) else f'{os.path.basename(path)}: {error}'
        return _named(path, findings) + [Finding('error', 'not_readable', message)], None
    return _named(path, findings), samples


def check_inputs(tool: Tool, countfile: str, metadata: Optional[str] = None,
                 second_countfile: Optional[str] = None, values: Optional[dict] = None) -> list:
    """
    Everything that would stop or distort a run of ``tool`` on these files: missing file, file format
    (:func:`moonstone.validation.check_file_format`), table not readable by this tool, and sample IDs that do
    not match the metadata (:func:`moonstone.validation.match_samples`).

    :param values: the values of the tool's options; used when the input is a folder of reads.
    :return: list of :class:`moonstone.validation.Finding`; a run should not start while one is an error.
    """
    if tool.input_kind == 'reads_dir':       # a folder of reads: the checks of `moonstone-reads run`, with its options
        return reads_run.check_inputs(countfile, values or {})[1]
    findings, tables = [], []
    for path in [countfile] + ([second_countfile] if second_countfile else []):
        found, samples = _check_table(path, tool.count_sep, tool.read_samples)
        findings += found
        if samples is not None:
            tables.append((path, samples))
    if not tool.needs_metadata:
        return findings
    found, metadata_ids = _check_table(metadata or '', ',', lambda path: _read_metadata(path).index)
    findings += found
    if metadata_ids is not None:
        for path, samples in tables:
            findings += _named(path, match_samples(samples, metadata_ids).findings())
    return findings


def metadata_columns(path: str) -> dict:
    """The variables of a metadata file: {column: sorted distinct values as text}; {} if it cannot be read."""
    try:
        metadata = _read_metadata(path)
    except Exception:
        return {}
    return {str(column): sorted(metadata[column].dropna().astype(str).unique()) for column in metadata.columns}


def summary(tool: Tool, countfile: str, values: Optional[dict] = None) -> str:
    """A one-line description of the input, e.g. '60 samples'; '' if the tool cannot read it."""
    try:
        if tool.input_kind == 'reads_dir':
            values = values or {}
            single = bool(values.get('single_end'))
            found = reads_run.samples.discover(countfile, paired=not single, match=values.get('match') or None)[0]
            names = [sample.name for sample in found]
            listed = ', '.join(names[:5]) + (' ...' if len(names) > 5 else '')
            return f"{len(names)} samples, {'single' if single else 'paired'}-end: {listed}" if names else ''
        return f'{len(tool.read_samples(countfile))} samples'
    except Exception:
        return ''
