"""Find read files, pair their mates, and give them the names the reads workflow expects."""
import csv
import os
import re
from dataclasses import dataclass
from typing import Optional

from moonstone.validation import Finding

# (suffix of mate 1, suffix of mate 2), longest first so that _R1_001 is not read as _001
MATE_SUFFIXES = (('_R1_001', '_R2_001'), ('_R1', '_R2'), ('_1', '_2'))
EXTENSIONS = ('fastq.gz', 'fq.gz', 'fastq', 'fq')
SAFE_NAME = re.compile(r'[A-Za-z0-9._-]+')
MAX_LISTED = 5


@dataclass(frozen=True)
class Sample:
    """A sample and its read file(s); ``read2`` is None for single-end reads."""
    name: str
    read1: str
    read2: Optional[str]


def extension(path: str) -> Optional[str]:
    """The read-file extension of a path (fastq.gz, fq.gz, fastq or fq), or None."""
    return next((ext for ext in EXTENSIONS if path.endswith('.' + ext)), None)


def _split_mate(stem: str):
    """(sample name, mate number) of a paired file name without its extension, or None."""
    for mates in MATE_SUFFIXES:
        for number, suffix in enumerate(mates, start=1):
            if stem.endswith(suffix) and len(stem) > len(suffix):
                return stem[:-len(suffix)], number
    return None


def _listing(items) -> str:
    items = sorted(items)
    return ', '.join(items[:MAX_LISTED]) + (f' ... ({len(items)} in all)' if len(items) > MAX_LISTED else '')


def discover(reads_dir: str, paired: bool = True, match: Optional[str] = None) -> tuple:
    """
    Find the read files under ``reads_dir`` (sub-folders included) and group them by sample.

    Paired files are recognised by the endings _1/_2, _R1/_R2 or _R1_001/_R2_001 before the extension.

    :param match: regular expression; only files whose name matches it (``re.search``) are used.
    :return: (list of :class:`Sample` sorted by name, list of findings). The samples can be used only when no
        finding is an error.
    """
    if not os.path.isdir(reads_dir):
        return [], [Finding('error', 'not_a_directory', f'{reads_dir}: not a directory.')]
    files = sorted(os.path.join(folder, name) for folder, _, names in os.walk(reads_dir, followlinks=True)
                   for name in names if extension(name) and (match is None or re.search(match, name)))
    if not files:
        return [], [Finding('error', 'no_reads', f'{reads_dir}: no read file found (expected extensions: '
                                                 f'{", ".join(EXTENSIONS)}' + (f'; names matching {match!r}' if match else '')
                                                 + ').')]
    findings, mates, unpaired = [], {}, []
    for path in files:
        stem = os.path.basename(path)[:-len(extension(path)) - 1]
        key = _split_mate(stem) if paired else (stem, 1)
        if key is None:
            unpaired.append(os.path.basename(path))
        else:
            mates.setdefault(key, []).append(path)
    if unpaired:
        findings.append(Finding('error', 'not_paired', f'{len(unpaired)} file(s) do not end with _1/_2, _R1/_R2 or '
                                                       f'_R1_001/_R2_001: {_listing(unpaired)}. Paired and single-end '
                                                       'reads cannot be mixed in one run.'))
    for (name, number), paths in sorted(mates.items()):
        if len(paths) > 1:
            findings.append(Finding('error', 'duplicate_sample', f'Sample {name} is found more than once: '
                                                                 f'{_listing(paths)}.'))
    names = sorted({name for name, _ in mates})
    unsafe = [name for name in names if not SAFE_NAME.fullmatch(name)]
    if unsafe:
        findings.append(Finding('error', 'sample_name', 'Sample names may only contain letters, digits, ".", "_" '
                                                        f'and "-": {_listing(unsafe)}.'))
    if paired:
        alone = [name for name in names if (name, 1) not in mates or (name, 2) not in mates]
        if alone:
            findings.append(Finding('error', 'mate_missing', f'{len(alone)} sample(s) have one mate without the '
                                                             f'other: {_listing(alone)}.'))
    extensions = {extension(path) for paths in mates.values() for path in paths}
    if len(extensions) > 1:
        findings.append(Finding('error', 'mixed_extensions', 'Read files have different extensions '
                                                             f'({", ".join(sorted(extensions))}); use one per run.'))
    samples = [Sample(name, mates[(name, 1)][0], mates[(name, 2)][0] if paired else None) for name in names
               if (name, 1) in mates and (not paired or (name, 2) in mates)]
    return samples, findings


def link(samples: list, target_dir: str) -> None:
    """
    Symlink the reads into ``target_dir`` as {sample}_1.<ext> and {sample}_2.<ext> (or {sample}.<ext>), the
    names the workflow expects. The original files are not touched.

    :raise FileExistsError: ``target_dir`` exists and is not empty.
    """
    if os.path.isdir(target_dir) and os.listdir(target_dir):
        raise FileExistsError(f'{target_dir} is not empty; it is never overwritten.')
    os.makedirs(target_dir, exist_ok=True)
    for sample in samples:
        reads = [(sample.read1, '_1' if sample.read2 else ''), (sample.read2, '_2')]
        for path, mate in reads[:2 if sample.read2 else 1]:
            os.symlink(os.path.abspath(path), os.path.join(target_dir, f'{sample.name}{mate}.{extension(path)}'))


def write_sheet(samples: list, path: str) -> None:
    """Write the sample sheet (tab-separated: sample, read1, read2), the record of which files a run used."""
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'read1', 'read2'])
        writer.writerows([sample.name, sample.read1, sample.read2 or ''] for sample in samples)


def read_sheet(path: str) -> list:
    """The samples of a sheet written by :func:`write_sheet`."""
    with open(path, newline='') as handle:
        return [Sample(row['sample'], row['read1'], row['read2'] or None) for row in csv.DictReader(handle, delimiter='\t')]
