"""Match the samples of a count table with the samples of a metadata file by sample ID.

A Moonstone project starts with two files: normalized counts and metadata. They describe the same
samples, but not necessarily all of them and not in the same order, so samples are paired by ID,
never by position or by number of rows.

:func:`match_samples` compares the IDs as text, exactly as written. It never renames, trims or
reorders anything: it reports what it finds and the caller decides what to do with it.
"""

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, List, Tuple

import pandas as pd

from moonstone.validation.format_check import Finding

MAX_LISTED = 10


def _listing(ids: List[str]) -> str:
    shown = ", ".join(ids[:MAX_LISTED])
    return f"{shown} ... ({len(ids)} in all)" if len(ids) > MAX_LISTED else shown


@dataclass(frozen=True)
class SampleMatch:
    """Result of :func:`match_samples`. All lists keep the order of the file they come from.

    :param matched: IDs found in both files, in count-table order.
    :param only_in_counts: IDs of the count table with no metadata row.
    :param only_in_metadata: IDs of the metadata with no column in the count table.
    :param duplicated_in_counts: IDs that appear more than once in the count table.
    :param duplicated_in_metadata: IDs that appear more than once in the metadata.
    :param blank_in_metadata: number of metadata rows with no sample ID.
    :param same_order: whether the matched samples come in the same order in both files.
    :param near_matches: (count ID, metadata ID) pairs left unmatched that differ only by case or spaces.
    """

    matched: List[str]
    only_in_counts: List[str]
    only_in_metadata: List[str]
    duplicated_in_counts: List[str]
    duplicated_in_metadata: List[str]
    blank_in_metadata: int
    same_order: bool
    near_matches: List[Tuple[str, str]]

    def findings(self) -> List[Finding]:
        """Explain the result: errors stop the analysis, warnings and info do not."""
        findings = []
        if self.duplicated_in_counts:
            findings.append(Finding("error", "duplicated_sample_in_counts",
                                    f"Sample IDs used more than once in the count table: "
                                    f"{_listing(self.duplicated_in_counts)}"))
        if self.duplicated_in_metadata:
            findings.append(Finding("error", "duplicated_sample_in_metadata",
                                    f"Sample IDs used more than once in the metadata: "
                                    f"{_listing(self.duplicated_in_metadata)}"))
        if not self.matched:
            findings.append(Finding("error", "no_common_samples",
                                    "No sample ID is found in both the count table and the metadata"))
        if self.blank_in_metadata:
            findings.append(Finding("warning", "blank_sample_id",
                                    f"{self.blank_in_metadata} metadata row(s) have no sample ID and are ignored"))
        if self.only_in_counts:
            findings.append(Finding("warning", "sample_without_metadata",
                                    f"{len(self.only_in_counts)} sample(s) of the count table have no metadata and "
                                    f"are left out: {_listing(self.only_in_counts)}"))
        if self.near_matches:
            pairs = [f"{c!r} / {m!r}" for c, m in self.near_matches]
            findings.append(Finding("warning", "near_match",
                                    f"Sample IDs that differ only by case or spaces (not matched): {_listing(pairs)}"))
        if self.only_in_metadata:
            findings.append(Finding("info", "sample_without_counts",
                                    f"{len(self.only_in_metadata)} sample(s) of the metadata are not in the count "
                                    f"table and are ignored: {_listing(self.only_in_metadata)}"))
        if self.matched and not self.same_order:
            findings.append(Finding("info", "different_order",
                                    "Samples are not in the same order in the two files; they are matched by ID"))
        return findings


def match_samples(count_samples: Iterable, metadata_samples: Iterable) -> SampleMatch:
    """Pair the sample IDs of a count table and a metadata file.

    :param count_samples: sample IDs of the count table (its columns or its index, depending on orientation).
    :param metadata_samples: sample IDs of the metadata (usually its index).
    """
    counts = [str(s) for s in count_samples]
    metadata_raw = list(metadata_samples)
    metadata = [str(s) for s in metadata_raw if not pd.isna(s)]
    blank = len(metadata_raw) - len(metadata)

    in_counts, in_metadata = set(counts), set(metadata)
    matched = list(dict.fromkeys(s for s in counts if s in in_metadata))
    only_in_counts = list(dict.fromkeys(s for s in counts if s not in in_metadata))
    only_in_metadata = list(dict.fromkeys(s for s in metadata if s not in in_counts))
    metadata_order = list(dict.fromkeys(s for s in metadata if s in in_counts))

    def key(s):
        return "".join(s.split()).casefold()

    unmatched_metadata = {}
    for s in only_in_metadata:
        unmatched_metadata.setdefault(key(s), s)
    near = [(s, unmatched_metadata[key(s)]) for s in only_in_counts if key(s) in unmatched_metadata]

    return SampleMatch(
        matched=matched,
        only_in_counts=only_in_counts,
        only_in_metadata=only_in_metadata,
        duplicated_in_counts=[s for s, n in Counter(counts).items() if n > 1],
        duplicated_in_metadata=[s for s, n in Counter(metadata).items() if n > 1],
        blank_in_metadata=blank,
        same_order=matched == metadata_order,
        near_matches=near,
    )
