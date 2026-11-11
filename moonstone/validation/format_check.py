"""First-pass format check for CSV, TSV and Excel input files.

Moonstone expects English-style input: UTF-8 text, comma (or tab) separators,
a dot as the decimal mark, and ISO dates (yyyy-mm-dd). Files saved by Excel in
a French locale often differ on every point: cp1252 encoding, ";" separators,
decimal commas, spaces as thousands separators and dd/mm/yyyy dates.

:func:`check_file_format` reads the file and reports what it finds. It never
changes, converts or guesses data: callers decide what to do with the findings.
"""

import codecs
import csv
import re
from dataclasses import dataclass
from typing import List, Optional

import pandas as pd

EXCEL_EXTENSIONS = ("xls", "xlsx", "xlsm", "xlsb", "ods", "odf")
CANDIDATE_SEPARATORS = (",", ";", "\t")
SEPARATOR_NAMES = {",": "comma", ";": "semicolon", "\t": "tab"}
SNIFF_LINES = 50

# English number: 12, -3.5, .5, 1e-3
_EN_NUMBER = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
# French decimal comma: 1,5 or 0,25 (not exactly three digits after the comma)
_DECIMAL_COMMA = re.compile(r"[+-]?\d+,(\d{1,2}|\d{4,})")
# "1,234": English thousands separator (1234) or French decimal comma (1.234)
_AMBIGUOUS_COMMA = re.compile(r"[+-]?\d{1,3},\d{3}")
# French grouping: 1 234 / 1 234,5 / 1.234.567 / 1.234,5 (space, no-break space, narrow no-break space or dot)
_FRENCH_GROUPED = re.compile(r"[+-]?\d{1,3}(( \d{3})+|( \d{3})+|( \d{3})+|(\.\d{3})+)(,\d+)?")
# English grouping that pandas reads as text: 1,234,567 or 1,234.5
_ENGLISH_GROUPED = re.compile(r"[+-]?\d{1,3}((,\d{3}){2,}(\.\d+)?|,\d{3}\.\d+)")
# dd/mm/yyyy or mm/dd/yyyy, with "/", "." or "-", optionally followed by a time
_SLASH_DATE = re.compile(r"(\d{1,2})[/.-](\d{1,2})[/.-](\d{2}|\d{4})([ T]\d{1,2}:\d{2}(:\d{2})?)?")


@dataclass(frozen=True)
class Finding:
    """One problem or observation about the file.

    :param level: ``"error"`` (the file will load wrongly or not at all),
        ``"warning"`` (likely problem, needs a human look) or ``"info"``.
    :param code: short stable identifier, e.g. ``"not_utf8"`` or ``"decimal_comma"``.
    :param message: explanation for the user.
    :param column: column name, when the finding is about one column.
    :param row: 1-based line number (CSV) or row number (Excel) of the first example.
    :param example: the first offending value, as written in the file.
    """

    level: str
    code: str
    message: str
    column: Optional[str] = None
    row: Optional[int] = None
    example: Optional[str] = None

    def __str__(self):
        where = []
        if self.column is not None:
            where.append(f"column {self.column!r}")
        if self.row is not None:
            where.append(f"row {self.row}")
        if self.example is not None:
            where.append(f"value {self.example!r}")
        suffix = f" ({', '.join(where)})" if where else ""
        return f"{self.level.upper()} [{self.code}] {self.message}{suffix}"


def check_file_format(path: str, sep: Optional[str] = None, sheet_name=0) -> List[Finding]:
    """Check the encoding, separator, number format and date format of an input file.

    :param path: CSV/TSV/text file, or Excel/OpenDocument spreadsheet.
    :param sep: separator the caller will use to read a text file. When ``None``
        the separator is detected and a ";" separator is reported.
    :param sheet_name: sheet to check in a spreadsheet (same as ``pandas.read_excel``).
    :return: list of findings, empty when nothing was found.
    """
    if str(path).lower().rsplit(".", 1)[-1] in EXCEL_EXTENSIONS:
        return _check_excel(path, sheet_name)
    return _check_text(path, sep)


def has_errors(findings: List[Finding]) -> bool:
    """True when at least one finding has level ``"error"``."""
    return any(f.level == "error" for f in findings)


def _check_text(path, sep):
    with open(path, "rb") as fh:
        raw = fh.read()
    findings, text = _check_encoding(raw)
    if text is None:
        return findings
    lines = text.splitlines()
    if not lines:
        findings.append(Finding("error", "empty_file", "The file is empty."))
        return findings

    detected = _detect_separator(lines[:SNIFF_LINES])
    if detected is None:
        if sep is None:
            findings.append(Finding(
                "warning", "separator_unknown",
                "Could not detect a comma, semicolon or tab separator; the file may have a single column."))
    elif sep is not None and detected != sep:
        findings.append(Finding(
            "error", "separator_mismatch",
            f"The file looks {SEPARATOR_NAMES[detected]}-separated but will be read with {sep!r}."))
    elif sep is None and detected == ";":
        findings.append(Finding(
            "warning", "semicolon_separator",
            "Semicolon separator, typical of CSV files saved by Excel in a French locale; "
            "Moonstone expects commas (or tabs)."))
    use_sep = sep if sep is not None else (detected or ",")

    reader = csv.reader(text.splitlines(keepends=True), delimiter=use_sep)
    header = next(reader, [])
    columns = _ColumnScanner(header)
    for record in reader:
        for i, cell in enumerate(record):
            columns.add(i, cell, reader.line_num)
    return findings + columns.findings()


def _check_excel(path, sheet_name):
    df = pd.read_excel(path, header=None, dtype=object, sheet_name=sheet_name)
    if df.empty:
        return [Finding("error", "empty_file", "The sheet is empty.")]
    columns = _ColumnScanner(["" if pd.isna(v) else str(v) for v in df.iloc[0]])
    for row_number, values in enumerate(df.iloc[1:].itertuples(index=False), start=2):
        for i, value in enumerate(values):
            # Real numbers and dates are stored as values and are fine; only text can be locale-formatted.
            if isinstance(value, str):
                columns.add(i, value, row_number)
    return columns.findings()


def _check_encoding(raw):
    """Return (findings, decoded text or None when the text cannot be inspected)."""
    if raw.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return [Finding(
            "error", "utf16",
            "The file is UTF-16 (Excel 'Unicode Text'); Moonstone expects UTF-8. "
            "Save it as 'CSV UTF-8' instead.")], None
    findings = []
    if raw.startswith(codecs.BOM_UTF8):
        findings.append(Finding(
            "info", "utf8_bom",
            "The file starts with a UTF-8 byte order mark (Excel 'CSV UTF-8'). pandas ignores it, "
            "but open(..., encoding='utf-8') keeps it in the first column name; use 'utf-8-sig'."))
        raw = raw[len(codecs.BOM_UTF8):]
    try:
        return findings, raw.decode("utf-8")
    except UnicodeDecodeError as err:
        bad = raw[err.start:err.start + 1]
        line = raw.count(b"\n", 0, err.start) + 1
        as_cp1252 = bad.decode("cp1252", errors="replace")
        findings.append(Finding(
            "error", "not_utf8",
            f"The file is not valid UTF-8: byte 0x{bad.hex().upper()} would be {as_cp1252!r} in cp1252 "
            "(Windows/Excel Western European). Re-save the file as 'CSV UTF-8'. "
            "The other checks below read the file as cp1252 for inspection only.",
            row=line))
        return findings, raw.decode("cp1252", errors="replace")


def _detect_separator(lines):
    """Pick the candidate giving the same number (>1) of fields on every sampled line."""
    best, best_fields = None, 1
    sample = [line for line in lines if line.strip()]
    for candidate in CANDIDATE_SEPARATORS:
        counts = {len(row) for row in csv.reader(sample, delimiter=candidate)}
        if len(counts) == 1:
            fields = counts.pop()
            if fields > best_fields:
                best, best_fields = candidate, fields
    return best


class _ColumnScanner:
    """Collect locale-format evidence per column, then turn it into findings."""

    KINDS = ("decimal_comma", "french_grouped", "ambiguous_comma", "english_grouped")

    def __init__(self, header):
        self.header = header
        self.first = {}       # (column index, kind) -> (row, value) of the first example
        self.counts = {}      # (column index, kind) -> number of cells
        self.day_first = {}   # column index -> (row, value) of a date whose first part is > 12
        self.month_first = {}  # column index -> (row, value) of a date whose second part is > 12
        self.dates = {}       # column index -> (row, value) of the first non-ISO date

    def _name(self, i):
        return self.header[i] if i < len(self.header) and self.header[i] else f"column {i + 1}"

    def _note(self, i, kind, row, value):
        self.first.setdefault((i, kind), (row, value))
        self.counts[(i, kind)] = self.counts.get((i, kind), 0) + 1

    def add(self, i, cell, row):
        value = cell.strip()
        if not value or value.isdigit() or _EN_NUMBER.fullmatch(value):
            return
        if _DECIMAL_COMMA.fullmatch(value):
            self._note(i, "decimal_comma", row, value)
        elif _AMBIGUOUS_COMMA.fullmatch(value):
            self._note(i, "ambiguous_comma", row, value)
        elif _FRENCH_GROUPED.fullmatch(value):
            self._note(i, "french_grouped", row, value)
        elif _ENGLISH_GROUPED.fullmatch(value):
            self._note(i, "english_grouped", row, value)
        else:
            match = _SLASH_DATE.fullmatch(value)
            if match:
                self.dates.setdefault(i, (row, value))
                first, second = int(match.group(1)), int(match.group(2))
                if first > 12:
                    self.day_first.setdefault(i, (row, value))
                if second > 12:
                    self.month_first.setdefault(i, (row, value))

    def findings(self):
        out = []
        columns = sorted({i for i, _ in self.first} | set(self.dates))
        for i in columns:
            name = self._name(i)
            french = sum(self.counts.get((i, k), 0) for k in ("decimal_comma", "french_grouped"))
            if (i, "decimal_comma") in self.first:
                row, value = self.first[(i, "decimal_comma")]
                out.append(Finding(
                    "error", "decimal_comma",
                    f"{self.counts[(i, 'decimal_comma')]} value(s) use a comma as the decimal mark; "
                    "Moonstone expects a dot, and pandas reads these as text.",
                    name, row, value))
            if (i, "french_grouped") in self.first:
                row, value = self.first[(i, "french_grouped")]
                out.append(Finding(
                    "error", "french_thousands",
                    f"{self.counts[(i, 'french_grouped')]} value(s) use a space or dot as thousands separator; "
                    "pandas reads these as text.",
                    name, row, value))
            if (i, "ambiguous_comma") in self.first:
                row, value = self.first[(i, "ambiguous_comma")]
                meaning = ("probably decimal commas, as other values in this column are French-formatted"
                           if french else "could be 1234 (English thousands) or 1.234 (French decimal)")
                out.append(Finding(
                    "error", "ambiguous_comma",
                    f"{self.counts[(i, 'ambiguous_comma')]} value(s) like '1,234' {meaning}; "
                    "numbers must not contain commas, and pandas reads these as text.",
                    name, row, value))
            if (i, "english_grouped") in self.first:
                row, value = self.first[(i, "english_grouped")]
                out.append(Finding(
                    "error", "english_thousands",
                    f"{self.counts[(i, 'english_grouped')]} value(s) use a comma as thousands separator; "
                    "numbers must not contain commas, and pandas reads these as text.",
                    name, row, value))
            if i in self.dates:
                out.append(self._date_finding(i, name))
        return out

    def _date_finding(self, i, name):
        advice = "Use ISO dates (yyyy-mm-dd)."
        if i in self.day_first and i in self.month_first:
            row, value = self.month_first[i]
            return Finding(
                "error", "mixed_date_order",
                f"Dates in both dd/mm and mm/dd order (e.g. {self.day_first[i][1]!r} and {value!r}). {advice}",
                name, row, value)
        if i in self.day_first:
            row, value = self.day_first[i]
            return Finding("warning", "day_first_dates", f"Dates are day first (dd/mm/yyyy). {advice}",
                           name, row, value)
        if i in self.month_first:
            row, value = self.month_first[i]
            return Finding("warning", "month_first_dates", f"Dates are month first (mm/dd/yyyy). {advice}",
                           name, row, value)
        row, value = self.dates[i]
        return Finding(
            "warning", "ambiguous_dates",
            f"Dates could be dd/mm or mm/dd: no day above 12 tells them apart. {advice}",
            name, row, value)
