"""Opt-in conversion of text input files to the conventions Moonstone expects.

Nothing in Moonstone calls these functions: they run only when the user asks,
after :func:`~moonstone.validation.format_check.check_file_format` has reported
a problem and the user has confirmed what the file's conventions are.

Both functions write a **new** file (they refuse to overwrite one), return every
change they made, and run the format check again on the output. A column is
converted only when every value in it fits the declared convention; otherwise
it is left as it is and reported, so that a half-converted column cannot happen.
"""

import codecs
import csv
import datetime
import os
import re
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import pandas as pd

from moonstone.validation.format_check import EXCEL_EXTENSIONS, Finding, _SLASH_DATE, check_file_format

DEFAULT_NA_VALUES = ("", "NA", "N/A", "NaN", "nan", "null", "NULL")
DATE_ORDERS = ("dmy", "mdy")


@dataclass(frozen=True)
class Change:
    """One value changed by a conversion.

    :param row: 1-based line number in the source file.
    :param column: column name, or ``None`` for a whole-line change (re-encoding).
    :param old: value as it was in the source file.
    :param new: value as written to the output file.
    """

    row: int
    column: Optional[str]
    old: str
    new: str


@dataclass
class ConversionResult:
    """What a conversion did.

    :param output: path of the file written.
    :param changes: every value that changed.
    :param notes: findings about the conversion itself, such as columns left unconverted.
    :param findings: result of :func:`check_file_format` on the output file.
    """

    output: str
    changes: List[Change] = field(default_factory=list)
    notes: List[Finding] = field(default_factory=list)
    findings: List[Finding] = field(default_factory=list)

    def changes_frame(self) -> pd.DataFrame:
        """Changes as a table (row, column, old, new), e.g. to save next to the output file."""
        return pd.DataFrame(
            [(c.row, c.column, c.old, c.new) for c in self.changes], columns=["row", "column", "old", "new"])


def convert_to_utf8(src: str, dst: str, source_encoding: str) -> ConversionResult:
    """Re-encode a text file to UTF-8.

    The source encoding must be given: cp1252, latin-1 and Mac Roman cannot be told
    apart reliably, so it is never guessed. Decoding is strict: a byte that does not
    exist in ``source_encoding`` stops the conversion instead of being replaced.

    :param src: file to convert.
    :param dst: new file to write (must not exist).
    :param source_encoding: encoding of ``src``, e.g. ``"cp1252"`` for Excel files
        saved in a Western European locale, or ``"utf-8-sig"`` to only remove a byte order mark.
    :raise FileExistsError: ``dst`` exists.
    :raise ValueError: ``src`` is not valid in ``source_encoding``, or is already UTF-8
        with non-ASCII characters (decoding it again would garble them).
    """
    _check_paths(src, dst)
    with open(src, "rb") as fh:
        raw = fh.read()
    if _is_non_ascii_utf8(raw) and source_encoding.replace("_", "-").lower() not in ("utf-8", "utf8", "utf-8-sig"):
        raise ValueError(
            f"{src} is already valid UTF-8 with non-ASCII characters; decoding it as {source_encoding} "
            "would garble them. Nothing was written.")
    try:
        text = raw.decode(source_encoding)
    except UnicodeDecodeError as err:
        line = raw.count(b"\n", 0, err.start) + 1
        raise ValueError(
            f"{src} is not valid {source_encoding}: byte 0x{raw[err.start:err.start + 1].hex().upper()} "
            f"on line {line}. Nothing was written.") from err

    result = ConversionResult(dst)
    codec = codecs.lookup(source_encoding).name
    if codec.startswith(("utf-16", "utf-32")):
        result.notes.append(Finding("info", "multibyte_source", f"Every line changed bytes ({codec} to UTF-8)."))
    elif not codec.startswith("utf-8"):
        for number, line in enumerate(text.splitlines(), start=1):
            if not line.isascii():
                old = line.encode(source_encoding).decode("ascii", errors="backslashreplace")
                result.changes.append(Change(number, None, old, line))
    if text.startswith("\ufeff"):
        # A BOM decoded by a non-"-sig" codec stays as U+FEFF; drop it so the output has none.
        text = text[1:]
        result.notes.append(Finding("info", "bom_removed", "The byte order mark was removed."))
    elif raw.startswith(codecs.BOM_UTF8):
        result.notes.append(Finding("info", "bom_removed", "The byte order mark was removed."))
    with open(dst, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
    result.notes.insert(0, Finding(
        "info", "reencoded",
        f"Re-encoded from {source_encoding} to UTF-8; {len(result.changes)} line(s) with non-ASCII characters."))
    result.findings = check_file_format(dst)
    return result


def convert_locale(src: str, dst: str, sep: str = ";", decimal: str = ",", thousands: Optional[str] = None,
                   date_order: Optional[str] = "dmy", out_sep: str = ",",
                   na_values: Sequence[str] = DEFAULT_NA_VALUES) -> ConversionResult:
    """Rewrite a UTF-8 text file with English conventions, from conventions the user declares.

    The output uses ``out_sep`` as separator, a dot as decimal mark, no thousands
    separator, and ISO dates (yyyy-mm-dd, keeping any time). The header row is never changed.

    A column is converted as numbers when every value (other than ``na_values``) is a
    number in the declared convention, and as dates when every value is a date in
    ``date_order``. A column where only some values fit is left unchanged and reported
    in ``notes``. Text columns such as sample names or taxonomy are left unchanged.

    With ``decimal=","`` declared, ``1,234`` means 1.234 and is converted to ``1.234``.

    :param src: UTF-8 file to convert (run :func:`convert_to_utf8` first if it is not).
    :param dst: new file to write (must not exist).
    :param sep: separator used in ``src``.
    :param decimal: decimal mark used in ``src`` (``","`` or ``"."``).
    :param thousands: characters used as thousands separators in ``src``. By default
        space, no-break space, narrow no-break space and the character that is not the decimal mark.
    :param date_order: ``"dmy"`` (dd/mm/yyyy), ``"mdy"`` (mm/dd/yyyy) or ``None`` to leave dates unchanged.
    :param out_sep: separator for the output file.
    :param na_values: values treated as missing; they are kept as written.
    :raise FileExistsError: ``dst`` exists.
    :raise ValueError: invalid arguments, a spreadsheet, or ``src`` not valid UTF-8.
    """
    _check_paths(src, dst)
    thousands = _locale_thousands(src, decimal, thousands, date_order)
    with open(src, "rb") as fh:
        raw = fh.read()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as err:
        raise ValueError(f"{src} is not valid UTF-8; run convert_to_utf8 first. Nothing was written.") from err

    reader = csv.reader(text.splitlines(keepends=True), delimiter=sep)
    header = next(reader, [])
    rows, line_numbers = [], []
    for record in reader:
        rows.append(record)
        line_numbers.append(reader.line_num)

    number = _number_pattern(decimal, thousands)
    na = set(na_values)
    result = ConversionResult(dst)
    width = max([len(header)] + [len(r) for r in rows])
    for i in range(width):
        name = header[i] if i < len(header) and header[i] else f"column {i + 1}"
        _convert_column(result, rows, line_numbers, i, name, number, decimal, thousands, date_order, na)

    with open(dst, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter=out_sep, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    if sep != out_sep:
        result.notes.insert(0, Finding("info", "separator_changed", f"Separator changed from {sep!r} to {out_sep!r}."))
    result.findings = check_file_format(dst, sep=out_sep)
    return result


def _check_paths(src, dst):
    if os.path.abspath(src) == os.path.abspath(dst):
        raise ValueError("The output file must be different from the input file.")
    if os.path.exists(dst):
        raise FileExistsError(f"{dst} already exists; it is never overwritten.")


def _locale_thousands(src, decimal, thousands, date_order):
    """Check convert_locale arguments and return the thousands separators to accept."""
    if str(src).lower().rsplit(".", 1)[-1] in EXCEL_EXTENSIONS:
        raise ValueError("Spreadsheets are not converted; save the sheet as 'CSV UTF-8' first.")
    if decimal not in (",", "."):
        raise ValueError("decimal must be ',' or '.'.")
    if date_order not in DATE_ORDERS + (None,):
        raise ValueError(f"date_order must be one of {DATE_ORDERS} or None.")
    if thousands is None:
        thousands = " \u00a0\u202f" + ("." if decimal == "," else ",")
    if decimal in thousands:
        raise ValueError("The decimal mark cannot also be a thousands separator.")
    return thousands


def _convert_column(result, rows, line_numbers, i, name, number, decimal, thousands, date_order, na):
    """Convert column ``i`` of ``rows`` in place when all its values fit, recording changes and notes."""
    values = [(line, r[i].strip()) for line, r in zip(line_numbers, rows) if i < len(r) and r[i].strip() not in na]
    new_values = _convert_numbers(values, number, decimal, thousands)
    kind = "numbers"
    if new_values is None and date_order is not None:
        new_values, kind = _convert_dates(values, date_order), "dates"
    if isinstance(new_values, Finding):
        result.notes.append(Finding(new_values.level, new_values.code, f"Column not converted: {new_values.message}",
                                    name, new_values.row, new_values.example))
        return
    if not new_values:
        return
    for (line, old), new in zip(values, new_values):
        if old != new:
            result.changes.append(Change(line, name, old, new))
    converted = dict(zip((line for line, _ in values), new_values))
    for line, r in zip(line_numbers, rows):
        if i < len(r) and line in converted:
            r[i] = converted[line]
    result.notes.append(Finding("info", f"converted_{kind}", f"Converted as {kind}.", name))


def _is_non_ascii_utf8(raw):
    if raw.isascii():
        return False
    try:
        raw.decode("utf-8")
    except UnicodeDecodeError:
        return False
    return True


def _number_pattern(decimal, thousands):
    group = "[" + re.escape(thousands) + "]"
    return re.compile(
        rf"([+-]?)(\d{{1,3}}(?:{group}\d{{3}})+|\d+)(?:{re.escape(decimal)}(\d+))?([eE][+-]?\d+)?")


def _convert_numbers(values, pattern, decimal, thousands):
    """New values, None when the column has no numbers to convert, or a Finding when only some values fit."""
    matches = [pattern.fullmatch(v) for _, v in values]
    marks = decimal + thousands
    if not any(m and any(c in v for c in marks) for m, (_, v) in zip(matches, values)):
        return None
    for m, (line, v) in zip(matches, values):
        if m is None:
            return Finding("warning", "partly_numeric",
                           f"some values are numbers with {decimal!r} as decimal mark but this one is not.",
                           row=line, example=v)
    out = []
    for m in matches:
        sign, integer, fraction, exponent = m.groups()
        integer = re.sub("[" + re.escape(thousands) + "]", "", integer)
        out.append(sign + integer + ("." + fraction if fraction is not None else "") + (exponent or ""))
    return out


def _convert_dates(values, date_order):
    """New ISO values, None when the column has no dates, or a Finding when a value cannot be converted."""
    matches = [_SLASH_DATE.fullmatch(v) for _, v in values]
    if not any(matches):
        return None
    out = []
    for m, (line, v) in zip(matches, values):
        if m is None:
            return Finding("warning", "partly_dates", "some values are dates but this one is not.", row=line, example=v)
        first, second, year, time = int(m.group(1)), int(m.group(2)), m.group(3), m.group(4) or ""
        if len(year) == 2:
            return Finding("warning", "two_digit_year", "two-digit years are ambiguous (which century?).",
                           row=line, example=v)
        day, month = (first, second) if date_order == "dmy" else (second, first)
        try:
            date = datetime.date(int(year), month, day)
        except ValueError:
            return Finding("warning", "invalid_date", f"not a valid date in {date_order} order.", row=line, example=v)
        out.append(date.isoformat() + (" " + time.strip(" T") if time else ""))
    return out
