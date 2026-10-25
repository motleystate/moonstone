import os
import tempfile
from unittest import TestCase

import pandas as pd

from moonstone.validation import Finding, check_file_format, has_errors


class TestCheckFileFormat(TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, content, encoding="utf-8"):
        path = os.path.join(self.tmp.name, name)
        with open(path, "wb") as fh:
            fh.write(content if isinstance(content, bytes) else content.encode(encoding))
        return path

    def codes(self, findings):
        return [f.code for f in findings]

    def by_code(self, findings, code):
        return next(f for f in findings if f.code == code)

    def test_clean_english_csv(self):
        path = self.write("ok.csv", "sample,date,count,ratio\nS1,2024-02-13,12,0.5\nS2,2024-03-01,3,1e-3\n")
        self.assertEqual(check_file_format(path), [])

    def test_clean_tsv(self):
        path = self.write("ok.tsv", "taxon\tS1\tS2\nE. coli\t1\t2\n")
        self.assertEqual(check_file_format(path), [])
        self.assertEqual(check_file_format(path, sep="\t"), [])

    def test_french_excel_csv(self):
        content = "Échantillon;Date;Abondance\nS1;13/02/2024;1,5\nS2;01/03/2024;1 234,75\n"
        path = self.write("fr.csv", content, encoding="cp1252")
        findings = check_file_format(path)
        self.assertTrue(has_errors(findings))
        self.assertEqual(
            self.codes(findings),
            ["not_utf8", "semicolon_separator", "day_first_dates", "decimal_comma", "french_thousands"])
        not_utf8 = self.by_code(findings, "not_utf8")
        self.assertEqual(not_utf8.row, 1)
        self.assertIn("0xC9", not_utf8.message)
        self.assertIn("'É'", not_utf8.message)
        decimal = self.by_code(findings, "decimal_comma")
        self.assertEqual((decimal.column, decimal.row, decimal.example), ("Abondance", 2, "1,5"))
        dates = self.by_code(findings, "day_first_dates")
        self.assertEqual((dates.column, dates.row, dates.example), ("Date", 2, "13/02/2024"))

    def test_file_is_not_modified(self):
        content = "a;b\n1,5;2\n".encode("cp1252")
        path = self.write("fr.csv", content)
        check_file_format(path)
        with open(path, "rb") as fh:
            self.assertEqual(fh.read(), content)

    def test_not_utf8_reports_line(self):
        path = self.write("late.csv", b"name,count\nok,1\ncaf\xe9,2\n")
        finding = self.by_code(check_file_format(path), "not_utf8")
        self.assertEqual(finding.row, 3)
        self.assertIn("'é'", finding.message)

    def test_utf8_bom_is_info(self):
        path = self.write("bom.csv", "﻿sample,count\nS1,3\n")
        findings = check_file_format(path)
        self.assertEqual(self.codes(findings), ["utf8_bom"])
        self.assertEqual(findings[0].level, "info")
        self.assertFalse(has_errors(findings))

    def test_utf16_stops(self):
        path = self.write("u16.txt", "sample\tcount\nS1\t3\n".encode("utf-16"))
        findings = check_file_format(path)
        self.assertEqual(self.codes(findings), ["utf16"])
        self.assertTrue(has_errors(findings))

    def test_separator_mismatch(self):
        path = self.write("semi.csv", "a;b;c\n1;2;3\n")
        finding = self.by_code(check_file_format(path, sep=","), "separator_mismatch")
        self.assertEqual(finding.level, "error")

    def test_semicolon_given_explicitly_is_not_reported(self):
        path = self.write("semi.csv", "a;b;c\n1;2;3\n")
        self.assertEqual(check_file_format(path, sep=";"), [])

    def test_single_column_separator_unknown(self):
        path = self.write("one.csv", "a\n1\n")
        self.assertEqual(self.codes(check_file_format(path)), ["separator_unknown"])

    def test_empty_file(self):
        path = self.write("empty.csv", "")
        self.assertEqual(self.codes(check_file_format(path)), ["empty_file"])

    def test_quoted_decimal_comma_in_comma_file(self):
        path = self.write("q.csv", 's,v\nS1,"1,5"\nS2,2\n')
        finding = self.by_code(check_file_format(path), "decimal_comma")
        self.assertEqual((finding.column, finding.row, finding.example), ("v", 2, "1,5"))

    def test_ambiguous_thousands_comma(self):
        path = self.write("amb.csv", 's,v\nS1,"1,234"\nS2,2\n')
        finding = self.by_code(check_file_format(path), "ambiguous_comma")
        self.assertEqual(finding.level, "warning")
        self.assertIn("could be 1234", finding.message)

    def test_ambiguous_comma_next_to_french_values(self):
        path = self.write("amb.csv", 's;v\nS1;1,234\nS2;2,5\n')
        finding = self.by_code(check_file_format(path, sep=";"), "ambiguous_comma")
        self.assertIn("probably decimal commas", finding.message)

    def test_english_thousands(self):
        path = self.write("en.csv", 's,v\nS1,"1,234,567"\nS2,"1,234.5"\n')
        finding = self.by_code(check_file_format(path), "english_thousands")
        self.assertEqual(finding.level, "warning")
        self.assertIn("2 value(s)", finding.message)

    def test_french_thousands_variants(self):
        for value in ["1 234", "1 234,5", "12 345", "1.234.567", "1.234,5"]:
            with self.subTest(value=value):
                path = self.write("t.csv", f"s;v\nS1;{value}\n")
                self.assertIn("french_thousands", self.codes(check_file_format(path, sep=";")))

    def test_english_numbers_pass(self):
        path = self.write("n.csv", "v\n1.234\n-0.5\n.25\n3E+2\n+7\n")
        self.assertEqual(check_file_format(path, sep=","), [])

    def test_date_orders(self):
        cases = {
            "month_first_dates": ["02/13/2024", "03/01/2024"],
            "ambiguous_dates": ["01/02/2024", "03/04/2024"],
            "mixed_date_order": ["13/02/2024", "02/13/2024"],
            "day_first_dates": ["31.12.2023 10:30", "01.01.2024"],
        }
        for code, dates in cases.items():
            with self.subTest(code=code):
                path = self.write("d.csv", "s,d\n" + "".join(f"S{i},{d}\n" for i, d in enumerate(dates)))
                findings = check_file_format(path)
                self.assertEqual(self.codes(findings), [code])
        path = self.write("d.csv", "s,d\nS1,13/02/2024\nS2,02/13/2024\n")
        self.assertTrue(has_errors(check_file_format(path)))

    def test_excel_text_cells_checked(self):
        path = os.path.join(self.tmp.name, "fr.xlsx")
        pd.DataFrame({
            "sample": ["S1", "S2"],
            "abundance_number": [1.5, 2.0],
            "abundance_text": ["1,5", "2"],
            "date_value": pd.to_datetime(["2024-02-13", "2024-03-01"]),
            "date_text": ["13/02/2024", "01/03/2024"],
        }).to_excel(path, index=False)
        findings = check_file_format(path)
        self.assertEqual(self.codes(findings), ["decimal_comma", "day_first_dates"])
        self.assertEqual(self.by_code(findings, "decimal_comma").column, "abundance_text")
        self.assertEqual(self.by_code(findings, "decimal_comma").row, 2)
        self.assertEqual(self.by_code(findings, "day_first_dates").column, "date_text")

    def test_clean_excel(self):
        path = os.path.join(self.tmp.name, "ok.xlsx")
        pd.DataFrame({"sample": ["S1"], "count": [3]}).to_excel(path, index=False)
        self.assertEqual(check_file_format(path), [])

    def test_finding_str(self):
        finding = Finding("error", "decimal_comma", "Comma decimal.", "v", 2, "1,5")
        self.assertEqual(str(finding), "ERROR [decimal_comma] Comma decimal. (column 'v', row 2, value '1,5')")
        self.assertEqual(str(Finding("info", "utf8_bom", "BOM.")), "INFO [utf8_bom] BOM.")
