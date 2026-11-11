import os
import tempfile
from unittest import TestCase

import pandas as pd

from moonstone.validation import check_file_format, convert_locale, convert_to_utf8, has_errors


class ConvertTestCase(TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def write(self, name, content, encoding="utf-8"):
        with open(self.path(name), "wb") as fh:
            fh.write(content if isinstance(content, bytes) else content.encode(encoding))
        return self.path(name)

    def read(self, name):
        with open(self.path(name), "rb") as fh:
            return fh.read().decode("utf-8")

    def codes(self, findings):
        return [f.code for f in findings]


class TestConvertToUtf8(ConvertTestCase):

    def test_cp1252_to_utf8(self):
        content = "Échantillon;Abondance\nS1;1,5\ncafé;2\n"
        src = self.write("fr.csv", content, encoding="cp1252")
        result = convert_to_utf8(src, self.path("out.csv"), "cp1252")
        self.assertEqual(self.read("out.csv"), content)
        self.assertEqual([(c.row, c.new) for c in result.changes],
                         [(1, "Échantillon;Abondance"), (3, "café;2")])
        self.assertEqual(result.changes[0].old, "\\xc9chantillon;Abondance")
        self.assertNotIn("not_utf8", self.codes(result.findings))
        self.assertIn("decimal_comma", self.codes(result.findings))
        with open(src, "rb") as fh:
            self.assertEqual(fh.read(), content.encode("cp1252"))

    def test_strict_decoding(self):
        # 0x81 does not exist in cp1252
        src = self.write("bad.csv", b"a,b\nx\x81,1\n")
        with self.assertRaisesRegex(ValueError, "0x81 on line 2"):
            convert_to_utf8(src, self.path("out.csv"), "cp1252")
        self.assertFalse(os.path.exists(self.path("out.csv")))

    def test_refuses_already_utf8(self):
        src = self.write("u.csv", "name\ncafé\n")
        with self.assertRaisesRegex(ValueError, "already valid UTF-8"):
            convert_to_utf8(src, self.path("out.csv"), "cp1252")
        self.assertFalse(os.path.exists(self.path("out.csv")))

    def test_never_overwrites(self):
        src = self.write("a.csv", "a\n1\n")
        dst = self.write("b.csv", "keep")
        with self.assertRaises(FileExistsError):
            convert_to_utf8(src, dst, "cp1252")
        self.assertEqual(self.read("b.csv"), "keep")
        with self.assertRaisesRegex(ValueError, "different"):
            convert_to_utf8(src, src, "cp1252")

    def test_remove_bom(self):
        src = self.write("bom.csv", "﻿sample,count\nS1,3\n")
        result = convert_to_utf8(src, self.path("out.csv"), "utf-8-sig")
        self.assertEqual(self.read("out.csv"), "sample,count\nS1,3\n")
        self.assertIn("bom_removed", self.codes(result.notes))
        self.assertEqual(result.findings, [])

    def test_utf16(self):
        src = self.write("u16.txt", "sample\tcount\nÉ1\t3\n".encode("utf-16"))
        result = convert_to_utf8(src, self.path("out.txt"), "utf-16")
        self.assertEqual(self.read("out.txt"), "sample\tcount\nÉ1\t3\n")
        self.assertIn("multibyte_source", self.codes(result.notes))
        self.assertEqual(result.findings, [])


class TestConvertLocale(ConvertTestCase):

    def test_french_file(self):
        src = self.write("fr.csv", (
            "Échantillon;Date;Abondance;Lecture;Note\n"
            "S1;13/02/2024;1,5;1 234;a, b\n"
            "S2;01/03/2024 10:30;1 234,75;12;NA\n"
            "S3;NA;1,234;7;c\n"
        ))
        result = convert_locale(src, self.path("out.csv"))
        self.assertEqual(self.read("out.csv"), (
            "Échantillon,Date,Abondance,Lecture,Note\n"
            "S1,2024-02-13,1.5,1234,\"a, b\"\n"
            "S2,2024-03-01 10:30,1234.75,12,NA\n"
            "S3,NA,1.234,7,c\n"
        ))
        self.assertEqual(result.findings, [])
        self.assertIn("separator_changed", self.codes(result.notes))
        self.assertEqual(
            [(c.row, c.column, c.old, c.new) for c in result.changes],
            [(2, "Date", "13/02/2024", "2024-02-13"), (3, "Date", "01/03/2024 10:30", "2024-03-01 10:30"),
             (2, "Abondance", "1,5", "1.5"), (3, "Abondance", "1 234,75", "1234.75"),
             (4, "Abondance", "1,234", "1.234"), (2, "Lecture", "1 234", "1234")])
        frame = result.changes_frame()
        self.assertEqual(list(frame.columns), ["row", "column", "old", "new"])
        self.assertEqual(len(frame), 6)
        df = pd.read_csv(self.path("out.csv"))
        self.assertEqual(df["Abondance"].tolist(), [1.5, 1234.75, 1.234])

    def test_month_first(self):
        src = self.write("d.csv", "s,d\nS1,02/13/2024\n")
        result = convert_locale(src, self.path("out.csv"), sep=",", date_order="mdy")
        self.assertEqual(self.read("out.csv"), "s,d\nS1,2024-02-13\n")
        self.assertEqual(result.findings, [])

    def test_dates_left_alone_without_order(self):
        src = self.write("d.csv", "s;d\nS1;13/02/2024\n")
        convert_locale(src, self.path("out.csv"), date_order=None)
        self.assertEqual(self.read("out.csv"), "s,d\nS1,13/02/2024\n")

    def test_partly_numeric_column_is_not_converted(self):
        src = self.write("m.csv", "s;v\nS1;1,5\nS2;1.5\n")
        result = convert_locale(src, self.path("out.csv"))
        self.assertEqual(self.read("out.csv"), "s,v\nS1,\"1,5\"\nS2,1.5\n")
        note = next(n for n in result.notes if n.code == "partly_numeric")
        self.assertEqual((note.column, note.row, note.example), ("v", 3, "1.5"))
        self.assertEqual(result.changes, [])
        self.assertTrue(has_errors(result.findings))

    def test_invalid_and_partial_dates_are_not_converted(self):
        cases = {
            "invalid_date": "S1;31/02/2024\n",
            "two_digit_year": "S1;13/02/24\n",
            "partly_dates": "S1;13/02/2024\nS2;soon\n",
        }
        for code, rows in cases.items():
            with self.subTest(code=code):
                src = self.write("d.csv", "s;d\n" + rows)
                dst = self.path(f"{code}.csv")
                result = convert_locale(src, dst)
                self.assertIn(code, self.codes(result.notes))
                self.assertEqual(result.changes, [])

    def test_text_columns_untouched(self):
        src = self.write("t.csv", "taxon;S1\nk__Bacteria;p__Firmicutes;3\n".replace(";p__", "|p__"))
        result = convert_locale(src, self.path("out.csv"))
        self.assertEqual(self.read("out.csv"), "taxon,S1\nk__Bacteria|p__Firmicutes,3\n")
        self.assertEqual(result.changes, [])

    def test_english_thousands_with_dot_decimal(self):
        src = self.write("e.csv", 's,v\nS1,"1,234,567"\nS2,"1,234.5"\n')
        result = convert_locale(src, self.path("out.csv"), sep=",", decimal=".")
        self.assertEqual(self.read("out.csv"), "s,v\nS1,1234567\nS2,1234.5\n")
        self.assertEqual(result.findings, [])

    def test_bom_input(self):
        src = self.write("b.csv", "﻿s;v\nS1;1,5\n")
        convert_locale(src, self.path("out.csv"))
        self.assertEqual(self.read("out.csv"), "s,v\nS1,1.5\n")

    def test_refusals(self):
        src = self.write("a.csv", "a;b\n1;2\n")
        dst = self.write("exists.csv", "keep")
        with self.assertRaises(FileExistsError):
            convert_locale(src, dst)
        with self.assertRaisesRegex(ValueError, "decimal must be"):
            convert_locale(src, self.path("o.csv"), decimal=";")
        with self.assertRaisesRegex(ValueError, "date_order"):
            convert_locale(src, self.path("o.csv"), date_order="ymd")
        with self.assertRaisesRegex(ValueError, "cannot also be"):
            convert_locale(src, self.path("o.csv"), thousands=",")
        with self.assertRaisesRegex(ValueError, "CSV UTF-8"):
            convert_locale(self.path("x.xlsx"), self.path("o.csv"))
        cp = self.write("cp.csv", "é;1\n", encoding="cp1252")
        with self.assertRaisesRegex(ValueError, "convert_to_utf8"):
            convert_locale(cp, self.path("o.csv"))
        self.assertFalse(os.path.exists(self.path("o.csv")))

    def test_full_pipeline_passes_check(self):
        content = "Échantillon;Date;Abondance\nS1;13/02/2024;1,5\nS2;01/03/2024;1 234,75\n"
        src = self.write("fr.csv", content, encoding="cp1252")
        self.assertTrue(has_errors(check_file_format(src)))
        convert_to_utf8(src, self.path("utf8.csv"), "cp1252")
        result = convert_locale(self.path("utf8.csv"), self.path("en.csv"))
        self.assertEqual(result.findings, [])
        self.assertEqual(check_file_format(self.path("en.csv")), [])
