#!/usr/bin/env python3
"""Write the run report (Markdown) of the Moonstone reads workflow from Kraken2/Bracken reports.

Kraken2 report columns: percent, clade reads, direct reads, rank, taxid, name.
"""
import argparse
import datetime
import json
import os
import subprocess


def read_report(path):
    """Return rows (clade_reads, rank, taxid) of a Kraken2-style report."""
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                rows.append((int(p[1]), p[3].strip(), p[4].strip()))
    return rows


def totals(rows):
    """(total reads, classified reads) of a report: unclassified is taxid 0, root is taxid 1."""
    unclassified = sum(r[0] for r in rows if r[2] == "0")
    classified = sum(r[0] for r in rows if r[2] == "1")
    return unclassified + classified, classified


def clade_reads(rows, taxid):
    return sum(r[0] for r in rows if r[2] == str(taxid))


def n_species(path):
    return sum(1 for reads, rank, _ in read_report(path) if rank == "S" and reads > 0)


def tool_version(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return ((out.stdout + out.stderr).strip().splitlines() or ["unknown"])[0]
    except Exception:
        return "not found"


def bracken_version():
    # "bracken -v" needs a dummy argument, and versions before 3.x have no -v flag at all
    v = tool_version(["bracken", "-v", "x"])
    return v if v.startswith("Bracken v") else "unknown (no -v flag in this version)"


def moonstone_version():
    try:
        from moonstone import __version__
        return __version__
    except Exception:
        return "unknown"


def table_columns(tsv):
    """Sample columns of a biom-converted TSV (read from the '#OTU ID' header line only)."""
    with open(tsv) as f:
        for line in f:
            if line.startswith("#OTU ID"):
                return line.rstrip("\n").split("\t")[1:]
            if not line.startswith("#"):
                break
    return []


def in_table(sample, columns):
    """kraken-biom names a column after the report file, e.g. 'S1.taxa' for S1.taxa.report."""
    return any(c == sample or c.startswith(sample + ".") for c in columns)


def build_report(a):
    settings = json.loads(a.settings)
    lines = ["# Moonstone reads pipeline report", ""]
    lines += ["## Run info", "",
              "- Date: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
              "- Moonstone: " + moonstone_version(),
              "- snakemake: " + tool_version(["snakemake", "--version"]),
              "- kraken2: " + tool_version(["kraken2", "--version"]),
              "- bracken: " + bracken_version(),
              "- trimmomatic: " + tool_version(["trimmomatic", "-version"]),
              "", "Settings:", ""]
    lines += ["- `{}`: {}".format(k, v) for k, v in settings.items()]
    host_taxids = settings.get("host_taxid", [9606])
    host_taxids = host_taxids if isinstance(host_taxids, list) else [host_taxids]

    header = ["Sample", "Reads before host removal", "% host", "Reads classified against DB",
              "% classified", "% unclassified"]
    if a.bracken_dir:
        header.append("Species (Bracken)")
    table = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    warnings = []
    columns = {tsv: table_columns(tsv) for tsv in a.tables}
    for s in a.samples:
        final = totals(read_report(os.path.join(a.final_dir, s + ".taxa.report")))
        if a.host_dir:
            host_rows = read_report(os.path.join(a.host_dir, s + ".taxa.report"))
            trimmed, _ = totals(host_rows)
            host_pct = "%.1f" % (100.0 * sum(clade_reads(host_rows, t) for t in host_taxids) / trimmed) if trimmed else "NA"
        else:
            trimmed, host_pct = "NA", "NA"
        total, classified = final
        pct = 100.0 * classified / total if total else 0.0
        row = [s, trimmed, host_pct, total, "%.1f" % pct, "%.1f" % (100 - pct if total else 0)]
        if a.bracken_dir:
            row.append(n_species(os.path.join(a.bracken_dir, s + ".brackenCorrected_taxa.report")))
        table.append("| " + " | ".join(str(x) for x in row) + " |")
        if total < a.min_reads:
            warnings.append("{}: only {} reads reach the final classification (< {})".format(s, total, a.min_reads))
        if total and pct < a.min_classified_pct:
            warnings.append("{}: {:.1f}% classified (< {}%)".format(s, pct, a.min_classified_pct))
        for tsv in a.tables:
            if not in_table(s, columns[tsv]):
                warnings.append("{}: not found in {}".format(s, os.path.basename(tsv)))
    lines += ["", "## Samples ({})".format(len(a.samples)), ""] + table
    lines += ["", "## Warnings", ""] + (["- " + w for w in warnings] or ["None."])
    lines += ["", "## Output files", ""] + ["- `{}`".format(t) for t in a.tables]
    lines += ["- `{}/` (per-sample classification reports)".format(a.final_dir)]
    if a.bracken_dir:
        lines += ["- `{}/` (per-sample Bracken reports)".format(a.bracken_dir)]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--final-dir", required=True)
    ap.add_argument("--host-dir", default="")
    ap.add_argument("--bracken-dir", default="")
    ap.add_argument("--tables", nargs="*", default=[])
    ap.add_argument("--settings", default="{}", help="JSON dict of settings to list in the report")
    ap.add_argument("--min-reads", type=int, default=100000)
    ap.add_argument("--min-classified-pct", type=float, default=20.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    with open(args.out, "w") as f:
        f.write(build_report(args))
