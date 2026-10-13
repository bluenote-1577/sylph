#!/usr/bin/env python3
"""Summarise how many stage-1 screen k-mers each genome ends up with.

Read from `sylph inspect` YAML by streaming (a 100k-genome database is a big
file and only two fields per genome are wanted, so a YAML parse would be pure
overhead). On an unfloored (`none`) build the `zero` column is the number that
motivates the whole exercise: references with no screen k-mers at all, which the
stage-1 screen can never return however much of them a sample contains.
"""

import argparse
import os
import re

import polars as pl

SUMMARY_FIELDS = ("c", "k", "screen_c", "effective_screen_c", "band_keys",
                  "band_genomes", "num_genomes", "file_size_bytes")
FIELD = re.compile(r"^\s*(\w+):\s*(\S+)\s*$")


def summarise(path):
    counts = []
    header = {}
    with open(path) as fh:
        for line in fh:
            m = FIELD.match(line)
            if not m:
                continue
            key, value = m.group(1), m.group(2)
            if key == "screen_kmers_num":
                counts.append(int(value))
            elif key in SUMMARY_FIELDS and key not in header:
                header[key] = int(value)
    counts.sort()
    n = len(counts)

    def pct(p):
        return counts[min(n - 1, int(p * n))] if n else None

    stem = os.path.basename(path)[: -len(".yaml")]
    db, mode = stem.rsplit(".", 1)
    return {
        "db": db,
        "mode": mode,
        **header,
        "genomes": n,
        "zero": sum(1 for c in counts if c == 0),
        "under_10": sum(1 for c in counts if c < 10),
        "under_50": sum(1 for c in counts if c < 50),
        "min": counts[0] if n else None,
        "p50": pct(0.5),
        "p90": pct(0.9),
        "max": counts[-1] if n else None,
        "total": sum(counts),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("yaml", nargs="+", help="sylph inspect YAML of each .syl2db")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    rows = [summarise(path) for path in args.yaml]
    df = pl.DataFrame(rows, infer_schema_length=None).sort(["db", "mode"])
    df.write_csv(args.output, separator="\t")
    print(df)


if __name__ == "__main__":
    main()
