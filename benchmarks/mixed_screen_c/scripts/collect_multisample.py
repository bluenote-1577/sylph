#!/usr/bin/env python3
"""Collect the multi-sample amortization measurements into one tidy TSV.

Same three-source approach as `collect_metrics.py` (sylph's own `bench_*` debug
lines, snakemake's `benchmark:` file, nothing from disk here), but per *process*
rather than per (database, arm, sample, rep): a multi-sample run is one `sylph
profile db -l list.txt` invocation covering N copies of one sample, opening the
database exactly once and screening it N times, so the interesting number is not
any single `bench_screen` line but the *sum* over all of them next to the *one*
`bench_open` line.
"""

import argparse
import glob
import os
import re
import sys

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_metrics import parse_bench_lines, read_snakemake_benchmark, grid_columns  # noqa: E402


def parse_name(stem):
    """`g2000_v100000.reader.marine.n100` -> its wildcards."""
    m = re.fullmatch(r"(g\d+_v\d+)\.([^.]+)\.([^.]+)\.n(\d+)", stem)
    if not m:
        raise ValueError(f"cannot parse multisample run from {stem!r}")
    db, arm, sample, n = m.groups()
    return db, arm, sample, int(n)


def collect(logdir, benchdir):
    rows = []
    for log in sorted(glob.glob(os.path.join(logdir, "multisample", "*.log"))):
        stem = os.path.basename(log)[: -len(".log")]
        db, arm, sample, n = parse_name(stem)
        opens = parse_bench_lines(log, "bench_open")
        screens = parse_bench_lines(log, "bench_screen")
        if not opens:
            raise SystemExit(f"{log} has no bench_open line: did profiling fail?")
        if len(screens) != n:
            raise SystemExit(
                f"{log} has {len(screens)} bench_screen line(s), expected n={n}: "
                "did profiling fail partway through the sample list?"
            )
        row = dict(db=db, arm=arm, sample=sample, n=n, **grid_columns(db))
        row.update({f"open_{k}": v for k, v in opens[-1].items() if k != "db"})
        # The per-process total each phase actually costs across all N samples --
        # what a real "score N samples against this database" job would pay.
        row["total_gather_s"] = sum(s["gather_s"] for s in screens)
        row["total_stage1_s"] = sum(s["stage1_s"] for s in screens)
        row["total_stage2_s"] = sum(s["stage2_s"] for s in screens)
        row["mean_gather_s"] = row["total_gather_s"] / n
        row.update(read_snakemake_benchmark(
            os.path.join(benchdir, "multisample", stem + ".tsv")
        ))
        rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--logdir", required=True)
    p.add_argument("--benchdir", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    df = collect(args.logdir, args.benchdir)
    df.write_csv(args.output, separator="\t")
    print(f"{df.height} multisample rows -> {args.output}")


if __name__ == "__main__":
    main()
