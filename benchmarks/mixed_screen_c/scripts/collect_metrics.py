#!/usr/bin/env python3
"""Collect the benchmark's cost measurements into two tidy TSVs.

Three sources, because no single one covers it:

  * `bench_build` / `bench_open` / `bench_screen` lines that sylph itself logs --
    the only place the *internal* breakdown is visible (index bytes, band size,
    the stage-1 inverted pass separately from the stage-2 dense decode).
  * snakemake `benchmark:` files -- whole-process wall time and peak RSS, which
    sylph cannot report about itself.
  * the databases on disk -- their actual size.

File names carry the grid point, so the wildcards are parsed back out of the
paths rather than threaded through as extra inputs.
"""

import argparse
import glob
import os
import re

import polars as pl

# `key=value` pairs, values possibly containing '/' or '.' (paths) but no spaces.
KV = re.compile(r"(\w+)=(\S+)")

INT_FIELDS = {
    "genomes", "densified", "band_keys", "band_owners", "dense_bytes", "index_bytes",
    "footer_bytes", "screen_c", "effective_screen_c", "min_sparse_kmers", "survivors",
    "dense_kept", "sample_kmers", "small_genomes", "small_kmers",
}
FLOAT_FIELDS = {"gather_s", "stage1_s", "stage2_s", "index_s", "open_s"}


def parse_bench_lines(path, tag):
    """Every `tag key=value ...` line in a log, as dicts with numbers typed."""
    out = []
    with open(path) as fh:
        for line in fh:
            idx = line.find(tag + " ")
            if idx < 0:
                continue
            row = {}
            for key, value in KV.findall(line[idx + len(tag):]):
                if key in INT_FIELDS:
                    row[key] = int(value)
                elif key in FLOAT_FIELDS:
                    row[key] = float(value)
                else:
                    row[key] = value
            out.append(row)
    return out


def read_snakemake_benchmark(path):
    """Wall seconds and peak RSS (MB) from a snakemake benchmark file.

    Snakemake writes "NA" for the resource columns when a job finishes before it
    can sample /proc, which is common for the smallest grid points.
    """
    if not os.path.exists(path):
        return {}
    df = pl.read_csv(path, separator="\t", infer_schema_length=None)
    if df.height == 0:
        return {}
    row = df.row(0, named=True)

    def num(key):
        value = row.get(key)
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    return {"wall_s": num("s"), "max_rss_mb": num("max_rss"), "cpu_s": num("cpu_time")}


def grid_columns(db_token):
    """`g2000_v10000` -> the two axis values."""
    m = re.fullmatch(r"g(\d+)_v(\d+)", db_token)
    if not m:
        raise ValueError(f"cannot parse grid point from {db_token!r}")
    return {"n_gtdb": int(m.group(1)), "n_metavr": int(m.group(2))}


def collect_build(logdir, benchdir, dbdir):
    rows = []
    for log in sorted(glob.glob(os.path.join(logdir, "convert", "*.log"))):
        stem = os.path.basename(log)[: -len(".log")]
        db_token, mode = stem.rsplit(".", 1)
        builds = parse_bench_lines(log, "bench_build")
        if not builds:
            raise SystemExit(f"{log} has no bench_build line: did the conversion fail?")
        # `bench_build` also carries `mode`; the log file name is authoritative.
        row = {"db": db_token, **grid_columns(db_token), **builds[-1], "mode": mode}
        row.update(read_snakemake_benchmark(os.path.join(benchdir, "convert", stem + ".tsv")))
        db_path = os.path.join(dbdir, stem + ".syl2db")
        row["db_bytes"] = os.path.getsize(db_path) if os.path.exists(db_path) else None
        rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=None)


def collect_screen(logdir, benchdir):
    rows = []
    for log in sorted(glob.glob(os.path.join(logdir, "profile", "*.log"))):
        stem = os.path.basename(log)[: -len(".log")]
        db_token, arm, sample, rep = stem.split(".")
        base = dict(
            db=db_token, arm=arm, sample=sample, rep=int(rep.removeprefix("rep")),
            **grid_columns(db_token),
        )
        base.update(read_snakemake_benchmark(os.path.join(benchdir, "profile", stem + ".tsv")))
        opens = parse_bench_lines(log, "bench_open")
        screens = parse_bench_lines(log, "bench_screen")
        if arm == "single":
            # No screen at all: only the whole-process numbers are meaningful.
            rows.append(base)
            continue
        if not screens:
            raise SystemExit(f"{log} has no bench_screen line: did profiling fail?")
        for screen in screens:
            row = dict(base)
            if opens:
                row.update({f"open_{k}": v for k, v in opens[-1].items() if k != "db"})
            # sylph names the sample by its sketch's internal file_name and the
            # database by path; the file name's wildcards are what the grid is
            # indexed by, so keep those and rename sylph's.
            screen = dict(screen)
            screen["sample_file"] = screen.pop("sample", None)
            screen["db_path"] = screen.pop("db", None)
            row.update(screen)
            # The dense stage is reported as a delta so the three phases add up.
            row["stage1_only_s"] = row["stage1_s"] - row["gather_s"]
            rows.append(row)
    return pl.DataFrame(rows, infer_schema_length=None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--logdir", required=True)
    p.add_argument("--benchdir", required=True)
    p.add_argument("--dbdir", required=True)
    p.add_argument("--build-out", required=True)
    p.add_argument("--screen-out", required=True)
    args = p.parse_args()

    build = collect_build(args.logdir, args.benchdir, args.dbdir)
    screen = collect_screen(args.logdir, args.benchdir)
    build.write_csv(args.build_out, separator="\t")
    screen.write_csv(args.screen_out, separator="\t")
    print(f"{build.height} build rows -> {args.build_out}")
    print(f"{screen.height} screen rows -> {args.screen_out}")


if __name__ == "__main__":
    main()
