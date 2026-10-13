#!/usr/bin/env python3
"""Run wgsim once per planned genome and pool the reads into one sample.

One rule rather than a per-genome rule: the plan (and so the genome list) is only
known after `prepare_simulation.py` has measured genome lengths, and each call is
a fraction of a second, so a fan-out here would cost more in scheduling than it
saves. The per-genome seed is derived from the base seed and the row index, so
the pooled sample is reproducible.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", required=True)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--error-rate", type=float, default=0.005)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--workdir", required=True)
    args = p.parse_args()

    if shutil.which("wgsim") is None:
        sys.exit("wgsim not found on PATH (it is in this benchmark's pixi environment)")
    os.makedirs(args.workdir, exist_ok=True)

    with open(args.plan) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{args.plan} has no genomes")

    with open(args.r1, "w") as out1, open(args.r2, "w") as out2:
        for i, row in enumerate(rows):
            r1 = os.path.join(args.workdir, f"{row['id']}.1.fq")
            r2 = os.path.join(args.workdir, f"{row['id']}.2.fq")
            cmd = [
                "wgsim",
                "-e", str(args.error_rate),
                "-r", "0", "-R", "0", "-X", "0",
                "-1", str(args.read_length), "-2", str(args.read_length),
                "-N", row["n_pairs"],
                "-S", str(args.seed + i),
                row["fasta"], r1, r2,
            ]
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                sys.exit(f"wgsim failed for {row['id']}: {res.stderr}")
            for src, dest in ((r1, out1), (r2, out2)):
                with open(src) as fh:
                    shutil.copyfileobj(fh, dest)
                os.remove(src)
            print(f"{row['id']}\t{row['kind']}\t{row['length']} bp\t"
                  f"{row['coverage']}x\t{row['n_pairs']} pairs")

    total = sum(int(r["n_pairs"]) for r in rows)
    print(f"pooled {total} read pairs from {len(rows)} genomes")


if __name__ == "__main__":
    main()
