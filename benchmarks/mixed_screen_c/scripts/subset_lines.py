#!/usr/bin/env python3
"""Deterministically take the first `n` lines of a seeded shuffle of a file.

Shuffling once with a fixed seed and then taking prefixes keeps the subsets
*nested*: the n=500 subset is contained in the n=2000 subset. The benchmark
relies on that, because the genomes it simulates reads from must be present at
every non-empty grid point.
"""

import argparse
import random
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    with open(args.input) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    if args.n > len(lines):
        sys.exit(f"asked for {args.n} lines but {args.input} has only {len(lines)}")
    random.Random(args.seed).shuffle(lines)
    with open(args.output, "w") as out:
        out.write("\n".join(lines[: args.n]) + "\n")
    print(f"wrote {args.n} of {len(lines)} lines to {args.output}")


if __name__ == "__main__":
    main()
