#!/usr/bin/env python3
"""Write a newline-delimited file listing one sample sketch N times.

`sylph profile db -l list.txt` treats each line as an independent sample -- the
database is opened once regardless -- so this is how the multi-sample
amortization benchmark gets "100 samples" without re-sketching or re-fetching
anything: the same sketch, scored N times.
"""

import argparse


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True, help="path to a *.sylsp sketch")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    with open(args.output, "w") as out:
        for _ in range(args.n):
            out.write(args.sample + "\n")
    print(f"wrote {args.n} copies of {args.sample} to {args.output}")


if __name__ == "__main__":
    main()
