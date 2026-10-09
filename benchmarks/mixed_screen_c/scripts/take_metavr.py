#!/usr/bin/env python3
"""Write the first N MetaVR vOTU representatives, round-robined across chunks.

Round-robin rather than contiguous blocks so every chunk gets a similar total
sequence length -- `sylph sketch` parallelises over input *files*, so uneven
chunks leave cores idle at the end.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metavr_stream import iter_records  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shards", required=True, help="directory of *.fna shards")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--chunks", type=int, required=True)
    p.add_argument("--out-prefix", required=True, help="writes {prefix}_{i}.fna")
    args = p.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)) or ".", exist_ok=True)
    handles = [open(f"{args.out_prefix}_{i}.fna", "w") for i in range(args.chunks)]
    total = 0
    bases = 0
    try:
        for i, (header, seq) in enumerate(iter_records(args.shards, limit=args.n)):
            fh = handles[i % args.chunks]
            fh.write(f">{header}\n")
            fh.write("\n".join(seq) + "\n")
            total += 1
            bases += sum(len(s) for s in seq)
    finally:
        for fh in handles:
            fh.close()
    if total < args.n:
        sys.exit(f"only {total} records available in {args.shards}, needed {args.n}")
    print(f"wrote {total} records ({bases} bases) across {args.chunks} chunks")


if __name__ == "__main__":
    main()
