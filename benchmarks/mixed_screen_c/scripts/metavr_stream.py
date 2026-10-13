#!/usr/bin/env python3
"""Streaming reader for the MetaVR vOTU representative FASTA shards.

Shared by `take_metavr.py` and `prepare_simulation.py` so both walk the shards in
exactly the same order: shard files sorted by name, records in file order. That
makes "the first N records" a well-defined, reproducible set, and makes the
subsets nested across N.
"""

import glob
import os
import sys


def shard_files(shards_dir):
    files = sorted(glob.glob(os.path.join(shards_dir, "*.fna")))
    if not files:
        sys.exit(f"no *.fna shards found in {shards_dir}")
    return files


def iter_records(shards_dir, limit=None):
    """Yield `(header, [sequence lines])` for the first `limit` records."""
    n = 0
    for path in shard_files(shards_dir):
        with open(path) as fh:
            header, seq = None, []
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if header is not None:
                        yield header, seq
                        n += 1
                        if limit is not None and n >= limit:
                            return
                    header, seq = line[1:], []
                elif line:
                    seq.append(line)
            if header is not None:
                yield header, seq
                n += 1
                if limit is not None and n >= limit:
                    return
