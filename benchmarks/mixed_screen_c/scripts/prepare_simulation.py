#!/usr/bin/env python3
"""Choose the genomes the simulated sample is built from, and plan their coverage.

The real metagenome sketches measure cost but carry no ground truth, so screen
*sensitivity* is measured on a sample simulated from references that are known to
be in the database. Two properties matter:

  * The chosen genomes are the *first* n of each subset stream, and both subset
    streams are nested by construction (see `subset_lines.py` and
    `metavr_stream.py`), so the same truth set is present at every non-empty
    grid point.
  * Viral coverages are spread over a range, because a screen with only a couple
    of k-mers per genome fails gradually: it is partial containment, not absence,
    that a too-sparse screen turns into a false negative.
"""

import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from metavr_stream import iter_records  # noqa: E402

COLUMNS = [
    "id",
    "kind",
    "fasta",
    "length",
    "coverage",
    "n_pairs",
    "contig_full",
    "contig_id",
]


def n_pairs(coverage, length, read_length):
    """Read pairs needed for `coverage`-fold coverage; at least one."""
    return max(1, round(coverage * length / (2 * read_length)))


def write_fasta(path, header, seq_lines):
    with open(path, "w") as out:
        out.write(f">{header}\n")
        out.write("\n".join(seq_lines) + "\n")


def gunzip_fasta(src, dest):
    """Decompress a GTDB representative (wgsim cannot read gzip) and measure it."""
    length = 0
    first_header = None
    with gzip.open(src, "rt") as fh, open(dest, "w") as out:
        for line in fh:
            out.write(line)
            if line.startswith(">"):
                if first_header is None:
                    first_header = line[1:].strip()
            else:
                length += len(line.strip())
    return first_header, length


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shards", required=True)
    p.add_argument("--n-viral", type=int, required=True)
    p.add_argument("--viral-coverages", required=True, help="comma separated")
    p.add_argument("--gtdb-list", default="", help="empty if no GTDB genomes are in the grid")
    p.add_argument("--n-bacterial", type=int, default=0)
    p.add_argument("--bacterial-coverage", type=float, default=0.5)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--outdir", required=True)
    p.add_argument("--plan", required=True)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    covs = [float(c) for c in args.viral_coverages.split(",") if c]
    rows = []

    # Viral: cycle the coverage levels so each level gets an equal share of
    # genomes, giving a recall-vs-coverage curve per arm.
    for i, (header, seq) in enumerate(iter_records(args.shards, limit=args.n_viral)):
        contig_id = header.split()[0] if header.split() else header
        fasta = os.path.join(args.outdir, f"viral_{i:05d}.fna")
        write_fasta(fasta, header, seq)
        length = sum(len(s) for s in seq)
        cov = covs[i % len(covs)]
        rows.append(
            dict(
                id=f"viral_{i:05d}",
                kind="viral",
                fasta=fasta,
                length=length,
                coverage=cov,
                n_pairs=n_pairs(cov, length, args.read_length),
                contig_full=header,
                contig_id=contig_id,
            )
        )
    if len(rows) < args.n_viral:
        sys.exit(f"only {len(rows)} MetaVR records available, needed {args.n_viral}")

    # Bacterial: the large-genome positives, at one (low) coverage.
    if args.gtdb_list and args.n_bacterial > 0:
        with open(args.gtdb_list) as fh:
            paths = [ln.strip() for ln in fh if ln.strip()][: args.n_bacterial]
        if len(paths) < args.n_bacterial:
            sys.exit(f"{args.gtdb_list} has fewer than {args.n_bacterial} genomes")
        for i, src in enumerate(paths):
            dest = os.path.join(args.outdir, f"bacterial_{i:05d}.fna")
            header, length = gunzip_fasta(src, dest)
            cov = args.bacterial_coverage
            rows.append(
                dict(
                    id=f"bacterial_{i:05d}",
                    kind="bacterial",
                    fasta=dest,
                    length=length,
                    coverage=cov,
                    n_pairs=n_pairs(cov, length, args.read_length),
                    contig_full=header,
                    contig_id=header.split()[0] if header else "",
                )
            )

    with open(args.plan, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")

    viral = [r for r in rows if r["kind"] == "viral"]
    print(f"planned {len(rows)} genomes: {len(viral)} viral "
          f"(median length {sorted(r['length'] for r in viral)[len(viral) // 2] if viral else 0} bp), "
          f"{len(rows) - len(viral)} bacterial")
    print(f"total simulated bases: {sum(r['n_pairs'] * 2 * args.read_length for r in rows)}")


if __name__ == "__main__":
    main()
