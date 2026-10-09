#!/usr/bin/env python3
"""Check the arms agree where they must, and measure what the others lose.

Three jobs, deliberately in one place because they use the same tables:

1. **Equivalence (a hard failure).** `loosen`, `band` and the reader-side screen
   are three ways of storing the *same* per-genome screen k-mer set, so they must
   produce the *same stage-1 survivor set* (`--screen-dump`). If they do not, one
   of them is wrong and the timing comparison is meaningless, so this exits
   non-zero and prints the differing genomes rather than writing a report that
   looks fine.

   The assertion is made on the survivor set and not on the final profile because
   the profile is not a function of the screen alone: `profile` assigns each
   shared k-mer to one genome, and a reference with many near-identical members
   (12.7 M MetaVR vOTUs) has ties that thread scheduling breaks differently from
   run to run. Two replicates of the *same* arm on the same database then differ
   by tens of genomes, so profile-level equality is not a property any arm can
   have at that size (measured: see the stability table below).
2. **Stability (a measurement, and a bounded failure).** Between-arm differences
   in the final profile, next to the same-arm between-replicate differences that
   bound what any comparison can resolve. Arms are only flagged when they differ
   by more than that noise floor.
3. **Sensitivity (a measurement).** Against the unscreened single-stage profile
   of the same references, and against the simulated sample's known truth broken
   down by genome kind and coverage, so "the screen loses small genomes" becomes
   a number per arm.
"""

import argparse
import glob
import os
import sys

import polars as pl

# Genomes are compared on the first token of the contig name: sylph stores the
# whole FASTA header line, and GTDB headers carry a description after the
# accession while the plan may have either form.
def first_token(col):
    return col.str.strip_chars().str.split(" ").list.first()


def parse_name(path):
    """`g2000_v10000.band.marine.rep1.tsv` -> its four wildcards."""
    stem = os.path.basename(path)[: -len(".tsv")]
    db, arm, sample, rep = stem.split(".")
    return db, arm, sample, int(rep.removeprefix("rep"))


def load_profiles(profile_dir):
    """All profile rows, plus the (db, arm, sample, rep) runs that produced them.

    The run list is kept separately because an arm that detected nothing
    contributes no rows -- and "detected nothing" is exactly one of the outcomes
    being measured, so it must not drop out of the tables.
    """
    frames = []
    runs = []
    for path in sorted(glob.glob(os.path.join(profile_dir, "*.tsv"))):
        db, arm, sample, rep = parse_name(path)
        runs.append((db, arm, sample, rep))
        df = pl.read_csv(path, separator="\t", infer_schema_length=None)
        keep = [c for c in ("Contig_name", "Adjusted_ANI", "Taxonomic_abundance",
                            "Genome_file", "Containment_ind") if c in df.columns]
        # An arm that detected nothing writes a header-only file, whose columns
        # come back as Null-typed; the string/float ops below need real dtypes.
        if df.height == 0:
            df = pl.DataFrame(schema={
                "Contig_name": pl.String, "Adjusted_ANI": pl.Float64,
                "Taxonomic_abundance": pl.Float64,
            })
            keep = list(df.columns)
        df = df.select(keep).with_columns(
            genome=first_token(pl.col("Contig_name").cast(pl.String)),
            db=pl.lit(db), arm=pl.lit(arm), sample=pl.lit(sample), rep=pl.lit(rep),
        )
        frames.append(df)
    if not frames:
        sys.exit(f"no profile TSVs found in {profile_dir}")
    return pl.concat(frames, how="diagonal_relaxed"), runs


def detected(profiles, db, arm, sample, rep=1):
    sel = profiles.filter(
        (pl.col("db") == db) & (pl.col("arm") == arm)
        & (pl.col("sample") == sample) & (pl.col("rep") == rep)
    )
    return set(sel["genome"].to_list())


def load_survivors(screen_dir):
    """`--screen-dump` survivor sets, keyed by (db, arm, sample, rep)."""
    out = {}
    for path in sorted(glob.glob(os.path.join(screen_dir, "*.tsv"))):
        df = pl.read_csv(path, separator="\t", infer_schema_length=None)
        out[parse_name(path)] = set(
            df.select(first_token(pl.col("Contig_name").cast(pl.String)))
            .to_series().to_list()
        )
    return out


def fail(header, problems):
    print(header, file=sys.stderr)
    for line in problems:
        print("  " + line, file=sys.stderr)
    sys.exit(1)


def describe_diff(label_a, set_a, label_b, set_b):
    only_a, only_b = set_a - set_b, set_b - set_a
    return (f"{len(only_a)} only in {label_a} (e.g. {sorted(only_a)[:3]}), "
            f"{len(only_b)} only in {label_b} (e.g. {sorted(only_b)[:3]})")


def check_survivor_equivalence(survivors, equivalent):
    """Exit non-zero if arms that store the same screen pass different genomes.

    This is the real invariant: the stage-1 survivor set is a function of the
    screen k-mer set alone, so storing that set three ways must not change it.
    """
    base_arm = equivalent[0]
    keys = {(db, sample, rep) for (db, _arm, sample, rep) in survivors}
    problems = []
    checked = 0
    for db, sample, rep in sorted(keys):
        if (db, base_arm, sample, rep) not in survivors:
            continue
        base = survivors[(db, base_arm, sample, rep)]
        for arm in equivalent[1:]:
            other = survivors.get((db, arm, sample, rep))
            if other is None:
                continue
            checked += 1
            if base != other:
                problems.append(f"{db} {sample} rep{rep}: {base_arm} vs {arm} -- "
                                + describe_diff(base_arm, base, arm, other))
    if problems:
        fail("EQUIVALENCE FAILURE: arms that store the same screen pass different "
             "stage-1 survivors:", problems)
    if not checked:
        return False
    print(f"stage-1 survivor sets identical for {equivalent} across "
          f"{checked} arm comparison(s)")
    # Stage 1 has no tie-breaking in it, so it should also be reproducible.
    for db, sample in sorted({(k[0], k[2]) for k in survivors}):
        for arm in sorted({k[1] for k in survivors}):
            reps = sorted(r for (d, a, s, r) in survivors if (d, a, s) == (db, arm, sample))
            first = survivors.get((db, arm, sample, reps[0])) if reps else None
            for rep in reps[1:]:
                if survivors[(db, arm, sample, rep)] != first:
                    print(f"WARNING: {db} {arm} {sample} survivors differ between "
                          f"rep{reps[0]} and rep{rep}", file=sys.stderr)
    return True


def stability_table(profiles, runs, equivalent, survivors_checked):
    """Between-arm profile differences against the same-arm replicate noise floor.

    `profile`'s output is not a function of the screen alone: shared k-mers are
    assigned to one genome, and a reference with near-identical members has ties
    that get broken differently from run to run. The floor measures that, so a
    between-arm difference can be read against something. With a deterministic
    reference the floor is 0 and this reduces to exact equality.
    """
    rows = []
    problems = []
    base_arm = equivalent[0]
    for db in sorted({r[0] for r in runs}):
        for sample in sorted({r[2] for r in runs}):
            reps = sorted({r[3] for r in runs if r[0] == db and r[2] == sample})
            sets = {(arm, rep): detected(profiles, db, arm, sample, rep)
                    for arm in equivalent for rep in reps}
            # Largest same-arm, different-replicate disagreement: the resolution
            # limit of any comparison on this (database, sample).
            floor = 0
            for arm in equivalent:
                for i, rep_a in enumerate(reps):
                    for rep_b in reps[i + 1:]:
                        floor = max(floor, len(sets[(arm, rep_a)] ^ sets[(arm, rep_b)]))
            for arm in equivalent[1:]:
                for rep in reps:
                    between = len(sets[(base_arm, rep)] ^ sets[(arm, rep)])
                    rows.append({
                        "db": db, "sample": sample, "arm": arm, "vs": base_arm, "rep": rep,
                        "n_detected": len(sets[(arm, rep)]),
                        "differing_genomes": between,
                        "replicate_noise_floor": floor,
                    })
                    # 50 % slack over the observed floor: it is itself an estimate
                    # from a handful of replicates. At floor 0 this is exact.
                    if between > floor * 1.5:
                        problems.append(
                            f"{db} {sample} rep{rep}: {base_arm} vs {arm} differ by {between} "
                            f"genomes, more than the {floor}-genome replicate noise floor -- "
                            + describe_diff(base_arm, sets[(base_arm, rep)], arm, sets[(arm, rep)])
                        )
    if problems:
        header = ("PROFILE DIFFERENCE beyond replicate noise between arms that store the "
                  "same screen:")
        if survivors_checked:
            # The survivor sets were checked directly and agreed, so this is a
            # downstream (reassignment) difference: report it, do not fail on it.
            print(header, file=sys.stderr)
            for line in problems:
                print("  " + line, file=sys.stderr)
        else:
            fail(header, problems)
    return pl.DataFrame(rows, infer_schema_length=None)


def numeric_drift(profiles, equivalent):
    """Largest per-genome ANI/abundance difference between equivalent arms.

    Not a failure: the detected sets being identical is the invariant, while the
    reported numbers can differ in the last bits when k-mer reassignment breaks a
    tie in a different order.
    """
    base_arm = equivalent[0]
    rows = []
    for arm in equivalent[1:]:
        joined = (
            profiles.filter((pl.col("arm") == base_arm) & (pl.col("rep") == 1))
            .join(
                profiles.filter((pl.col("arm") == arm) & (pl.col("rep") == 1)),
                on=["db", "sample", "genome"], how="inner", suffix="_other",
            )
        )
        for col in ("Adjusted_ANI", "Taxonomic_abundance"):
            if col in joined.columns and f"{col}_other" in joined.columns:
                diff = (joined[col].cast(pl.Float64) - joined[f"{col}_other"].cast(pl.Float64)).abs()
                rows.append({
                    "arm": arm, "vs": base_arm, "column": col,
                    "max_abs_diff": float(diff.max() or 0.0), "n_compared": joined.height,
                })
    return pl.DataFrame(rows, infer_schema_length=None)


def agreement_table(profiles, runs, reference_arm):
    rows = []
    arms = sorted({r[1] for r in runs})
    for db in sorted({r[0] for r in runs}):
        for sample in sorted({r[2] for r in runs}):
            ref = detected(profiles, db, reference_arm, sample)
            for arm in arms:
                got = detected(profiles, db, arm, sample)
                rows.append({
                    "db": db, "sample": sample, "arm": arm,
                    "n_detected": len(got),
                    "n_reference": len(ref),
                    "shared": len(got & ref),
                    "missed_vs_reference": len(ref - got),
                    "extra_vs_reference": len(got - ref),
                    "recall_vs_reference": (len(got & ref) / len(ref)) if ref else None,
                })
    return pl.DataFrame(rows, infer_schema_length=None)


def recall_table(profiles, runs, plan_path):
    """Truth-based recall on the simulated sample, by genome kind and coverage."""
    plan = pl.read_csv(plan_path, separator="\t", infer_schema_length=None).with_columns(
        genome=first_token(pl.col("contig_full"))
    )
    rows = []
    dbs = sorted({r[0] for r in runs})
    arms = sorted({r[1] for r in runs})
    for db in dbs:
        for arm in arms:
            got = detected(profiles, db, arm, "simulated")
            for (kind, coverage), group in plan.group_by(["kind", "coverage"], maintain_order=True):
                truth = set(group["genome"].to_list())
                found = truth & got
                rows.append({
                    "db": db, "arm": arm, "kind": kind, "coverage": coverage,
                    "n_truth": len(truth), "n_found": len(found),
                    "recall": len(found) / len(truth) if truth else None,
                })
    return pl.DataFrame(rows, infer_schema_length=None)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile-dir", required=True)
    p.add_argument("--screen-dir", help="`--screen-dump` survivor TSVs, if collected")
    p.add_argument("--plan", required=True)
    p.add_argument("--equivalent", required=True, help="comma separated arm names")
    p.add_argument("--reference-arm", default="single")
    p.add_argument("--agreement-out", required=True)
    p.add_argument("--recall-out", required=True)
    p.add_argument("--stability-out", required=True)
    args = p.parse_args()

    profiles, runs = load_profiles(args.profile_dir)
    equivalent = [a for a in args.equivalent.split(",") if a]
    present = {r[1] for r in runs}
    missing = [a for a in equivalent if a not in present]
    if missing:
        sys.exit(f"arms {missing} have no profiles in {args.profile_dir}")

    survivors = load_survivors(args.screen_dir) if args.screen_dir else {}
    survivors_checked = check_survivor_equivalence(survivors, equivalent) if survivors else False
    if not survivors_checked:
        print("no stage-1 survivor dumps found; asserting equivalence on the profiles instead")
    stability = stability_table(profiles, runs, equivalent, survivors_checked)
    stability.write_csv(args.stability_out, separator="\t")
    drift = numeric_drift(profiles, equivalent)
    if drift.height:
        print("numeric drift between equivalent arms (informational):")
        print(drift)

    reference = args.reference_arm if args.reference_arm in present else equivalent[0]
    if reference != args.reference_arm:
        print(f"reference arm {args.reference_arm} not profiled; using {reference}")
    agreement = agreement_table(profiles, runs, reference)
    agreement = agreement.with_columns(reference_arm=pl.lit(reference))
    agreement.write_csv(args.agreement_out, separator="\t")
    recall_table(profiles, runs, args.plan).write_csv(args.recall_out, separator="\t")
    print(f"wrote {args.agreement_out} and {args.recall_out}")


if __name__ == "__main__":
    main()
