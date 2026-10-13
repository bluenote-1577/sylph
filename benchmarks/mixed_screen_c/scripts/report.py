#!/usr/bin/env python3
"""Render the collected tables into RESULTS.md.

Deliberately a summary, not a dump: the question is which of the two options to
implement, so the report leads with the per-sample screen cost against the number
of small genomes (where the crossover is), then the build/index cost, then the
sensitivity each arm buys or loses.
"""

import argparse

import polars as pl

pl.Config.set_tbl_rows(500)
pl.Config.set_tbl_cols(30)
pl.Config.set_tbl_width_chars(200)
pl.Config.set_fmt_str_lengths(50)
pl.Config.set_tbl_hide_dataframe_shape(True)
pl.Config.set_tbl_hide_column_data_types(True)
pl.Config.set_tbl_formatting("ASCII_MARKDOWN")


def md(df):
    return str(df) if df.height else "_(no rows)_"


# Columns that are numbers but can come back all-null (and so Null-typed) when a
# measurement was unavailable for every row -- e.g. snakemake reports "NA" for
# peak RSS on jobs too short to sample.
NUMERIC = (
    "wall_s", "cpu_s", "max_rss_mb", "gather_s", "stage1_s", "stage1_only_s",
    "stage2_s", "open_open_s", "open_index_s", "db_bytes", "index_bytes",
    "band_keys", "densified", "survivors", "dense_kept", "effective_screen_c",
    "differing_genomes", "replicate_noise_floor",
    "n", "total_gather_s", "total_stage1_s", "total_stage2_s", "mean_gather_s",
)


def read_table(path):
    df = pl.read_csv(path, separator="\t", infer_schema_length=None)
    return df.with_columns(
        [pl.col(c).cast(pl.Float64, strict=False) for c in NUMERIC if c in df.columns]
    )


def screen_summary(screen):
    """Per (grid point, arm, sample): the fastest replicate of each phase."""
    if "gather_s" not in screen.columns:
        return pl.DataFrame()
    cols = ["n_gtdb", "n_metavr", "sample", "arm"]
    return (
        screen.filter(pl.col("gather_s").is_not_null())
        .group_by(cols)
        .agg(
            pl.col("gather_s").min().round(3).alias("gather_s"),
            # Replicate spread: wall-clock on a shared cluster is noisy, so a
            # difference between arms is only real if it clears this.
            (pl.col("gather_s").max() - pl.col("gather_s").min()).round(3)
                .alias("gather_spread_s"),
            pl.col("stage1_only_s").min().round(3).alias("finalise_s"),
            pl.col("stage2_s").min().round(3).alias("stage2_s"),
            pl.col("survivors").max(),
            pl.col("dense_kept").max(),
            pl.col("wall_s").min().round(1).alias("wall_s"),
            pl.col("max_rss_mb").max().round(0).alias("rss_mb"),
            pl.col("open_open_s").max().round(2).alias("open_s"),
            pl.col("effective_screen_c").max(),
        )
        .sort(cols)
    )


def wall_by_arm(screen):
    """Whole-process wall time, including the `single` (unscreened) arm."""
    return (
        screen.group_by(["n_gtdb", "n_metavr", "sample", "arm"])
        .agg(pl.col("wall_s").min().round(1), pl.col("max_rss_mb").max().round(0))
        .sort(["n_gtdb", "n_metavr", "sample", "arm"])
    )


def headline(screen, build, agreement, reference):
    """The comparison the benchmark exists to make, at its largest grid point.

    Per-sample stage-1 gather time by arm, with each option's speedup over the
    incumbent, plus the sensitivity each arm keeps. Computed, not written down, so
    it cannot drift from the tables below.
    """
    if "gather_s" not in screen.columns or screen.height == 0:
        return ["_(no screen measurements)_"]
    largest = (
        screen.filter(pl.col("gather_s").is_not_null())
        .select("n_gtdb", "n_metavr").unique()
        .sort(["n_metavr", "n_gtdb"], descending=True).row(0, named=True)
    )
    point = screen.filter(
        (pl.col("n_gtdb") == largest["n_gtdb"]) & (pl.col("n_metavr") == largest["n_metavr"])
    )
    fastest = (
        point.group_by(["sample", "arm"]).agg(pl.col("gather_s").min())
        .pivot(on="arm", index="sample", values="gather_s")
        .sort("sample")
    )
    # `single` has no stage-1 at all, so its column is entirely null.
    fastest = fastest.drop(
        [c for c in fastest.columns if c != "sample" and fastest[c].null_count() == fastest.height]
    )
    arms = [c for c in fastest.columns if c != "sample"]
    if "loosen" in arms:
        fastest = fastest.with_columns(
            [(pl.col("loosen") / pl.col(a)).round(1).alias(f"{a}_speedup")
             for a in arms if a != "loosen"]
        )
    index_mb = (
        build.filter(
            (pl.col("n_gtdb") == largest["n_gtdb"]) & (pl.col("n_metavr") == largest["n_metavr"])
        )
        .select("mode", (pl.col("index_bytes") / 1e6).round(1).alias("index_mb"))
        .sort("mode")
    )
    worst_recall = (
        agreement.filter(pl.col("recall_vs_reference").is_not_null())
        .group_by("arm").agg(pl.col("recall_vs_reference").min().round(3))
        .sort("arm")
    )
    # What option B actually charges per sample for its small genomes: the gap
    # between screening them reader-side and not screening them at all.
    marginal = []
    if {"reader", "nofloor"} <= set(arms):
        n_small = point.filter(pl.col("arm") == "reader")["open_small_genomes"].max()
        for row in fastest.iter_rows(named=True):
            delta = (row["reader"] or 0) - (row["nofloor"] or 0)
            marginal.append(
                f"* `{row['sample']}`: +{delta * 1000:.1f} ms over screening no small "
                f"genomes at all ({n_small} of them)"
            )

    return [
        "## Summary",
        "",
        f"At the largest grid point ({largest['n_gtdb']} GTDB representatives + "
        f"{largest['n_metavr']} MetaVR vOTUs), fastest replicate of the stage-1 gather "
        "per sample, in seconds:",
        "",
        md(fastest),
        "",
        *(["Option B's marginal per-sample cost at this point:", ""] + marginal + [""]
          if marginal else []),
        "Stage-1 index size at the same point, by build mode "
        "(`none` is the unfloored baseline both options are measured against):",
        "",
        md(index_mb),
        "",
        f"Worst-case detection recall over every grid point and sample, against the "
        f"`{reference}` arm's profile:",
        "",
        md(worst_recall),
        "",
    ]


def build_summary(build):
    return build.select(
        "n_gtdb", "n_metavr", "mode",
        pl.col("wall_s").round(1).alias("build_wall_s"),
        pl.col("max_rss_mb").round(0).alias("build_rss_mb"),
        (pl.col("index_bytes") / 1e6).round(1).alias("index_mb"),
        (pl.col("db_bytes") / 1e6).round(1).alias("db_mb"),
        "effective_screen_c", "densified", "band_keys",
    ).sort(["n_gtdb", "n_metavr", "mode"])


def multisample_summary(multi):
    """One process, N copies of one sample: does the reader's open cost wash out?

    `wall_s` is the whole `sylph profile db -l list.txt` invocation (one process,
    one database open, N screens); `open_s` is that one-time open cost;
    `total_gather_s` is what the stage-1 pass alone cost, summed over all N
    samples. The gap between `wall_s` and `open_s + total_gather_s +
    total_stage2_s` is real work sylph does not separately instrument: each of
    the N lines in `-l list.txt` is read and bincode-deserialized from disk
    fresh (measured: this, not logging, dominates the gap -- a database whose
    stage-1 survivor *count* differs 3x between arms shows almost the same gap),
    plus the final per-candidate ANI/coverage pass after stage 2. Both are
    per-sample and reference-independent, so they shrink each arm's *relative*
    speedup in `wall_s` versus the pure `gather_s` numbers above, without
    changing which arm is fastest.
    """
    if multi.height == 0:
        return pl.DataFrame()
    return multi.select(
        "n_gtdb", "n_metavr", "sample", "arm", "n",
        pl.col("open_open_s").round(2).alias("open_s"),
        pl.col("total_gather_s").round(2),
        pl.col("mean_gather_s").round(4),
        pl.col("total_stage2_s").round(2),
        pl.col("wall_s").round(1),
        pl.col("max_rss_mb").round(0).alias("rss_mb"),
        (pl.col("wall_s") / pl.col("n")).round(3).alias("wall_s_per_sample"),
    ).sort(["n_gtdb", "n_metavr", "sample", "arm"])


def multisample_headline(multi):
    """Whole-process wall time by arm, at whatever N the table was run at --
    the number that answers "does the reader's fixed open cost pay for itself
    once samples stop being scored one at a time".
    """
    if multi.height == 0:
        return ["_(no multi-sample measurements)_"]
    lines = []
    for (ng, nv, sample, n), group in multi.group_by(
        ["n_gtdb", "n_metavr", "sample", "n"], maintain_order=True
    ):
        pivot = (
            group.select("arm", "wall_s").sort("arm")
            .pivot(on="arm", index=None, values="wall_s")
        )
        arms = pivot.columns
        row = pivot.row(0, named=True)
        lines.append(
            f"* {ng} GTDB + {nv} MetaVR, `{sample}` x {int(n)}: "
            + ", ".join(f"{a}={row[a]:.1f}s" for a in arms if row[a] is not None)
        )
    return lines


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build", required=True)
    p.add_argument("--screen", required=True)
    p.add_argument("--agreement", required=True)
    p.add_argument("--recall", required=True)
    p.add_argument("--stability", required=True)
    p.add_argument("--screen-kmers", required=True)
    p.add_argument("--multisample", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    build, screen = read_table(args.build), read_table(args.screen)
    agreement = read_table(args.agreement)
    recall, kmers = read_table(args.recall), read_table(args.screen_kmers)
    stability = read_table(args.stability)
    multi = read_table(args.multisample)

    reference = (agreement["reference_arm"][0] if "reference_arm" in agreement.columns
                 and agreement.height else "?")

    parts = [
        "# Mixed screen-c scalability: loosen vs value-banded tiers vs reader-side screening",
        "",
        "Generated by `scripts/report.py`; see the Snakefile docstring for the design.",
        "Arms: **loosen** (incumbent, widen the pooled MPHF), **band** (option A,",
        "value-banded tiers), **reader** (option B, `--screen-small-genomes` on an",
        "unfloored build), **nofloor** (no small-genome handling at all), **single**",
        "(no two-stage screen: every genome profiled densely).",
        "",
        *headline(screen, build, agreement, reference),
        "## Stage-1 screen cost",
        "",
        "`gather_s` is the inverted pass over the sample plus, for `reader`, the",
        "small-genome intersection -- the number the options exist to reduce.",
        "`finalise_s` is the per-survivor ANI pass, `stage2_s` the dense decode.",
        "",
        md(screen_summary(screen)),
        "",
        "## Whole-process cost, all arms",
        "",
        md(wall_by_arm(screen)),
        "",
        "## Amortized cost over many samples",
        "",
        "One `sylph profile db -l list.txt` process, N copies of one sample: the",
        "database opens once and screens N times, so this is where a one-time open",
        "cost (`band`'s pooled-index deserialize, `reader`'s dense-block read-back)",
        "either washes out against N cheap per-sample screens, or doesn't.",
        "",
        *multisample_headline(multi),
        "",
        md(multisample_summary(multi)),
        "",
        "## Build and index cost",
        "",
        md(build_summary(build)),
        "",
        "## Screen k-mers per genome, as stored",
        "",
        "`zero` is the number of references the nominal `--screen-c` leaves with no",
        "stage-1 k-mers at all (invisible to the screen); on `loosen`/`band` builds the",
        "floor removes them.",
        "",
        md(kmers.select([c for c in ("db", "mode", "genomes", "zero", "under_10",
                                     "under_50", "p50", "p90", "screen_c",
                                     "effective_screen_c", "band_keys", "total")
                         if c in kmers.columns])),
        "",
        "## Profile stability",
        "",
        "The three equivalent arms are asserted to pass exactly the same stage-1",
        "survivors (`compare_profiles.py` fails the run otherwise). The *final*",
        "profile is not a function of the screen alone: `profile` gives each shared",
        "k-mer to one genome, and near-identical references make that a tie which",
        "thread scheduling breaks differently each run. `replicate_noise_floor` is",
        "the largest same-arm, different-replicate disagreement, i.e. the resolution",
        "limit of any between-arm comparison here.",
        "",
        md(stability),
        "",
        f"## Detections vs the `{reference}` arm",
        "",
        md(agreement.sort(["db", "sample", "arm"])),
        "",
        "## Truth-based recall on the simulated sample",
        "",
        md(recall.sort(["db", "kind", "coverage", "arm"])),
        "",
    ]
    with open(args.output, "w") as out:
        out.write("\n".join(parts))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
