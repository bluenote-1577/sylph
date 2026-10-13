# Mixed screen-c scalability benchmark

Does a `.syl2db` stage-1 screen that mixes several `screen_c` values (denser for
small genomes, coarser for large ones) have to cost what the incumbent
`--small-genome-screen loosen` costs? This harness measures the two candidate
answers against it, on GTDB representatives (large genomes) with MetaVR vOTU
representatives (small genomes) mixed in at increasing numbers:

| arm | build | profile | what it is |
| --- | --- | --- | --- |
| `loosen` | `--small-genome-screen loosen` | — | incumbent: widen the pooled index to the densest rate any genome needed |
| `band` | `--small-genome-screen band` | — | **option A**: value-banded tiers; pooled index stays at `--screen-c` |
| `reader` | `--small-genome-screen none` | `--screen-small-genomes 50` | **option B**: no index change; screen small genomes directly at open time |
| `nofloor` | `--small-genome-screen none` | — | no small-genome handling at all: the speed ceiling and the sensitivity it costs |
| `single` | (plain `.syldb`) | — | no two-stage screen: what a perfect screen would let through |

`loosen`, `band` and `reader` screen the *same* per-genome k-mer set, so they must
pass the same genomes out of stage 1; `scripts/compare_profiles.py` fails the run
if they do not. That makes this pipeline a correctness regression test as well as a
benchmark. It compares the `--screen-dump` survivor sets rather than the final
profiles, because a profile also depends on how reassignment breaks ties between
near-identical references, and at MetaVR scale that differs between two runs of the
same arm (quantified in the report's stability table).

See the Snakefile docstring for the design, `config.yaml` for every knob, and
`RESULTS.md` / `RESULTS_FULL.md` (generated) for the measurements.

## Running

```bash
pixi run snakemake --configfile config_smoke.yaml --cores 4       # smoke test, ~1 min
pixi run snakemake --profile aqua                                 # the size sweep
pixi run snakemake --profile aqua --configfile config_full.yaml   # all of GTDB R226 + MetaVR
```

The smoke test runs every rule and every assertion on a handful of genomes, with
all of its output under `smoke_work/` so it cannot be confused with a real run.

The sweep answers where the crossover between the options is; the full run answers
whether the winner survives the size the MetaVR scan actually needs (12.8 M
references), which is where option B's O(small references) cost and option A's
band-1 lookup are most exposed. It writes to its own workdir, `results_full/` and
`RESULTS_FULL.md`.
