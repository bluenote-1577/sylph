# Unreleased

- `convert-db-two-screen --small-genome-screen` chooses how a `.syl2db` stores the mixture of stage-1 screen rates that `--min-sparse-kmers` creates for small (viral/plasmid/short-contig) genomes. All three modes screen the same per-genome k-mer set, so they detect the same genomes; they differ in what that costs.
    - `loosen` (default, the previous behaviour): widen the single pooled screen index to the densest rate any genome needed. One 8 kb genome drags the whole database's effective screen rate toward the dense `-c`.
    - `band`: keep the pooled index at exactly `--screen-c` and put the extra keys of small genomes in an appended value band. The band block is appended after the existing index, so an older sylph reads such a database as a plain nominal-rate one instead of failing.
    - `none`: no floor at all, i.e. what a build without small-genome handling writes.
- `profile`/`query --screen-small-genomes N` screens genomes with fewer than N stored stage-1 k-mers directly, reading their smallest N dense k-mers back when the database is opened. This recovers small-genome sensitivity on a `--small-genome-screen none` database (including one built by an older sylph) at a per-sample cost proportional to the number of such genomes.
- `sylph inspect` reports a `.syl2db`'s `effective_screen_c` (the densest rate anything in it is screened at) alongside `screen_c`, plus the band-1 key/genome counts.
- Per-genome conversion diagnostics are now aggregated into one counted warning with examples, instead of one line per genome: a viral database densifies millions of genomes, and the old output buried everything else.
- Fixed: `profile --debug` without `-u/--estimate-unknown` aborted on the first sample (it unwrapped an estimated-identity value that only exists with `-u`).
- The hidden debug flag `--screen-dump` now identifies each stage-1 survivor by genome id and contig name, not only by source file (a database built from multi-record FASTA has one file for millions of genomes), and writes its rows in genome order so two runs' dumps can be compared directly.
- Benchmark harness for the above under `benchmarks/mixed_screen_c` (snakemake; MetaVR + GTDB, with a smoke test and a full-corpus configuration).

# sylph v1.0.0: Major performance improvements with new optional database formats - 09-10-2026

- Added a new database format (courtesy of Ben Woodcroft's [weebill](https://github.com/wwood/weebill) fork): `.syl2db`. This uses a two-stage profiling approach: first stage with really sparse k-mers, and second stage with dense k-mers (default `-c 200` as before). This improves speed and memory usage by several times for *prokaryotic genomes* at little to no sensitive cost. 
    - This new `syl2db` format requires sylph v1.0.0; older versions can not use it
    - This new format can be obtained via `sylph convert-db-two-screen`. 
- Made read sketching multi-threaded (thanks to Ben Woodcroft for showing preliminary results)
- Changed CLI help interface. 
- Added the `-d` option to specify databases for `sylph profile`. 
- Default sensitivity increased for small genomes. Min number of k-mers for a genome (`-M`) is set to 10 now, but 7 k-mer hits are required at minimum

# sylph v0.9.0: 10-13-2025

- Added an option `--estimate-read-count` to VERY ROUGHLY output estimated read counts in the "Sequence_abundance" column instead of an actual sequence abundance. This forces `-u`. Only works for short reads right now. 
- Record ids could previously contain tabs, causing sylph's tsv output to not be a true tsv. This is fixed (thanks Donovan Parks). 

## sylph v0.8.1: Small fixes

### Major

- There was a bug in the defn of the "e" constant; fixed now (thanks to Benjamin Lieser). Results may change by about ~2%. 

### Minor
- Allowing `-r` option for `profile` for single-ended reads (thanks again to Florian Plaza Onate for this) 

## sylph v0.8.0: default output behaviour change + more efficient `inspect`.

### Major
* `inspect` option takes much less memory (thanks to Martin Larralde/@althonos for this)
* BREAKING: Changed default behaviour of sylph to write a TSV header even when no genomes are detected (thanks to Florian Plaza Onate for this suggestion)


## sylph v0.7.0: `inspect` option - 2024-11-06

### Major
* Added `inspect` option for inspecting `.syldb/sylsp`.
* Removed native compile flag. 

## sylph v0.6.1: improved automatic detection of sequencing error for -u option

### Major

* Improved automatic estimation of sequencing error for estimating unknown abundances/coverages.

Explanation:

The -u option estimates the % of sequences that are "unknown" i.e. not captured by the database and the "true" coverage. This requires knowledge of sequencing error. Previous versions failed when the sample was too diverse compared to sequencing depth (e.g. low-throughput sequencing or complex (ocean/soil) metagenomes). 

**New fallback added**: For short-reads only, if the diversity is too high relative to sequencing depth, (avg k-mer depth < 3) then 99.5% is used as a fallback sequence identity estimate. 

## sylph v0.6.0 release: New output column, lazy raw paired fastq profiling: 2024-04-06 

### Major

* A new column called `kmers_reassigned` is now in the profile output. This states how many k-mers are lost due to reassignment for that particular genome. 
* `-1, -2` options are now available for `sylph profile`. You can now do `sylph profile database.syldb -1 1.fq -2 2.fq ...`

## sylph v0.5.1 release: **Memory improvement and bug fixes** : Dec 27 2023

### Major

* Scalable cuckoo filters are now used for read deduplication for memory savings. 
* Deduplication algorithm improved. v0.5.0 worked poorly on highly (>15%) duplicated read sets. 
* Shorter reads can be sketched now. Down to 32bp instead of 63 bp before.

## sylph v0.5.0 release: **Big improvements on real illumina data** : Dec 23 2023

### Major

**In previous versions, sylph was underperforming on real illumina data sets**. See https://github.com/bluenote-1577/sylph/issues/5 

This is because many real illumina datasets have a non-trivial number of duplicate reads. Duplicate reads mess up sylph's statistical model.

For the single and paired sketching options, a new deduplication routine has been added. This will be described in version 2 of our preprint. 

**This increases sketching memory by 3-4x but greatly increases performance on real datasets with > 1-2% of duplication, especially for low-abundance genomes**. 

For paired-end illumina reads with non-trivial (> 1% duplication), sylph can now 

1. detect up to many more species low-abundance species below 0.3x coverage
2. give better coverage/abundance estimates for low-abundance species 

### BREAKING

- sequence sketches (sylsp) have changed formats. Sequences will need to be re-sketched.
- `--read-length` option removed and incorporated into the sketches by default. (suggested by @fplaza)

### Other changes

- New warning when `-o` specified and only reads are sketched (https://github.com/bluenote-1577/sylph/issues/7)
- You can now rename sylph samples by specifing a sample naming file with `--sample-names` or `--lS` (suggested by @jolespin)
- Newline delimited files are available in `profile` and `query` now (suggested by @jolespin)


## sylph v0.4.1 release: getting ready for preprinting

### Minor

- small changes for help text, options, and output texts. 

## sylph v0.4.0 release: major interface changes

### BREAKING

- renamed `sylph contain` to `sylph query`. 
- methods for sketching are drastically different now. E.g. we use `-g genome1.fa genom2.fa` for specifying genomes and `-r read1.fa read2.fq` for specifying reads when sketching. 

### Major

- `-u` or `--estimate-unknown` options are now present for estimating unknown organisms in the sample. 
- When using `-u`, associated options `--read-seq-id` and `--read-len` are available for calculating true coverages with sylph, i.e., coverages concordant with read mapping

### Minor

- Coverage calculation is slightly different now.

## sylph v0.3.0 release: first class support for pseudotax, now called "profile" - 2023-10-01

Continuing development of sylph taxonomic profiling. 

### BREAKING

- `--pseudotax` option in previous version is now a new command called `profile`.
- Databases are enabled for profiling by default. 
- Changed file suffices to `syldb` and `sylsp`.

### Major
- Default parameter changes. --min-spacing is set to 30 now. 
- Made profiling faster with some algorithmic tweaks. 
- Coverage calculated slightly differently
- Many small software changes with respect to threading and outputs

## sylph v0.2.0 release: pseudotax improved - 2023-09-19

### BREAKING
- Sylph's *.sylqueries are no longer compatible with older versions of sylph (< v0.2). Files will need to be resketched. 

### Major
- Fixed a major bug for the `--pseudotax` option that required redesigning file formats. Please use `--enable-pseudotax` when using using `contain --pseudotax` from now on.
- `--pseudotax` option gives relative abundances now. We are gaining some confidence that this approach gives a rough, but surprisingly decent taxonomic classification.  
- Changed how `Eff_cov` is calculated. We just use the median coverage now, except when we apply coverage-adjustment 

### Minor
- Fixed command line ambiguity for sketching outputs. `-s` has been replaced with `-d` for `sylph sketch`.
- Sylph outputs the results after processing every sample, instead of batching results, now


## sylph v0.1.0 release - 2023-09-03

### Major

- Added `--pseudotax` option, similar to the `-w` option in mash screen, where k-mers are assigned to the highest ANI genome so redundancy is removed. The output is a very rough taxonomic classification of the sample. 

### Minor

- Some fixes and parameter changes from the v0.0.x releases. 
