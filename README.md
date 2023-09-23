# sylph -  ultrafast, coverage-aware genome querying against shotgun metagenomes 

## Introduction

**sylph** is a program that quickly estimates the average nucleotide identity (ANI) of genomes against their nearest neighbour genome in your metagenomic shotgun sample. An experimental taxonomic classifier option is available. 

That is, if we query an E. coli genome and sylph gives an estimate of 97% ANI, then there is a genome in your sample with approximately 97% ANI compared to the queried E. coli genome. 

sylph uses a k-mer containment method. sylph's novelty lies in **using a statistical technique to correct ANI for low coverage genomes** within the sample, allowing accurate ANI queries for even low abundance genomes or shallow depth samples. sylph offers:

1. **Accurate ANI queries for genomes of down to 0.1x coverage**: for bacterial species-level ANI queries (> 95%), sylph can give accurate ANI estimates down to 0.1x coverage.

2. **Ultrafast, multithreaded runtimes**: speed is on the scale of Mash or sourmash, but indexing is faster and querying is multithreaded. Entire databases of > 100,000 genomes can be queried against even high-depth samples shotgun samples in **seconds**.
   
3. **(EXPERIMENTAL) Fast taxonomic classification**: sylph can be turned into a true taxonomic classifier (with relative abundances) with the `--pseudotax` option. This is experimental, but seems to work quite well from brief testing.

## WARNING EARLY DEVELOPMENT

sylph is being developed rapidly. It has not been officially released yet. I am planning on releasing sylph officially in the next 1-3 months (October-December 2023).  

I have confidence in sylph's results right now, and I believe it works quite well. But be aware that 

- I will have no qualms about making breaking changes until the official release.
- Parameters will change. 

##  Install (current version v0.2.0)

#### Option 1: conda install 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/version.svg)](https://anaconda.org/bioconda/sylph)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/latest_release_date.svg)](https://anaconda.org/bioconda/sylph)

```sh
conda install -c bioconda sylph
```

#### Option 2: Pre-built x86-64 linux statically compiled executable

We offer a pre-built statically compiled executable for x86-64 linux systems. That is, if you're on a x86-64 linux system, you can just download the binary and run it without installing anything. 

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), possibly impacting performance (slightly).

#### Option 3: Build from source

Requirements:
1. [rust](https://www.rust-lang.org/tools/install) (version > 1.63) programming language and associated tools such as cargo are required and assumed to be in PATH.
2. A c compiler (e.g. GCC)
3. make
4. cmake

Building takes a few minutes (depending on # of cores).

```sh
git clone https://github.com/bluenote-1577/sylph
cd sylph

# If default rust install directory is ~/.cargo
cargo install --path . --root ~/.cargo
sylph contain test_files/*

# If ~/.cargo doesn't exist use below commands instead
#cargo build --release
#./target/release/sylph -h
```

## Quick start

```sh
# one fastQ -> one *.sylsample; fastQ  assumed to be samples (reads)
# all fastA -> one *.sylqueries; fastA  assumed to be queries (genomes)
sylph sketch reads1.fq reads2.fq.gz genome1.fa genome2.fa.gz
sylph contain *.sylqueries *.sylsample -t (threads) > all-to-all.tsv

# do pseudotaxonomic classification instead of nearest neighbour ANI
sylph sketch *.fa --enable-pseudotax -o pseudotax-enabled-db
sylph contain query.fq pseudotax-enabled-db.sylqueries --pseudotax 
```

Below shows different ways of sketching/indexing samples (reads) or queries (genomes)

```sh
# lazy contain without pre-sketching (convenient, not recommended for large files)
sylph contain reads.fq genome.fa

# paired end reads mode, multi-sample allowed
sylph sketch -1 pairA_1.fq pairB_1.fq -2 pairA_2.fq pairB_2.fq -s sample_prefix

# sketch file with fasta files line-by-line
sylph sketch -l database_file_names.txt -o genomes_output_prefix

# fasta reads
sylph sketch --sample-force fasta_reads.fa

# query contigs instead of genomes
sylph sketch my_contigs.fa -i -o contig_queries

```

See [Pre-sketched databases](#pre-databases) below to download pre-indexed databases. 

## Tutorials and manuals

### Tutorials

1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://github.com/bluenote-1577/sylph/wiki/5%E2%80%90minute-sylph-tutorial)

### Manuals

Manuals forthcoming...

## Output format

```sh
Sample_file	Query_file	Adjusted_ANI	Naive_ANI	ANI_5-95_percentile	Eff_cov	Eff_lambda	Lambda_5-95_percentile	Median_cov	Mean_cov_geq1	Containment_ind Contig_name
test_files/o157_reads.fastq	test_files/e.coli-o157.fasta	99.64	96.02	99.51-99.85	0.374	0.374	0.35-0.39	1	1.188	5845/20554	NZ_CP017438.1 Escherichia coli O157:H7 strain 2159 chromosome, complete genome
```

- Sample_file: the filename of the sample.
- Query_file: the filename of the query.
- (*Only if `--pseudotax` enabled*) Relative_abundance: relative abundance as a percentage
- Adjusted_ANI: nearest neighbour ANI. **Most important value**.
    * If coverage adjustment is possible (cov is < 3x cov): returns coverage-adjusted ANI
    * If coverage is too low/high: returns Naive_ANI (see below)
- Naive_ANI: nearest neighbour ANI without coverage adjustment.
- ANI_5-95_percentile: [5%,95%] confidence intervals. **Not always a decimal number**.
   * If coverage adjustment is possible: `float-float` e.g. `98.52-99.55`
   * If coverage is too low/high: `NA-NA` is given. 
- Eff_cov: estimate of the coverage. **Always a decimal number.** 
    * If coverage adjustment is possible: this is Eff_lambda
    * If coverage is too low/high: this is `Median_cov`
- Eff_lambda: estimate of the effective coverage parameter. **Not always a decimal number**. 
    * If coverage adjustment is possible: lambda estimate is given
    * If coverage is too low/high: `LOW` or `HIGH` is output
- Lambda_5-95_percentile: [5%, 95%] confidence intervals for lambda. Same format rules as ANI_5-95_percentile.
- Median_cov: median k-mer multiplicity for k-mers with >= 1 multiplicity.
- Mean_cov_geq1: mean k-mer multiplicity for k-mers with >= 1 multiplicity.
- Containment_ind: `int/int` showing the containment index, e.g. `959/1053`.
- Contig_name: name of the first contig in the fasta or just the contig name for contig queries.

<a name="pre-databases"></a>
## Pre-sketched databases

We have some pre-sketched databases available for download below. 

### Pre-sketched GTDB r214 database (85,202 genomes)
The databases with `*-tax-*` in the title indicate `--pseudotax` can be used. The results are otherwise the same. These sketches work only for sylph v0.2 and above. 

1. `-c 100`, more sensitive database (20 GB): https://storage.googleapis.com/sylph-stuff/v0.2-tax-c100-gtdb-r214.sylqueries
2. `-c 100`, more sensitive database without `--pseudotax` capabilities (8 GB): https://storage.googleapis.com/sylph-stuff/v0.2-standard-c100-gtdb-r214.sylqueries
3. `-c 1000` more efficient, less sensitive database (2 GB): https://storage.googleapis.com/sylph-stuff/v0.2-tax-c1000-gtdb-r214.sylqueries

Quick usage example

```sh
# faster, less sensitive database
wget https://storage.googleapis.com/sylph-stuff/v0.2-tax-c1000-gtdb-r214.sylqueries
sylph contain your_sample.sylsample v0.2-tax-c1000-gtdb-r214.sylqueries -t 30 > results.tsv
```

## Citing sylph

Jim Shaw and Yun William Yu. Ultrafast, coverage-corrected genome similarity queries for metagenomic shotgun samples with sylph (Preprint to be released soon). 
