# sylph -  ANI genome searching and taxonomic profiling for shotgun metagenomes 

## Introduction

**sylph** is a program that can perform ultrafast **nearest neighbour ANI search** or **taxonomic profiling** for metagenomic shotgun samples. sylph can search tens of thousands of genomes against gigabases of reads in seconds.

**Nearest neighbour containment search (sylph contain)**: sylph can search a genome, e.g. E. coli, against your sample. If sylph gives an estimate of 97% ANI, then a genome is contained in your sample with 97% ANI to the queried E. coli genome. 

**ANI-based taxonomic profiling (sylph profile)**: sylph can determine what species are in your sample as well as their abundances (like Kraken or MetaPhlAn), while also **giving ANI estimates for the classifications**. 

### Why sylph?

1. **Accurate ANI-based genome search down to 0.1x coverage**: for bacterial species-level ANI queries (> 95%), sylph can give accurate ANI estimates down to 0.1x coverage and often even lower.

2. **Precise ANI-based taxonomic profiling**: our preliminary results show sylph is as precise and sensitive as MetaPhlAn4 with better abundance estimates. Database choice is flexible and even viruses/eukaryotes can be profiled.  

3. **Ultrafast, multithreaded runtimes**: sylph is **50x faster than MetaPhlAn** and **10x faster than Kraken**. sylph only takes 10GB of RAM for classifying against the entire GTDB-R214 database (85k genomes). 

sylph uses a k-mer containment method, similar to sourmash or Mash. sylph's novelty lies in **using a statistical technique to correct ANI for low coverage genomes** within the sample, allowing accurate ANI queries for even low abundance genomes or shallow depth samples.

## WARNING EARLY DEVELOPMENT

sylph is being developed rapidly. It has not been officially released yet. I am planning on releasing sylph officially in the next 1-3 months (October-December 2023).  

I have confidence in sylph's results right now, and I believe it works quite well. But be aware that I will have no qualms about making breaking changes until the official release.

That is:
   - any sketches you use may not work by the next release.
   - the command line will change.
   - Parameters will change. 

##  Install (current version v0.3.0)

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
# one fastq -> one *.sylsp; fastq are assumed to be samples (reads)
# all fasta -> one *.syldb; fasta are assumed to be genomes
sylph sketch reads1.fq reads2.fq.gz genome1.fa genome2.fa.gz -o database

# nearest neighbour containment search 
sylph contain database.syldb *.sylsample -t (threads) > containment.tsv

# taxonomic profiling 
sylph profile database.syldb *.sylsample -t (threads) > profiling.tsv
```

Below shows different ways of sketching/indexing samples (reads) or constructing databases

```sh
# lazy inputs without pre-sketching (convenient, not recommended for large files)
sylph profile reads.fq genome1.fa genome2.fa

# paired end reads mode, multi-sample allowed
sylph sketch -1 pairA_1.fq pairB_1.fq -2 pairA_2.fq pairB_2.fq -d read_sketches

# sketch file with fasta files line-by-line
sylph sketch -l database_file_names.txt -o my_database

# fasta reads need to be forced to be samples
sylph sketch --sample-force fasta_reads.fa

# database of contigs instead of genomes
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
Sample_file     Genome_file      Taxonomic_abundance     Sequence_abundance      Adjusted_ANI    Eff_cov ANI_5-95_percentile     Eff_lambda      Lambda_5-95_percentile  Median_cov      Mean_cov_geq1   Containment_ind Naive_ANI       Contig_name
parks_bench_data/ani95_cLOW_stFalse_r8_R1.fq.gz release89/bacteria/RS_GCF_000178875.2_genomic.fna.gz    78.1242 81.8234 97.53   264.000 NA-NA   HIGH    NA-NA   264     264.143 10281/22299     97.53   NC_016901.1 Shewanella baltica OS678, complete genome
```

- Sample_file: the filename of the sample.
- Genome_file: the filename of the genome.
- (*Not present for `contain`*) Taxonomic_abundance: normalized taxnomic abundance as a percentage (i.e. not scaled by sequence length, same as MetaPhlAn)
- (*Not present for `contain`*) Sequence_abundance: normalized sequence abundance as a percentage (i.e. scaled by sequence length, same as Kraken)
- Adjusted_ANI: nearest neighbour ANI. **Most important value**.
    * If coverage adjustment is possible (cov is < 3x cov): returns coverage-adjusted ANI
    * If coverage is too low/high: returns Naive_ANI (see below)
- Eff_cov: estimate of the coverage. *This is an underestimate of the true coverage*. **Always a decimal number.** 
    * If coverage adjustment is possible: this is Eff_lambda
    * If coverage is too low/high: this is `Median_cov`
- ANI_5-95_percentile: [5%,95%] confidence intervals. **Not always a decimal number**.
   * If coverage adjustment is possible: `float-float` e.g. `98.52-99.55`
   * If coverage is too low/high: `NA-NA` is given. 
- Eff_lambda: estimate of the effective coverage parameter. **Not always a decimal number**. 
    * If coverage adjustment is possible: lambda estimate is given
    * If coverage is too low/high: `LOW` or `HIGH` is output
- Lambda_5-95_percentile: [5%, 95%] confidence intervals for lambda. Same format rules as ANI_5-95_percentile.
- Median_cov: median k-mer multiplicity for k-mers with >= 1 multiplicity.
- Mean_cov_geq1: mean k-mer multiplicity for k-mers with >= 1 multiplicity.
- Containment_ind: `int/int` showing the containment index (number of k-mers contained divided by total k-mers), e.g. `959/1053`.
- Naive_ANI: nearest neighbour ANI without coverage adjustment.
- Contig_name: name of the first contig in the fasta or just the contig name for the -i option.

<a name="pre-databases"></a>
## Pre-sketched databases

We have some pre-sketched databases available for download below. 

### Pre-sketched GTDB r214 database (85,202 genomes)

1. `-c 200`, more sensitive database (10 GB): https://storage.googleapis.com/sylph-stuff/v0.3-c200-gtdb-r214.syldb
3. `-c 1000` more efficient, less sensitive database (2 GB): https://storage.googleapis.com/sylph-stuff/v0.2-c100-gtdb-r214.syldb

Quick usage example

```sh
# faster, less sensitive database
wget https://storage.googleapis.com/sylph-stuff/v0.3-tax-c1000-gtdb-r214.syldb
sylph profile reads.fq v0.3-c200-gtdb-r214.syldb -t 30 > results.tsv
```

## Citing sylph

Jim Shaw and Yun William Yu. Ultrafast, coverage-corrected genome similarity queries for metagenomic shotgun samples with sylph (Preprint to be released soon). 

