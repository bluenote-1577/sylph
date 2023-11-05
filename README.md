# sylph -  ANI genome querying and taxonomic profiling for shotgun metagenomes 

## Introduction

**sylph** is a program that can perform ultrafast (1) **ANI querying** or (2) **taxonomic profiling** for metagenomic shotgun samples. 

**ANI querying**: sylph can search a genome, e.g. E. coli, against your sample. If sylph gives an estimate of 97% ANI, then a genome is contained in your sample with 97% ANI to the queried E. coli genome. 

**Taxonomic profiling**: Just like e.g. Kraken or MetaPhlAn, sylph can determine what species are in your sample and their abundances, as well as their _ANI to the database_.

### Why sylph?

1. **Accurate ANIs down to 0.1x coverage**: for bacterial ANI queries of > 90% ANI, sylph can give accurate ANI estimates down to 0.1x coverage and often even lower.

2. **Precise, flexible taxonomic profiling**: Our tests show that sylph is as precise and sensitive as MetaPhlAn4, but with better abundance estimates. Compared to MetaPhlAn4, database choice and read technology (e.g. nanopore) are flexible. Even viruses/eukaryotes can be profiled.  

3. **Ultrafast, multithreaded, multi-sample**: sylph is > 100x faster than MetaPhlAn. sylph only takes 10GB of RAM for profiling against the entire GTDB-R214 database (85k genomes).

4. **Easily customized databases**: sylph does not require taxonomic information, so anything you can profile against custom metagenome-assembled genomes (MAGs), viruses, even assembled contigs, etc.

### How does sylph work?

sylph uses a k-mer containment method, similar to sourmash or Mash. sylph's novelty lies in **using a statistical technique to correct ANI for low coverage genomes** within the sample, allowing accurate ANI queries for even low abundance genomes.

## WARNING EARLY DEVELOPMENT

sylph is being developed rapidly. It has not been officially released yet. I am planning on releasing sylph officially in the next 1-3 months (October-December 2023).  

The following may change:
   - any sketches you use may not work by the next release
   - the command line options
   - Parameters will change 

##  Install (current version v0.4.0)

#### Option 1: conda install 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/version.svg)](https://anaconda.org/bioconda/sylph)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/latest_release_date.svg)](https://anaconda.org/bioconda/sylph)

```sh
conda install -c bioconda sylph
```

#### Option 2: Build from source

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
sylph query test_files/*
```
#### Option 3: Pre-built x86-64 linux statically compiled executable

If you're on an x86-64 system, you can download the binary and use it without any installation. 

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), possibly impacting performance. 

## Quick start

```sh
# all fasta -> one *.syldb; fasta are assumed to be genomes
sylph sketch genome1.fa genome2.fa -o database
#EQUIVALENT: sylph sketch -g genome1.fa genome2.fa -o database

# multi-sample sketching of paired reads
sylph sketch -1 A_1.fq B_1.fq -2 A_2.fq B_2.fq -d output_read_sketch_folder

# multi-sample sketching for single end reads, fastq are assumed to be reads
sylph sketch reads.fq 
#EQUIVALENT: sylph sketch -r reads.fq

# ANI querying 
sylph query database.syldb *.sylsp -t (threads) > ani_queries.tsv

# taxonomic profiling 
sylph profile database.syldb *.sylsp -t (threads) > profiling.tsv
```

See [Pre-sketched databases](#pre-databases) below to download pre-indexed databases. 

## Tutorials and manuals

### [Cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook)

For common use-cases and fast explanations, see the above [cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook). 

### Tutorials

1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://github.com/bluenote-1577/sylph/wiki/5%E2%80%90minute-sylph-tutorial)

### Manuals

1. #### [Incoporating taxonomy to get CAMI-like or MetaPhlAn-like outputs for GTDB (and custom taxonomy)](https://github.com/bluenote-1577/sylph/wiki/MetaPhlAn-or-CAMI%E2%80%90like-output-with-the-GTDB-database)

## Output format

```sh
Sample_file   Genome_file   Taxonomic_abundance   Sequence_abundance   Adjusted_ANI   Eff_cov   ANI_5-95_percentile   Eff_lambda   Lambda_5-95_percentile   Median_cov   Mean_cov_geq1   Containment_ind   Naive_ANI   Contig_name
reads.fq   genome.fa   78.1242   81.8234   97.53   264.000   NA-NA   HIGH   NA-NA   264   264.143   10281/22299   97.53   NC_016901.1 Shewanella baltica OS678, complete genome
```

- Sample_file: the filename of the sample.
- Genome_file: the filename of the genome.
- (*Not present for `query`*) Taxonomic_abundance: normalized taxnomic abundance as a percentage (i.e. not scaled by sequence length; same as MetaPhlAn)
- (*Not present for `query`*) Sequence_abundance: normalized sequence abundance as a percentage (i.e. scaled by sequence length; same as Kraken)
- Adjusted_ANI: nearest neighbour ANI.
    * If coverage adjustment is possible (cov is < 3x cov): returns coverage-adjusted ANI
    * If coverage is too low/high: returns Naive_ANI (see below)
- Eff_cov/True_cov: an estimate of the effective, or if `-u` specified, the true coverage. **Always a decimal number.** 
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
- Contig_name: name of the first contig in the fast (or just the contig name for the -i option).

<a name="pre-databases"></a>
## Pre-sketched databases

We have some pre-sketched databases available for download below. 

### Pre-sketched GTDB r214 database (85,202 genomes). Works with v0.3.0 - current

1. `-c 200`, more sensitive database (10 GB): https://storage.googleapis.com/sylph-stuff/v0.3-c200-gtdb-r214.syldb
3. `-c 1000` more efficient, less sensitive database (2 GB): https://storage.googleapis.com/sylph-stuff/v0.3-c1000-gtdb-r214.syldb

Quick usage example

```sh
# faster, less sensitive database
wget https://storage.googleapis.com/sylph-stuff/v0.3-c1000-gtdb-r214.syldb
sylph profile reads.fq v0.3-c200-gtdb-r214.syldb -t 30 > results.tsv
```

## Citing sylph

Jim Shaw and Yun William Yu. Ultrafast, coverage-corrected genome similarity queries for metagenomic shotgun samples with sylph (Preprint to be released soon). 

