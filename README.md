# sylph - fast and precise species-level metagenomic profiling with ANIs 

> [!IMPORTANT]
> All documentation for sylph is migrating to https://sylph-docs.github.io/.
>
> The GitHub versions of the README/Wikis/Manuals are now OUTDATED. 

## Introduction

**sylph** is a program that performs ultrafast (1) **ANI querying** or (2) **metagenomic profiling** for metagenomic shotgun samples. 

**Containment ANI querying**: sylph can search a genome, e.g. E. coli, against your sample. If sylph outputs an estimate of 97% ANI, your sample contains an E. coli with 97% ANI to the queried genome.

**Metagenomic profiling**: sylph can determine the species/taxa in your sample and their abundances, just like [Kraken](https://ccb.jhu.edu/software/kraken/) or [MetaPhlAn](https://github.com/biobakery/MetaPhlAn).

<p align="center"><img src="assets/sylph.gif?raw=true"/></p>
<p align="center">
   <i>
   Profiling 1 Gbp of mouse gut reads against 85,205 genomes in a few seconds 
   </i>
</p>


### Why sylph?

1. **Precise species-level profiling**: sylph has less false positives than Kraken and is about as precise and sensitive as marker gene methods (MetaPhlAn, mOTUs). 

2. **Ultrafast, multithreaded, multi-sample**: sylph can be > 50x faster than other methods. Sylph only takes ~15GB of RAM for profiling against the entire GTDB-R220 database (110k genomes).

3. **Accurate (containment) ANI information**: sylph can give accurate **ANI estimates** between reference genomes and your metagenome sample down to 0.1x coverage.

4. **Customizable databases and pre-built databases**: We offer pre-built databases of [prokaryotes, viruses, eukaryotes](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases). Custom databases (e.g. using your own MAGs) are easy to build.  

5. **Short or long reads**: Sylph was also the most accurate method [on Oxford Nanopore's independent benchmarks](https://nanoporetech.com/resource-centre/genomic-and-epigenomic-insights-into-microbial-biology-with-nanopore-metagenomic-and-isolate-sequencing).

### How does sylph work?

sylph uses a k-mer containment method. sylph's novelty lies in **using a statistical technique to estimate k-mer containment for low coverage genomes** , giving accurate results for low abundance organisms. See [here for more information on what sylph can and can not do](https://github.com/bluenote-1577/sylph/wiki/Introduction:-what-is-sylph-and-how-does-it-work%3F). 

## Very quick start

#### Profile metagenome sample against [GTDB-R220](https://gtdb.ecogenomic.org/) (113,104 bacterial/archaeal species representative genomes) 

```sh
conda install -c bioconda sylph

# download GTDB-R220 pre-built database (~13 GB)
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb

# multi-sample paired-end profiling (sylph version >= 0.6)
sylph profile gtdb-r220-c200-dbv1.syldb -1 *_1.fastq.gz -2 *_2.fastq.gz -t (threads) > profiling.tsv

# multi-sample single-end profiling
sylph profile gtdb-r220-c200-dbv1.syldb *.fastq -t (threads) > profiling.tsv
```

##  Install 

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
sylph profile test_files/*
```
#### Option 3: Pre-built x86-64 linux statically compiled executable

If you're on an x86-64 system, you can download the binary and use it without any installation. 

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), probably impacting performance. 

## Tutorials, manuals, and pre-built databases

### [Pre-built databases](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases)

The pre-built databases [available here](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases) can be downloaded and used with sylph for profiling and containment querying. 

### [Cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook)

For common use cases and fast explanations, see the above [cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook).

### Tutorials
1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://github.com/bluenote-1577/sylph/wiki/5%E2%80%90minute-sylph-tutorial)
2. #### [Taxonomic profiling against GTDB database with MetaPhlAn-like output format](https://github.com/bluenote-1577/sylph/wiki/Taxonomic-profiling-with-the-GTDB%E2%80%90R214-database)

### Manuals
1. #### [Output format (TSV) and containment ANI explanation](https://github.com/bluenote-1577/sylph/wiki/Output-format)
2. #### [Taxonomic integration and custom taxonomies](https://github.com/bluenote-1577/sylph/wiki/Incorporating-taxonomic-information-into-sylph-with-sylph%E2%80%90tax)

### [sylph-tax](https://github.com/bluenote-1577/sylph-tax) 

To incorporate *taxonomy* into sylph's outputs, see the [sylph-tax repository](https://github.com/bluenote-1577/sylph-tax). 

> [!TIP] 
> The new [sylph-tax](https://github.com/bluenote-1577/sylph-tax) program replaces the old [sylph-utils](https://github.com/bluenote-1577/sylph-utils) repository. 

## Changelog

#### Version v0.8.0 - 2024-12-12. 

* Made the `inspect` option much less memory intensive. Slightly changed outputs when no genomes are found.

See the [CHANGELOG](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md) for complete details.

## Citing sylph

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.

