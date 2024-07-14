# sylph - fast and precise species-level metagenomic profiling with ANIs 

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

1. **Precise species-level profiling**: Our tests show that sylph is more precise than Kraken and about as precise and sensitive as marker gene methods (MetaPhlAn, mOTUs). 

2. **Ultrafast, multithreaded, multi-sample**: sylph can be > 50x faster than MetaPhlAn for multi-sample processing. sylph only takes ~15GB of RAM for profiling against the entire GTDB-R220 database (110k genomes).

3. **Accurate (containment) ANIs down to 0.1x effective coverage**: for bacterial ANI queries of > 90% ANI, sylph can often give accurate ANI estimates down to 0.1x coverage.

4. **Easily customized databases**: sylph can profile against [metagenome-assembled genomes (MAGs), viruses, eukaryotes](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases), and more. Taxonomic information can be incorporated downstream for traditional profiling reports. 

### How does sylph work?

sylph uses a k-mer containment method, similar to sourmash or Mash. sylph's novelty lies in **using a statistical technique to correct ANI for low coverage genomes** within the sample, giving accurate results for low abundance genomes. See [here for more information on what sylph can and can not do](https://github.com/bluenote-1577/sylph/wiki/Introduction:-what-is-sylph-and-how-does-it-work%3F). 

## Very quick start

#### Profile metagenome sample against [GTDB-R220](https://gtdb.ecogenomic.org/) (113,104 bacterial/archaeal species representative genomes) 

```sh
# download GTDB-R220 pre-built database (~13 GB)
wget https://storage.googleapis.com/sylph-stuff/gtdb-r220-c200-dbv1.syldb

# multi-sample paired-end profiling (sylph version >= 0.6)
sylph profile gtdb-r220-c200-dbv1.syldb -1 *_1.fastq.gz -2 *_2.fastq.gz -t (threads) > profiling.tsv

# multi-sample single-end profiling
sylph profile gtdb-r220-c200-dbv1.syldb *.fastq -t (threads) > profiling.tsv
```

##  Install (current version v0.6.1)

#### Option 1: conda install 
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/version.svg)](https://anaconda.org/bioconda/sylph)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sylph/badges/latest_release_date.svg)](https://anaconda.org/bioconda/sylph)

```sh
conda install -c bioconda sylph
```

> [!WARNING]
> conda install may break if AVX2 instructions are not available on your CPU. See the [issue here](https://github.com/bluenote-1577/sylph/issues/2). The binary and source install still work. 

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

Note: the binary is compiled with a different set of libraries (musl instead of glibc), probably impacting performance. 

## Standard usage

#### Sketching reads/genomes (indexing)

```sh
# all fasta -> one *.syldb; fasta are assumed to be genomes
sylph sketch genome1.fa genome2.fa -o database
#EQUIVALENT: sylph sketch -g genome1.fa genome2.fa -o database

# multi-sample sketching of paired reads
sylph sketch -1 A_1.fq B_1.fq -2 A_2.fq B_2.fq -d read_sketch_folder

# multi-sample sketching for single end reads, fastq are assumed to be reads
sylph sketch reads.fq 
#EQUIVALENT: sylph sketch -r reads.fq
```

#### Profiling or querying with sketch files
```sh
# ANI querying 
sylph query database.syldb read_sketch_folder/*.sylsp -t (threads) > ani_queries.tsv

# taxonomic profiling 
sylph profile database.syldb read_sketch_folder/*.sylsp -t (threads) > profiling.tsv
```

## Tutorials, manuals, and pre-built databases

### [Pre-built databases](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases)

The pre-built databases [available here](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases) can be downloaded and used with sylph for profiling and containment querying. 

### [Cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook)

For common use cases and fast explanations, see the above [cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook).

### Tutorials
1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://github.com/bluenote-1577/sylph/wiki/5%E2%80%90minute-sylph-tutorial)
2. #### [Taxonomic profiling against GTDB database with MetaPhlAn output format](https://github.com/bluenote-1577/sylph/wiki/Taxonomic-profiling-with-the-GTDB%E2%80%90R214-database)

### Manuals
1. #### [Output format (TSV) and containment ANI explanation](https://github.com/bluenote-1577/sylph/wiki/Output-format)
2. #### [Incoporating custom taxonomies to get CAMI-like or MetaPhlAn-like outputs](https://github.com/bluenote-1577/sylph/wiki/Integrating-taxonomic-information-with-sylph)

### [sylph-utils](https://github.com/bluenote-1577/sylph-utils) 

For incorporating taxonomy and manipulating output formats, see the [sylph-utils repository](https://github.com/bluenote-1577/sylph-utils).

## Changelog

#### Version v0.6.1 - 2024-04-29. 

* Made unknown estimation (-u) more robust for low-depth short-read sequencing. 

See the [CHANGELOG](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md) for complete details.

## Citing sylph

Jim Shaw and Yun William Yu. Metagenome profiling and containment estimation through abundance-corrected k-mer sketching with sylph (2023). bioRxiv.

