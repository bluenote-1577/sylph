# sylph - fast and precise species-level metagenomic profiling with ANIs 

## Introduction

**sylph** is a program that can perform ultrafast (1) **ANI querying** or (2) **metagenomic profiling** for metagenomic shotgun samples. 

**Containment ANI querying**: sylph can search a genome, e.g. E. coli, against your sample. If sylph gives an estimate of 97% ANI, then a genome is contained in your sample with 97% ANI to the queried E. coli genome. 

**Metagenomic profiling**: Like e.g. Kraken or MetaPhlAn, sylph can determine what species are in your sample and their abundances, as well as their _containment ANI in your metagenome_.

<p align="center"><img src="assets/sylph.gif?raw=true"/></p>
<p align="center">
   <i>
   Profiling 1 Gbp of mouse gut reads against 85,205 genomes in a few seconds 
   </i>
</p>

### Why sylph?

1. **Accurate (containment) ANIs down to 0.1x effective coverage**: for bacterial ANI queries of > 90% ANI, sylph can often give accurate ANI estimates down to 0.1x coverage.

2. **Precise species-level profiling**: Our tests show that sylph is more precise than Kraken and about as precise and sensitive as marker gene methods (MetaPhlAn, mOTUs). 

3. **Ultrafast, multithreaded, multi-sample**: sylph can be > 100x faster than MetaPhlAn for multi-sample processing. sylph only takes 13GB of RAM for profiling against the entire GTDB-R214 database (85k genomes).

4. **Easily customized databases**: sylph does not require taxonomic information, so you can profile against [metagenome-assembled genomes (MAGs), viruses, eukaryotes](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases), even assembled contigs, etc. Taxonomic information can be incorporated downstream for traditional profiling reports. 

### How does sylph work?

sylph uses a k-mer containment method, similar to sourmash or Mash. sylph's novelty lies in **using a statistical technique to correct ANI for low coverage genomes** within the sample, allowing accurate ANI queries for even low abundance genomes.

## Changelog

### Version v0.5.0 - Dec 23, 2023. Major breaking update.

#### IMPORTANT

* Big sensitivity boost for real illumina profiling in v0.5 versus v0.4.
* Breaking change: *.sylsp files are now a new format. Old sketches will no longer work.

See the [CHANGELOG](https://github.com/bluenote-1577/sylph/blob/main/CHANGELOG.md) for complete details.


##  Install (current version v0.5.0)

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

Note: the binary is compiled with a different set of libraries (musl instead of glibc), probably impacting performance. 

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

## [Pre-built databases](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases)

The pre-built databases [available here](https://github.com/bluenote-1577/sylph/wiki/Pre%E2%80%90built-databases) can be downloaded and used with sylph for profiling and containment querying. 

## Tutorials and manuals

### [Cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook)

For common use cases and fast explanations, see the above [cookbook](https://github.com/bluenote-1577/sylph/wiki/sylph-cookbook). 

### Tutorials

1. #### [Introduction: 5-minute sylph tutorial outlining basic usage](https://github.com/bluenote-1577/sylph/wiki/5%E2%80%90minute-sylph-tutorial)
2. #### [Taxonomic profiling against GTDB-R214 database with MetaPhlAn output format](https://github.com/bluenote-1577/sylph/wiki/Taxonomic-profiling-with-the-GTDB%E2%80%90R214-database)

### Manuals
1. #### [Sylph's TSV output and containment ANI explanation](https://github.com/bluenote-1577/sylph/wiki/Output-format)
2. #### [Incoporating custom taxonomies to get CAMI-like or MetaPhlAn-like outputs](https://github.com/bluenote-1577/sylph/wiki/Integrating-taxonomic-information-with-sylph)

## Citing sylph

Jim Shaw and Yun William Yu. Metagenome profiling and containment estimation through abundance-corrected k-mer sketching with sylph (2023). bioRxiv.

