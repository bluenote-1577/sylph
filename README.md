# sylph -  ultrafast genome querying against shotgun metagenomes

## Introduction

**sylph** is a program that quickly estimates the average nucleotide identity (ANI) of a collection of genomes against their nearest neighbour genome in your metagenomic shotgun sample.

That is, if we query an E. coli genome and sylph gives an estimate of 97% ANI, then there is a genome in your sample with approximately 97% ANI compared to the queried E. coli genome. 

sylph uses a k-mer containment method, similar to sourmash. However, **sylph uses a statistical adjustment to correct for low coverage genomes** within the sample, allowing accurate ANI queries for even low abundance genomes or shallow depth samples. sylph offers

1. **Accurate ANI queries for genomes of down to 0.1x coverage**: for bacterial species-level ANI queries (> 95%), sylph can give accurate ANI estimates down to 0.1x coverage.

2. **Ultrafast, multithreaded runtimes**: speed is on the scale of Mash or sourmash, but indexing is faster and querying is multithreaded. Entire databases of > 100,000 genomes can be queried against even high-depth samples shotgun samples in **seconds**.

##  Install

#### Option 1: Build from source

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

#### Option 2: Pre-built x86-64 linux statically compiled executable

We offer a pre-built statically compiled executable for x86-64 linux systems. That is, if you're on a x86-64 linux system, you can just download the binary and run it without installing anything. 

```sh
wget https://github.com/bluenote-1577/sylph/releases/download/latest/sylph
chmod +x sylph
./sylph -h
```

Note: the binary is compiled with a different set of libraries (musl instead of glibc), possibly impacting performance (slightly).

#### Option 3: conda install

Forthcoming TODO. 

## Quick start

```sh
# fastq files are _assumed_ to be samples (reads)
# fasta files are _assumed_ to be queries (genomes)
sylph sketch reads1.fq reads2.fq.gz genome1.fa genome2.fa.gz
sylph contain *.sylqueries *.sylsample -t (threads) > all-to-all.tsv

# lazy contain without pre-sketching (convenient, not recommended for large files)
sylph contain reads.fq genome.fa

# paired end reads mode
sylph sketch -1 pair_1.fq -2 pair_2.fq -s sample_prefix

# sketch file with fasta files line-by-line
sylph sketch -l database_file_names.txt -o genomes_output_prefix

# fasta reads
sylph sketch --sample-force fasta_reads.fa

# query contigs instead of genomes
sylph sketch my_contigs.fa -i -o contig_queries
```

## 1 minute tutorial

After installation, clone this repository if you have not done so and run the following.

```sh 
git clone https://github.com/bluenote-1577/sylph
cd sylph

# install sylph. use the below command is rust is installed; see above otherwise
cargo install --path . --root ~/.cargo

# sketch reads and genomes. fastq -> samples, fasta -> queries
sylph sketch test_files/o157_reads.fastq test_files/e.coli*.fa

# sketch all files in line-by-line text file
sylph sketch -l database_fasta_files.txt -o sketch_database_prefix

# query for ANI
sylph contain o157_reads.fastq.sylsample sylph_queries.sylqueries

# ALTERNATIVE: lazy containment without pre-sketching also works 
sylph contain test_files/*
```

You'll see something like the following

```sh
Sample_file     Query_file      Adjusted_ANI    Naive_ANI       ANI_5-95_percentile     Eff_cov Eff_lambda      Lambda_5-95_percentile  Median_cov      Mean_cov_geq1   Containment_ind Contig_name
test_files/o157_reads.fastq     test_files/e.coli-o157.fasta    99.64   96.02   99.51-99.85     0.374   0.374   0.35-0.39       1       1.188   5845/20554      NZ_CP017438.1 Escherichia coli O157:H7 strain 2159 chromosome, complete genome
test_files/o157_reads.fastq     test_files/e.coli-EC590.fasta   98.33   94.38   98.08-98.60     0.321   0.321   0.29-0.35       1       1.168   2987/17970      NZ_CP016182.2 Escherichia coli strain EC590 chromosome, complete genome
test_files/o157_reads.fastq     test_files/e.coli-K12.fasta     98.15   94.30   97.93-98.38     0.334   0.334   0.31-0.35       1       1.171   2938/18124      NC_007779.1 Escherichia coli str. K-12 substr. W3110, complete sequence
```

The o157_reads.fastq is a "sample" containg only E. coli O157 with 1x coverage, 95% identity reads. We queried the E. coli fasta files against this "sample" and got the above results.

1. The 3rd column gives the coverage adjusted ANI between the genome and this sample. **This is the main number you should care about**.
2. The 4th column is the Naive ANI -- what you would approximately get without sylph's statistical adjustment (i.e. if you used Mash or Sourmash).

Notice the big difference between 1. and 2. This is because the reads are only 1x coverage. **MinHash k-mer methods like mash screen and sourmash give biased ANI when coverage is low**. 

However, the Eff_cov gives smaller than 1x: this is because Eff_cov takes into account sequencing error. Sequencing error reduces the k-mer based coverages (errors can cause mutate k-mers). 

Here are the ANIs computed by [skani](https://github.com/bluenote-1577/skani) between the three genomes:

```sh
test_files/e.coli-EC590.fasta	100.00	99.39	98.14
test_files/e.coli-K12.fasta	99.39	100.00	98.09
**test_files/e.coli-o157.fasta	98.14	98.09	100.00**
```

So the ANIs should be 98.14, 98.09, and 100.0 for EC590, K12, and O157 respectively against the sample. As you can see, sylphs estimates are quite good and much more reasonable than the Naive ANI.  

## Output format

```sh
Sample_file	Query_file	Adjusted_ANI	Naive_ANI	ANI_5-95_percentile	Eff_cov	Eff_lambda	Lambda_5-95_percentile	Median_cov	Mean_cov_geq1	Containment_ind Contig_name
test_files/o157_reads.fastq	test_files/e.coli-o157.fasta	99.64	96.02	99.51-99.85	0.374	0.374	0.35-0.39	1	1.188	5845/20554	NZ_CP017438.1 Escherichia coli O157:H7 strain 2159 chromosome, complete genome
```

- Sample_file: the filename of the sample.
- Query_file: the filename of the query.
- Adjusted_ANI: nearest neighbour ANI with coverage adjustment if coverage is low.
- Naive_ANI: nearest neighbour ANI without coverage adjustment.
- ANI_5-95_percentile: [5%,95%] confidence intervals. **Not always a decimal number**.
   * If coverage adjustment is possible: `float-float` e.g. `98.52-99.55`
   * If coverage is too low/high: `NA-NA` is given. 
- Eff_cov: estimate of the coverage. **Always a decimal number.** 
    * If coverage adjustment is possible: this is Eff_lambda
    * If coverage is too low/high: this is Mean_cov_geq1
- Eff_lambda: estimate of the effective coverage parameter. **Not always a decimal number**. 
    * If coverage adjustment is possible: lambda estimate is given
    * If coverage is too low/high: `LOW` or `HIGH` is output
- Lambda_5-95_percentile: [5%, 95%] confidence intervals for lambda. Same format rules as ANI_5-95_percentile.
- Median cov: median k-mer multiplicity for k-mers with >= 1 multiplicity.
- Mean_cov_geq1: mean k-mer multiplicity for k-mers with >= 1 multiplicity.
- Containment_ind: `int/int` showing the containment index, e.g. `959/1053`.
- Contig_name: name of the first contig in the fasta or just the contig name for contig queries.

## Citing sylph

Jim Shaw and Yun William Yu. Ultrafast genome similarity queries from low-abundance metagenomic samples with sylph (Preprint to be released soon). 
