# sylph - fast and precise species-level metagenomic + taxonomic profiling for shotgun sequencing 

<p align="center"><img src="assets/jay_logo.png" width = 450 /></p>
 
> [!IMPORTANT]
> Documentation for sylph has moved to https://sylph-docs.github.io/. All GitHub documentation (e.g., Wikis) are out of date. 

**Sylph** is a tool for ultrafast taxonomic profiling of metagenomic shotgun reads. Sylph can detect species and their abundances against massive databases such as GTDB (~200k genomes) in seconds. 

### Why sylph?

1. **Ultrafast, multithreaded, multi-sample**: sylph can be > 50x faster than other methods. Sylph only takes ~5 GB of RAM for profiling against the entire GTDB-R232 database (200k genomes).

2. **Precise species-level profiling**: sylph has less false positives than Kraken and is about as precise and sensitive as marker gene methods (MetaPhlAn, mOTUs). 

3. **Accurate (containment) ANI information**: sylph can give accurate **ANI estimates** between reference genomes and your metagenome sample down to 0.1x coverage.

4. **Customizable databases and pre-built databases**: We offer pre-built databases of [prokaryotes, viruses, eukaryotes](pre‐built-databases.md). Custom databases (e.g. using your own MAGs) are easy to build.  

5. **Short or long reads**: Sylph works for nanopore, PacBio, and short reads. Sylph is the most accurate method [on Oxford Nanopore's independent benchmarks](https://nanoporetech.com/resource-centre/genomic-and-epigenomic-insights-into-microbial-biology-with-nanopore-metagenomic-and-isolate-sequencing).

Documentation, installation, and usage information is available at https://sylph-docs.github.io/.

## Citing sylph

Jim Shaw and Yun William Yu. Rapid species-level metagenome profiling and containment estimation with sylph (2024). Nature Biotechnology.

