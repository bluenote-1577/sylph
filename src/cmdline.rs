use clap::{Args, Parser, Subcommand};
use crate::constants::*;
use crate::twostage_db::SmallGenomeMode;

#[derive(Parser)]
#[clap(author, version, about = "Ultrafast genome ANI queries and taxonomic profiling for metagenomic shotgun samples.\n\n--- Preparing inputs by sketching (indexing)\n## fastq (reads) and fasta (genomes all at once)\n## *.sylsp found in -d; *.syldb given by -o\nsylph sketch -t 5 sample1.fq sample2.fq genome1.fa genome2.fa -o genome1+genome2 -d sample_dir\n\n## paired-end reads\nsylph sketch -1 a_1.fq b_1.fq -2 b_2.fq b_2.fq -d paired_sketches\n\n--- Taxonomic profiling with relative abundances and ANI\nsylph profile *.syldb *.sylsp > all-to-all-profile.tsv\n\n--- Direct profiling against database with raw reads\nsylph profile *.syldb -1 sampleA_1.fq -2 sampleA_2.fq", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb. 
    #[clap(display_order = 1)]
    Sketch(SketchArgs),
    /// Coverage-adjusted ANI querying between databases and samples without abundances.
    #[clap(display_order = 3)]
    Query(ContainArgs),
    ///Species-level taxonomic profiling with abundances and ANIs. 
    #[clap(display_order = 2)]
    Profile(ContainArgs),
    ///Inspect sketched .syldb and .sylsp files.
    #[clap(arg_required_else_help = true, display_order = 4)]
    Inspect(InspectArgs),
    #[clap(arg_required_else_help = true, display_order = 5)]
    /// Convert a standard database (.syldb) into a two-stage seekable database (.syl2db), automatically used by `query`/`profile` when given as input. Much faster for genomes >~200 kbp with no accuracy change. Don't use for small genomes. 
    ConvertDbTwoScreen(DbConvertArgs),
}

#[derive(Args)]
pub struct DbConvertArgs {
    #[clap(multiple=true, help = "Standard genome database sketches (*.syldb) to convert")]
    pub files: Vec<String>,
    #[clap(short='o', long="output", help = "Output two-stage database name (.syl2db appended)")]
    pub output: String,
    #[clap(long="screen-c", default_value_t = SCREEN_C_DEFAULT, help = "Subsampling rate -c of the small in-memory stage-1 SCREEN index (the bincoded sparse hashes). Must be >= the database -c. A coarser (larger) value gives a smaller/faster screen index. The dense per-genome blocks always keep every k-mer at the database -c.")]
    pub screen_c: usize,
    #[clap(long="min-sparse-kmers", default_value_t = SPARSE_TARGET_MIN_DEFAULT, help = "Minimum stage-1 sparse/screen k-mers per genome; genomes whose nominal --screen-c subsample would fall short use a denser, genome-specific screen rate to reach this floor (or all of their dense k-mers if they have fewer than this to begin with). Must be >= 1.")]
    pub min_sparse_kmers: usize,
    #[clap(long="small-genome-screen", arg_enum, default_value_t = SmallGenomeMode::Loosen,
        help = "How the mixture of screen rates that --min-sparse-kmers creates is stored. `loosen`: widen the single pooled index to the densest rate any genome needed (one small genome slows the screen for the whole database). `band`: keep the pooled index at --screen-c and put the extra keys of small genomes in a separate value band. `none`: no floor at all, so the file is identical to a build without small-genome handling; use `profile --screen-small-genomes` to recover them reader-side.")]
    pub small_genome_screen: SmallGenomeMode,
    #[clap(long="min-contain",default_value_t = 7, help_heading = "ALGORITHM", help = "Throw away genomes with fewer than this many dense k-mers (they could never pass `profile`/`query`'s hit threshold at the matching default anyway, or are likely erroneous/fragmentary genomes). Set to 7 in line with default profiling options. ")]
    pub min_contain: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help_heading = "INPUT", help = "fasta/fastq files; gzip optional. Default: fastq file produces a sample sketch (*.sylsp) while fasta files are combined into a database (*.syldb).")]
    pub files: Vec<String>,
    #[clap(short='o',long="out-name-db", default_value = "database", help_heading = "OUTPUT", help = "Output name for database sketch (with .syldb appended)")]
    pub db_out_name: String,
    #[clap(short='d',long="sample-output-directory", default_value = "./", help_heading = "OUTPUT", help = "Output directory for sample sketches")]
    pub sample_output_dir: String,
    #[clap(short,long="individual-records", help_heading = "GENOME INPUT", help = "Use individual records (contigs) for database construction")]
    pub individual: bool,
    #[clap(multiple=true,short,long="reads", help_heading = "SINGLE-END INPUT", help = "Single-end fasta/fastq reads")]
    pub reads: Option<Vec<String>>,
    #[clap(multiple=true,short='g', long="genomes", help_heading = "GENOME INPUT", help = "Genomes in fasta format")]
    pub genomes: Option<Vec<String>>,
    #[clap(short,long="list", help_heading = "INPUT", help = "Newline delimited file with inputs; fastas -> database, fastq -> sample")]
    pub list_sequence: Option<String>,
    #[clap(long="rl", hidden=true, help_heading = "SINGLE-END INPUT", help = "Newline delimited file; inputs assumed reads")]
    pub list_reads: Option<String>,
    #[clap(long="gl", help_heading = "GENOME INPUT", help = "Newline delimited file; inputs assumed genomes")]
    pub list_genomes: Option<String>,
    #[clap(long="l1", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are first pair of PE reads")]
    pub list_first_pair: Option<String>,
    #[clap(long="l2", help_heading = "PAIRED-END INPUT", help = "Newline delimited file; inputs are second pair of PE reads")]
    pub list_second_pair: Option<String>,
    #[clap(long="lS", help_heading = "INPUT", help = "Newline delimited file; read sketches are renamed to given sample names")]
    pub list_sample_names: Option<String>,
    #[clap(multiple=true, short='S', long="sample-names", help_heading = "INPUT", help = "Read sketches are renamed to given sample names")]
    pub sample_names: Option<Vec<String>>,

    #[clap(short, default_value_t = 31,help_heading = "ALGORITHM", help ="Value of k. Only k = 21, 31 are currently supported")]
    pub k: usize,
    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='s', long="sample-threads", help = "Number of input files to sketch concurrently (out of the shared -t/--threads budget); each concurrently-sketched file's own pipeline gets threads / this value threads. Default: min(number of input files, --threads)")]
    pub sample_threads: Option<usize>,
    #[clap(long="ram-barrier", help = "Stop multi-threaded read sketching when (virtual) RAM is past this value (in GB). Does NOT guarantee max RAM limit", hidden=true)]
    pub max_ram: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,


    #[clap(long="no-dedup", help_heading = "ALGORITHM", help = "Disable read deduplication procedure. Reduces memory; not recommended for illumina data")]
    pub no_dedup: bool,
    #[clap(long="disable-profiling", help_heading = "ALGORITHM", help = "Disable sylph profile usage for databases; may decrease size and make sylph query slightly faster", hidden=true)]
    pub no_pseudotax: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the genomes")]
    pub min_spacing_kmer: usize,
    #[clap(long="fpr", default_value_t = DEFAULT_FPR, help_heading = "ALGORITHM", help = "False positive rate for read deduplicate hashing; valid values in [0,1).")]
    pub fpr: f64,
    #[clap(short='1',long="first-pairs", multiple=true, help_heading = "PAIRED-END INPUT", help = "First pairs for paired end reads")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pairs", multiple=true, help_heading = "PAIRED-END INPUT", help = "Second pairs for paired end reads")]
    pub second_pair: Vec<String>,

    #[clap(long="sketch-batch-size", hidden=true, default_value_t = DEFAULT_SKETCH_BATCH_SIZE,
        help = "Reads per batch handed from the I/O thread to worker threads in multi-threaded read sketching.")]
    pub sketch_batch_size: usize,
    #[clap(long="sketch-batch-max-bytes", hidden=true, default_value_t = DEFAULT_SKETCH_BATCH_MAX_BYTES,
        help = "Maximum accumulated sequence bytes per batch handed from the I/O thread to worker threads (safety cap alongside --sketch-batch-size; matters mainly for long-read data where a fixed record count can otherwise produce very large batches).")]
    pub sketch_batch_max_bytes: usize,
    #[clap(long="sketch-channel-depth", hidden=true, default_value_t = DEFAULT_SKETCH_CHANNEL_DEPTH,
        help = "In-flight batches buffered per worker thread between the I/O thread and worker threads (memory/backpressure knob).")]
    pub sketch_channel_depth: usize,
    #[clap(long="sketch-shards", hidden=true,
        help = "Number of dedup/count shards for multi-threaded read sketching. Default: automatically scaled to worker thread count.")]
    pub sketch_shards: Option<usize>,
    #[clap(long="no-sketch-pipeline", hidden=true,
        help = "Disable the multi-threaded per-file read-sketching pipeline; always use the legacy single-threaded-per-file path.")]
    pub no_sketch_pipeline: bool,
}

#[derive(Args, Clone)]
pub struct ContainArgs {
    #[clap(multiple=true, help = "Pre-sketched *.syldb/*.sylsp files. Raw single-end fastq/fasta are allowed and will be automatically sketched to .sylsp/.syldb")]
    pub files: Vec<String>,

    
    #[clap(short='l',long="list", help = "Newline delimited file of file inputs",help_heading = "INPUT/OUTPUT")]
    pub file_list: Option<String>,

    #[clap(short='d', long="databases", multiple=true, help = "Explicitly specify database files (*.syldb/*.syl2db) instead of/in addition to positional input", help_heading = "INPUT/OUTPUT")]
    pub databases: Vec<String>,

    #[clap(long,default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum k-mer multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(short='M',long,default_value_t = 10., help_heading = "ALGORITHM", help = "Discard genomes with fewer than this many sampled k-mers")]
    pub min_number_kmers: f64,
    #[clap(long="min-contain",default_value_t = 7, help_heading = "ALGORITHM", help = "Minimum number of contained k-mers required for a hit")]
    pub min_contain: usize,
    #[clap(short, long="minimum-ani", help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 90 for query and 95 for profile. Smaller than 95 for profile will give inaccurate results." )]
    pub minimum_ani: Option<f64>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='s', long="sample-threads", help = "Number of samples to be processed concurrently (out of the shared -t/--threads budget); each concurrently-processed sample's own sketching pipeline gets threads / this value threads. Default: min(number of input samples, --threads)")]
    pub sample_threads: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(long="estimate-read-counts", help_heading = "ALGORITHM", help = "Very roughly estimate read counts in the 'Sequence_abundance' column instead of relative abundance. This forces `-u`, which may have caveats for long reads and complex environments.")]
    pub estimate_read_counts: bool,

    #[clap(short='u', long="estimate-unknown", help_heading = "ALGORITHM", help = "Estimate true coverage and scale sequence abundance in `profile` by estimated unknown sequence percentage" )]
    pub estimate_unknown: bool,
    
    #[clap(short='I',long="read-seq-id", help_heading = "ALGORITHM", help = "Sequence identity (%) of reads. Only used in -u option and overrides automatic detection. ")]
    pub seq_id: Option<f64>,

    //#[clap(short='l', long="read-length", help_heading = "ALGORITHM", help = "Read length (single-end length for pairs). Only necessary for short-read coverages when using --estimate-unknown. Not needed for long-reads" )]
    //pub read_length: Option<usize>,

    #[clap(short='R', long="redundancy-threshold", help_heading = "ALGORITHM", help = "Removes redundant genomes up to a rough ANI percentile when profiling", default_value_t = 99.0, hidden=true)]
    pub redundant_ani: f64,

    #[clap(short='r',long="reads", multiple=true, help = "Single-end raw reads (fastx/gzip)", display_order = 1, help_heading = "SKETCHING")]
    pub reads: Vec<String>,

    #[clap(short='1', long="first-pairs", multiple=true, help = "First pairs for raw paired-end reads (fastx/gzip)", help_heading = "SKETCHING")]
    pub first_pair: Vec<String>,

    #[clap(short='2', long="second-pairs", multiple=true, help = "Second pairs for raw paired-end reads (fastx/gzip)", help_heading = "SKETCHING")]
    pub second_pair: Vec<String>,

    #[clap(short, default_value_t = 200, help_heading = "SKETCHING", help = "Subsampling rate. Does nothing for pre-sketched files")]
    pub c: usize,
    #[clap(short, default_value_t = 31, help_heading = "SKETCHING", help = "Value of k. Only k = 21, 31 are currently supported. Does nothing for pre-sketched files")]
    pub k: usize,
    #[clap(short,long="individual-records", help_heading = "SKETCHING", help = "Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-sketched files")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "SKETCHING", help = "Minimum spacing between selected k-mers on the database genomes. Does nothing for pre-sketched files")]
    pub min_spacing_kmer: usize,

    #[clap(long="sketch-batch-size", hidden=true, default_value_t = DEFAULT_SKETCH_BATCH_SIZE, help_heading = "SKETCHING",
        help = "Reads per batch handed from the I/O thread to worker threads in multi-threaded read sketching.")]
    pub sketch_batch_size: usize,
    #[clap(long="sketch-batch-max-bytes", hidden=true, default_value_t = DEFAULT_SKETCH_BATCH_MAX_BYTES, help_heading = "SKETCHING",
        help = "Maximum accumulated sequence bytes per batch handed from the I/O thread to worker threads (safety cap alongside --sketch-batch-size; matters mainly for long-read data where a fixed record count can otherwise produce very large batches).")]
    pub sketch_batch_max_bytes: usize,
    #[clap(long="sketch-channel-depth", hidden=true, default_value_t = DEFAULT_SKETCH_CHANNEL_DEPTH, help_heading = "SKETCHING",
        help = "In-flight batches buffered per worker thread between the I/O thread and worker threads (memory/backpressure knob).")]
    pub sketch_channel_depth: usize,
    #[clap(long="sketch-shards", hidden=true, help_heading = "SKETCHING",
        help = "Number of dedup/count shards for multi-threaded read sketching. Default: automatically scaled to worker thread count.")]
    pub sketch_shards: Option<usize>,
    #[clap(long="no-sketch-pipeline", hidden=true, help_heading = "SKETCHING",
        help = "Disable the multi-threaded per-file read-sketching pipeline; always use the legacy single-threaded-per-file path.")]
    pub no_sketch_pipeline: bool,

    #[clap(short='o',long="output-file", help = "Output to this file (TSV format). [default: stdout]", help_heading="INPUT/OUTPUT")]
    pub out_file_name: Option<String>,
    #[clap(long="log-reassignments", help = "Output information for how k-mers for genomes are reassigned during `profile`. Caution: can be verbose and slows down computation.")]
    pub log_reassignments: bool,

    #[clap(long="screen-ani", default_value_t = SCREEN_MIN_ANI_DEFAULT, help_heading = "TWO-STAGE PROFILING", help = "Two-stage databases (.syl2db) only: minimum adjusted ANI (0-100) for a genome to pass the first-stage screen. Deliberately permissive; the dense stage recovers specificity.")]
    pub screen_ani: f64,
    #[clap(long="screen-dump", hidden=true, help_heading = "TWO-STAGE PROFILING", help = "Debug: write a TSV of every stage-1 screen survivor (genome, matched/total screen k-mers, naive/adjusted ANI, median coverage) to this file.")]
    pub screen_dump: Option<String>,
    #[clap(long="screen-small-genomes", default_value_t = 0, help_heading = "TWO-STAGE PROFILING", help = "Two-stage databases (.syl2db) only: screen genomes with fewer than this many stored stage-1 k-mers (small genomes: viruses, plasmids, short contigs) directly instead, by reading back their smallest this-many dense k-mers when the database is opened. Recovers small-genome sensitivity from a database built with --small-genome-screen none, at a per-sample cost proportional to the number of such genomes. 0 disables.")]
    pub screen_small_genomes: usize,


    //Hidden options that are embedded in the args but no longer used... 
    #[clap(short, hidden=true, long="pseudotax", help_heading = "ALGORITHM", help = "Pseudo taxonomic classification mode. This removes shared k-mers between species by assigning k-mers to the highest ANI species. Requires sketches with --enable-pseudotax option" )]
    pub pseudotax: bool,
    #[clap(long="ratio", hidden=true)]
    pub ratio: bool,
    #[clap(long="mme", hidden=true)]
    pub mme: bool,
    #[clap(long="mle", hidden=true)]
    pub mle: bool,
    #[clap(long="nb", hidden=true)]
    pub nb: bool,
    #[clap(long="no-ci", help = "Do not output confidence intervals", hidden=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hidden=true)]
    pub no_adj: bool,
    #[clap(long="mean-coverage", help_heading = "ALGORITHM", help = "Use the robust mean coverage estimator instead of median estimator", hidden=true )]
    pub mean_coverage: bool,

}

#[derive(Args)]
pub struct InspectArgs {
    #[clap(multiple=true, help = "Pre-sketched *.syldb/*.sylsp files.")]
    pub files: Vec<String>,
    #[clap(short='o',long="output-file", help = "Output to this file (YAML format). [default: stdout]")]
    pub out_file_name: Option<String>,

}
    
