use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "Ultrafast genome similarity ANI queries for genomes against metagenomic shotgun samples.\n\nsylph sketch sample1.fq sample2.fq genome1.fa genome2.fa\nsylph contain *.sylqueries *.sylsamples > all-to-all-results.tsv", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch sequences into samples (reads) and queries (genomes). 
    Sketch(SketchArgs),
    /// Calculate coverage-adjusted containment ANI between queries and samples.
    Contain(ContainArgs),
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help = "fasta or fastq file; gzip optional. fastq files will be considered samples (*.sylsample) and fasta files as queries (*.sylqueries) by default")]
    pub files: Vec<String>,
    #[clap(short='o',long="query-output-prefix", default_value = "sylph_queries", help_heading = "OUTPUT", help = "Prefix for query sketches")]
    pub query_prefix: String,
    #[clap(short,long="sample-output-prefix", default_value = "", help_heading = "OUTPUT", help = "Prefix for sample sketches")]
    pub sample_prefix: String,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (e.g. contigs) as queries instead")]
    pub individual: bool,
    #[clap(long="sample-force", help_heading = "INPUT", help = "Ignore fasta/fastq and force all inputs to be samples")]
    pub sample_force: bool,
    #[clap(long="query-force", help_heading = "INPUT", help = "Ignore fasta/fastq and force all inputs to be queries. Does not work for paired reads (-1, -2 options)")]
    pub query_force: bool,
    #[clap(short,long="list", help_heading = "INPUT", help = "Use files in a newline delimited text file as inputs")]
    pub list_sequence: Option<String>,
    #[clap(short, default_value_t = 31,help_heading = "ALGORITHM", help ="Value of k. Only k = 21, 31 are currently supported")]
    pub k: usize,
    #[clap(short, default_value_t = 100, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output for debugging")]
    pub trace: bool,
    #[clap(long="min-spacing", default_value_t = 150, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the queries")]
    pub min_spacing_kmer: usize,
    #[clap(short='1',long="first-pair", multiple=true, help_heading = "INPUT", help = "First pairs in paired end reads e.g. S1_1.fq S2_1.fq")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pair", multiple=true, help_heading = "INPUT", help = "Second pairs in paried end reads e.g. S1_2.fq S2_2.fq")]
    pub second_pair: Vec<String>,
}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(multiple=true, help = "Pre-sketched query or sample files. Raw fastq/fasta also allowed but presketching is recommended; see sylph sketch for more info")]
    pub files: Vec<String>,
    #[clap(short, default_value_t = 31, help_heading = "ALGORITHM", help = "Value of k. Only k = 21, 31 are currently supported. Only for raw fasta/fastq")]
    pub k: usize,
    #[clap(short, default_value_t = 100, help_heading = "ALGORITHM", help = "Subsampling rate. Only for raw fasta/fastq")]
    pub c: usize,
    #[clap(short, long="minimum-ani", default_value_t = 78., help_heading = "OUTPUT", help = "Minimum adjusted ANI to output (0-100)" )]
    pub minimum_ani: f64,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="trace", help = "Trace output for debugging")]
    pub trace: bool,
    #[clap(long="ratio", hidden=true)]
    pub ratio: bool,
    #[clap(long="mme", hidden=true)]
    pub mme: bool,
    #[clap(long="mle", hidden=true)]
    pub mle: bool,
    #[clap(long="nb", hidden=true)]
    pub nb: bool,
    #[clap(long="no-ci", help_heading = "OUTPUT", help = "Do not output confidence intervals")]
    pub no_ci: bool,
    #[clap(long="no-adjust", hidden=true)]
    pub no_adj: bool,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (e.g. contigs) as queries instead. Only for raw fasta/fastq")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 150, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the queries. Only for raw fasta/fastq")]
    pub min_spacing_kmer: usize
}
