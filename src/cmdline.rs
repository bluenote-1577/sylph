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
    /// Ex: sylph sketch reads.fq genomes.fa
    Sketch(SketchArgs),
    /// Calculate coverage-adjusted containment ANI between queries and samples. Can also operate
    /// on raw fastq/fasta files. 
    /// Ex: sylph contain reads.fq.sylsample sylph_queries.sylqueries
    Contain(ContainArgs),
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help = "fasta or fastq file; gzip optional. fastq files will be considered samples (*.sylsample) and fasta files as queries (*.sylqueries) by default")]
    pub files: Vec<String>,
    #[clap(short='o',long="query-out-name", default_value = "sylph_queries", help_heading = "OUTPUT", help = "Output name for query sketch. All queries (i.e. genomes) will be aggregated into one file")]
    pub query_out_name: String,
    #[clap(short='d',long="sample-output-directory", default_value = "", help_heading = "OUTPUT", help = "Output directory for sample sketches. Each sample (i.e. reads) is written to its own file")]
    pub sample_prefix: String,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (e.g. contigs) as queries instead")]
    pub individual: bool,
    #[clap(short, long="sample-force", help_heading = "INPUT", help = "Ignore fasta/fastq extension and force all inputs to be samples (i.e. reads)")]
    pub sample_force: bool,
    #[clap(short, long="query-force", help_heading = "INPUT", help = "Ignore fasta/fastq extension and force all inputs to be queries (i.e. genomes). Does not work for paired reads (-1, -2 options)")]
    pub query_force: bool,
    #[clap(short,long="list", help_heading = "INPUT", help = "Use files in a newline delimited text file as inputs")]
    pub list_sequence: Option<String>,
    #[clap(short, default_value_t = 31,help_heading = "ALGORITHM", help ="Value of k. Only k = 21, 31 are currently supported")]
    pub k: usize,
    #[clap(short, default_value_t = 100, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,

    #[clap(long="ram-barrier", help = "Stop multi-threaded read sketching when (virtual) RAM is past this value (in GB). Does NOT guarantee max RAM limit", hidden=true)]
    pub max_ram: Option<usize>,

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
    #[clap(multiple=true, help = "Pre-sketched query or sample files. Raw fastq/fasta also allowed but presketching is recommended; see sylph sketch for more info on file extension handling")]
    pub files: Vec<String>,
    #[clap(short, default_value_t = 31, help_heading = "ALGORITHM", help = "Value of k. Only k = 21, 31 are currently supported. -k does nothing for pre-sketched files")]
    pub k: usize,
    #[clap(short, default_value_t = 100, help_heading = "ALGORITHM", help = "Subsampling rate. -c does nothing for pre-sketched files.")]
    pub c: usize,
    #[clap(long,default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum k-mer multiplicity needed for coverage correction. Higher gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(short, long="minimum-ani", help_heading = "OUTPUT", help = "Minimum adjusted ANI to output (0-100). Default is 90; if --pseudotax is enabled then default is 96" )]
    pub minimum_ani: Option<f64>,
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
    #[clap(long="no-ci", help_heading = "OUTPUT", help = "Do not output confidence intervals", hidden=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hidden=true)]
    pub no_adj: bool,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (e.g. contigs) as queries instead. Only for raw fasta/fastq")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 150, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the queries. Only for raw fasta/fastq")]
    pub min_spacing_kmer: usize,
    #[clap(short, long="pseudotax", help_heading = "ALGORITHM", help = "Pseudo taxonomic classification mode. This removes shared k-mers between species by assigning k-mers to the highest ANI species" )]
    pub pseudotax: bool,
}
