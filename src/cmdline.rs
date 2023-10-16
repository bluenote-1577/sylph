use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about = "Ultrafast genome ANI queries and taxonomic profiling for genomes against metagenomic shotgun samples.\n\n--- Preparing inputs by sketching (indexing)\nsylph sketch sample1.fq sample2.fq genome1.fa genome2.fa -o genomes -d sample_dir \nls sample_dir/sample1.fq.sylsp sample_dir/sample2.fq.sylsp genomes.syldb \n\n--- Nearest neighbour containment ANI\nsylph contain *.syldb *.sylsp > all-to-all-contain.tsv\n\n--- Taxonomic profiling with relative abundances and ANI\nsylph profile *.syldb *.sylsp > all-to-all-profile.tsv", arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Sketch sequences into samples (reads) and databases (genomes). Each sample.fq -> sample.sylsp. All *.fa -> *.syldb. 
    Sketch(SketchArgs),
    /// Calculate coverage-adjusted containment ANI between databases and samples.
    Contain(ContainArgs),
    ///Species-level taxonomic profiling with relative abundance and coverage-adjusted ANI output. 
    Profile(ContainArgs),
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true, help = "fasta/fastq files; gzip optional. Each fastq file produces a sample sketch (*.sylsp) while fasta files are combined into a database (*.syldb)")]
    pub files: Vec<String>,
    #[clap(short='o',long="out-name-db", default_value = "database", help_heading = "OUTPUT", help = "Output name for database sketch (with .syldb appended)")]
    pub db_out_name: String,
    #[clap(short='d',long="sample-output-directory", default_value = "./", help_heading = "OUTPUT", help = "Output directory for sample sketches")]
    pub sample_output_dir: String,
    #[clap(short,long="individual-records", help_heading = "INPUT", help = "Use individual records (contigs) for database construction")]
    pub individual: bool,
    #[clap(short,long="sample-force", help_heading = "INPUT", help = "Ignore fasta/fastq extension and consider inputs as samples (reads)")]
    pub sample_force: bool,
    #[clap(short='g', long="db-force", help_heading = "INPUT", help = "Ignore fasta/fastq extension and consider inputs as genomes for a database")]
    pub db_force: bool,
    #[clap(short,long="list", help_heading = "INPUT", help = "Use files in a newline delimited text file as inputs")]
    pub list_sequence: Option<String>,
    #[clap(short, default_value_t = 31,help_heading = "ALGORITHM", help ="Value of k. Only k = 21, 31 are currently supported")]
    pub k: usize,
    #[clap(short, default_value_t = 200, help_heading = "ALGORITHM", help = "Subsampling rate")]
    pub c: usize,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(long="ram-barrier", help = "Stop multi-threaded read sketching when (virtual) RAM is past this value (in GB). Does NOT guarantee max RAM limit", hidden=true)]
    pub max_ram: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(long="disable-profiling", help_heading = "ALGORITHM", help = "Disable profiling capabilities for databases; may decrease size and make containment slightly faster")]
    pub no_pseudotax: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "ALGORITHM", help = "Minimum spacing between selected k-mers on the genomes")]
    pub min_spacing_kmer: usize,
    #[clap(short='1',long="first-pair", multiple=true, help_heading = "INPUT", help = "First pairs in paired end reads e.g. S1_1.fq S2_1.fq")]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pair", multiple=true, help_heading = "INPUT", help = "Second pairs in paired end reads e.g. S1_2.fq S2_2.fq")]
    pub second_pair: Vec<String>,
}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(multiple=true, help = "Pre-sketched *.syldb/*.sylsp files. Raw fastq/fasta are allowed and will be automatically sketched to .sylsp/.syldb")]
    pub files: Vec<String>,
    #[clap(long,default_value_t = 3., help_heading = "ALGORITHM", help = "Minimum k-mer multiplicity needed for coverage correction. Higher values gives more precision but lower sensitivity")]
    pub min_count_correct: f64,
    #[clap(long,default_value_t = 50., help_heading = "ALGORITHM", help = "Exclude genomes with less than this number of sampled k-mers")]
    pub min_number_kmers: f64,
    #[clap(short, long="minimum-ani", help_heading = "ALGORITHM", help = "Minimum adjusted ANI to consider (0-100). Default is 90 for contain and 95 for profile" )]
    pub minimum_ani: Option<f64>,
    #[clap(short, default_value_t = 3, help = "Number of threads")]
    pub threads: usize,
    #[clap(short='r', long="sample-threads", help = "Number of samples to be processed concurrently. Default: (# of total threads / 3) + 1 for profile, 1 for contain")]
    pub sample_threads: Option<usize>,
    #[clap(long="trace", help = "Trace output (caution: very verbose)")]
    pub trace: bool,
    #[clap(long="debug", help = "Debug output")]
    pub debug: bool,

    #[clap(short='u', long="estimate-unknown", help_heading = "ALGORITHM", help = "Estimates true coverage and scales sequence abundance in `profile` by estimated unknown sequence percentage" )]
    pub estimate_unknown: bool,

    #[clap(short, default_value_t = 200, help_heading = "SKETCHING", help = "Subsampling rate. Does nothing for pre-sketched files")]
    pub c: usize,
    #[clap(short, default_value_t = 31, help_heading = "SKETCHING", help = "Value of k. Only k = 21, 31 are currently supported. Does nothing for pre-sketched files")]
    pub k: usize,
    #[clap(short,long="individual-records", help_heading = "SKETCHING", help = "Use individual records (e.g. contigs) for database construction instead. Does nothing for pre-sketched files")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 30, help_heading = "SKETCHING", help = "Minimum spacing between selected k-mers on the database genomes. Does nothing for pre-sketched files")]
    pub min_spacing_kmer: usize,

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
    #[clap(long="no-ci", help_heading = "OUTPUT", help = "Do not output confidence intervals", hidden=true)]
    pub no_ci: bool,
    #[clap(long="no-adjust", hidden=true)]
    pub no_adj: bool,

}
