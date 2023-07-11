use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None, arg_required_else_help = true, disable_help_subcommand = true)]
pub struct Cli {
    #[clap(subcommand,)]
    pub mode: Mode,
}

#[derive(Subcommand)]
pub enum Mode {
    /// Adds files to myapp
    Sketch(SketchArgs),
    Contain(ContainArgs),
}


#[derive(Args, Default)]
pub struct SketchArgs {
    #[clap(multiple=true)]
    pub files: Vec<String>,
    #[clap(short='o',long="output-genome-prefix", default_value = "sylph_queries")]
    pub genome_prefix: String,
    #[clap(short,long="reads-output-prefix", default_value = "")]
    pub read_prefix: String,
    #[clap(short,long="individual-records")]
    pub individual: bool,
    #[clap(long="read-force")]
    pub read_force: bool,
    #[clap(long="genome-force")]
    pub genome_force: bool,
    #[clap(short,long="list-sequence")]
    pub list_sequence: Option<String>,
    #[clap(short, default_value_t = 31)]
    pub k: usize,
    #[clap(short, default_value_t = 100)]
    pub c: usize,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,
    #[clap(long="trace")]
    pub trace: bool,
    #[clap(long="min-spacing", default_value_t = 150)]
    pub min_spacing_kmer: usize,
    #[clap(short='1',long="first-pair", multiple=true)]
    pub first_pair: Vec<String>,
    #[clap(short='2',long="second-pair", multiple=true)]
    pub second_pair: Vec<String>,
}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(multiple=true)]
    pub files: Vec<String>,
    #[clap(short, default_value_t = 31)]
    pub k: usize,
    #[clap(short, default_value_t = 100)]
    pub c: usize,
    #[clap(short, long="minimum-ani", default_value_t = 0.78)]
    pub minimum_ani: f64,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,
    #[clap(long="trace")]
    pub trace: bool,
    #[clap(long="ratio", hidden=true)]
    pub ratio: bool,
    #[clap(long="mme", hidden=true)]
    pub mme: bool,
    #[clap(long="mle", hidden=true)]
    pub mle: bool,
    #[clap(long="nb", hidden=true)]
    pub nb: bool,
    #[clap(long="no-ci")]
    pub no_ci: bool,
    #[clap(long="no-adjust")]
    pub no_adj: bool,
    #[clap(short,long="individual-records")]
    pub individual: bool,
    #[clap(long="min-spacing", default_value_t = 100)]
    pub min_spacing_kmer: usize


}
