use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct Cli {
    #[clap(subcommand)]
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
    #[clap(short, multiple=true)]
    pub sequences: Option<Vec<String>>,
    #[clap(short, multiple=true)]
    pub genomes: Option<Vec<String>>,
    #[clap(long="sequence-prefix", default_value = "")]
    pub sequence_output_prefix: String,
    #[clap(long="genome-prefix", default_value = "")]
    pub genome_output_prefix: String,
    #[clap(short, default_value_t = 21)]
    pub k: usize,
    #[clap(short, default_value_t = 500)]
    pub c: usize,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,

}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(short, multiple=true, required=true)]
    pub sequence_sketches: Vec<String>,
    #[clap(short, multiple=true, required=true)]
    pub genome_sketches: Vec<String>,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,
}
