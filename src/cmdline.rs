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
    #[clap(long="so", default_value = "")]
    pub sequence_output_prefix: String,
    #[clap(long="go", default_value = "")]
    pub genome_output_prefix: String,
    #[clap(long="gl")]
    pub genome_list: Option<String>,
    #[clap(long="sl")]
    pub sequence_list: Option<String>,
    #[clap(short, default_value_t = 31)]
    pub k: usize,
    #[clap(short, default_value_t = 1000)]
    pub c: usize,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,
    #[clap(long="trace")]
    pub trace: bool

}

#[derive(Args)]
pub struct ContainArgs {
    #[clap(short, multiple=true)]
    pub sequence_sketches: Option<Vec<String>>,
    #[clap(short, multiple=true)]
    pub genome_sketches: Option<Vec<String>>,
    #[clap(long="gf")]
    pub genome_folder: Option<String>,
    #[clap(long="sf")]
    pub sequence_folder: Option<String>,
    #[clap(short, default_value_t = 31)]
    pub k: usize,
    #[clap(short, default_value_t = 1000)]
    pub c: usize,
    #[clap(short, default_value_t = 3)]
    pub threads: usize,
    #[clap(long="trace")]
    pub trace: bool,
    #[clap(long="ratio")]
    pub ratio: bool,
    #[clap(long="mme")]
    pub mme: bool,
    #[clap(long="nb")]
    pub nb: bool,
    #[clap(long="ci")]
    pub ci: bool,


}
