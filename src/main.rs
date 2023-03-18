use clap::Parser;
use prita::cmdline::*;
use prita::sketch;
use prita::contain;

fn main() {
    simple_logger::SimpleLogger::new().
        with_level(log::LevelFilter::Trace)
        .init().unwrap();
    
    let cli = Cli::parse();
    match cli.mode {
        Mode::Sketch(sketch_args) => sketch::sketch(sketch_args),
        Mode::Contain(contain_args) => contain::contain(contain_args),
    }
}
