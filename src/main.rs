use clap::Parser;
use sylph::cmdline::*;
use sylph::sketch;
use sylph::contain;
use std::panic::set_hook;

//Use this allocator when statically compiling
//instead of the default
//because the musl statically compiled binary
//uses a bad default allocator which makes the
//binary take 60% longer!!! Only affects
//static compilation though. 
#[cfg(target_env = "musl")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    set_hook(Box::new(|info| {
        if let Some(s) = info.payload().downcast_ref::<String>() {
            log::error!("{}", s);
        }
    }));
    let cli = Cli::parse();
    match cli.mode {
        Mode::Sketch(sketch_args) => sketch::sketch(sketch_args),
        Mode::Contain(contain_args) => contain::contain(contain_args),
    }
}
