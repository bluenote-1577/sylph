use crate::types::*;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use log::*;
use crate::constants::*;
use crate::cmdline::*;
use serde::{Deserialize, Serialize};

fn pipe_write(text: &str, writer: &mut Box<dyn Write + Send>){
    let result = write!(writer, "{}", text);
    match result {
        Err(e) if e.kind() == std::io::ErrorKind::BrokenPipe => {},
        _other => {},
    }
}

#[derive(Default, Deserialize, Serialize, Debug, PartialEq)]
struct SequencesSketchInspect{
    pub file_name: String,
    pub c: usize,
    pub k: usize,
    pub num_sketched_kmers: usize,
    pub approximate_number_bases: f32,
    pub mean_read_length: f64,
    pub sample_name: Option<String>,
    pub paired: bool,
}

impl SequencesSketchInspect{
    pub fn new(
        sk: SequencesSketch
    ) -> Self {
        SequencesSketchInspect{
            file_name: sk.file_name,
            num_sketched_kmers: sk.kmer_counts.len(),
            c: sk.c,
            k: sk.k,
            approximate_number_bases: (sk.mean_read_length + sk.k as f64 - 1.) as f32 / (sk.mean_read_length) as f32 * sk.c as f32 * sk.kmer_counts.len() as f32,
            sample_name: sk.sample_name,
            paired: sk.paired,
            mean_read_length: sk.mean_read_length,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct GenomeSketchInspect{
    pub file_name: String,
    pub genome_kmers_num: usize,
    pub first_contig_name: String,
    pub genome_size: usize,
}

impl GenomeSketchInspect{
    pub fn new(
        sk: GenomeSketch
    ) -> Self {
        GenomeSketchInspect{
            genome_kmers_num: sk.genome_kmers.len(),
            file_name: sk.file_name,
            first_contig_name: sk.first_contig_name,
            genome_size: sk.gn_size,
        }
    }
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Hash, PartialOrd, Eq, Ord, Default, Clone)]
pub struct DatabaseSketch{
    pub database_file: String,
    pub c: usize,
    pub k: usize,
    pub min_spacing_parameter: usize,
    pub genome_files: Vec<GenomeSketchInspect>,
}


pub fn inspect(args: InspectArgs){
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    let mut read_sketch_files = Vec::new();
    let mut genome_sketch_files = Vec::new();

    for file in args.files.iter(){
        let mut genome_sketch_good_suffix = false;
        for suff in QUERY_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                genome_sketch_good_suffix = true;
                break
            }
        }

        let mut sample_sketch_good_suffix = false;
        for suff in SAMPLE_FILE_SUFFIX_VALID{
            if file.ends_with(suff){
                sample_sketch_good_suffix = true;
                break
            }
        }

        if genome_sketch_good_suffix{
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix{
            read_sketch_files.push(file);
        } else {
            warn!(
                "{} file is not a .sylsp or .syldb file. Skipping...",
                &file
            );
        }
    }

    let mut out_writer = match args.out_file_name {
        Some(ref x) => {
            Box::new(BufWriter::new(File::create(&x).unwrap())) as Box<dyn Write + Send>
        }
        None => Box::new(BufWriter::new(std::io::stdout())) as Box<dyn Write + Send>,
    };

    let mut db_sketches_inspect = Vec::new();
    for file in genome_sketch_files.iter(){
        let db_sketch = get_db_sketch_inspect(file);
        db_sketches_inspect.push(db_sketch);
    }
    let yaml = serde_yaml::to_string(&db_sketches_inspect).unwrap();
    if !db_sketches_inspect.is_empty(){
        pipe_write(&yaml, &mut out_writer);
    }

    let mut seq_sketches_inspect = Vec::new();

    for file in read_sketch_files.iter(){
        let seq_sketch = get_seq_sketch_inspect(file);
        seq_sketches_inspect.push(seq_sketch);
    }
    let yaml = serde_yaml::to_string(&seq_sketches_inspect).unwrap();
    if !seq_sketches_inspect.is_empty(){
        pipe_write(&yaml, &mut out_writer);
    }
}

fn get_db_sketch_inspect(
    genome_sketch_file: &String,
)  -> DatabaseSketch{

    let file = File::open(genome_sketch_file).expect(&format!("The sketch `{}` could not be opened. Exiting", genome_sketch_file));
    let genome_reader = BufReader::with_capacity(10_000_000, file);
    let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
        .expect(&format!(
            "The database sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
            &genome_sketch_file
        ));
    if genome_sketches_vec.is_empty(){
        warn!(
            "The database sketch `{}` is empty. Skipping...",
            &genome_sketch_file
        );
        return DatabaseSketch::default();
    }

    info!(
        "Database file {} processed with {} genomes",
        genome_sketch_file, genome_sketches_vec.len()
    );

    let c = genome_sketches_vec.first().unwrap().c;
    let k = genome_sketches_vec.first().unwrap().k;
    let min_spacing = Some(genome_sketches_vec.first().unwrap().min_spacing);

    let inspect_sketch_vec = genome_sketches_vec.into_iter().map(|sk| GenomeSketchInspect::new(sk)).collect();

    
    let database_inspect = DatabaseSketch{
        database_file: genome_sketch_file.clone(),
        c,
        k,
        min_spacing_parameter: min_spacing.unwrap(),
        genome_files: inspect_sketch_vec,
    };

    return database_inspect;
}

fn get_seq_sketch_inspect(
    read_file: &String,
) -> SequencesSketchInspect{
    let file = File::open(read_file).expect(&format!("The sketch `{}` could not be opened. Exiting", read_file));
    let seq_reader = BufReader::with_capacity(10_000_000, file);
    let seq_sketch: SequencesSketch = bincode::deserialize_from(seq_reader)
        .expect(&format!(
            "The sequence sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
            &read_file
        ));
    info!(
        "Sequence file {} processed",
        read_file,
    );
    return SequencesSketchInspect::new(seq_sketch);
}
