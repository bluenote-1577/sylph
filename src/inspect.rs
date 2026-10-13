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

impl From<SequencesSketch> for SequencesSketchInspect{
    fn from(
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

impl From<GenomeSketch> for GenomeSketchInspect {
    fn from(
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

#[derive(Deserialize, Serialize, Debug, PartialEq, Default, Clone)]
pub struct Syl2DbGenomeInspect{
    pub file_name: String,
    pub first_contig_name: String,
    pub gn_size: usize,
    pub min_spacing: usize,
    pub has_pseudotax: bool,
    // Sparse/screen k-mer count -- NOT the dense count (unlike .syldb's
    // genome_kmers_num above): cheap to report without decoding the dense
    // Golomb-Rice block, which is the whole point of the format.
    pub screen_kmers_num: usize,
}

#[derive(Deserialize, Serialize, Debug, PartialEq, Default, Clone)]
pub struct Syl2DbInspect{
    pub database_file: String,
    pub file_size_bytes: u64,
    pub c: usize,
    pub k: usize,
    pub screen_c: usize,
    /// Densest stage-1 rate anything in the database is screened at: equal to
    /// `screen_c` unless small genomes forced a denser one (`loosen` collapses
    /// `screen_c` itself; `band` keeps `screen_c` nominal and reports the denser
    /// rate only here).
    pub effective_screen_c: usize,
    /// Distinct band-1 keys, and genomes owning them: 0 unless the database was
    /// built with `--small-genome-screen band`.
    pub band_keys: usize,
    pub band_genomes: usize,
    pub num_genomes: usize,
    pub genomes: Vec<Syl2DbGenomeInspect>,
}

#[derive(Debug, Default)]
struct DatabaseVisitor {
    c: Option<usize>,
    k: Option<usize>,
    min_spacing: Option<usize>,
    sketches: Vec<GenomeSketchInspect>,
}

impl<'de> serde::de::Visitor<'de> for DatabaseVisitor {
    type Value = Self;
    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("sequence of struct GenomeSketch")
    }

    fn visit_seq<S>(mut self, mut seq: S) -> Result<Self, S::Error>
        where
            S: serde::de::SeqAccess<'de>,
        {
            while let Some(value) = seq.next_element::<GenomeSketch>()? {
                self.c.get_or_insert(value.c);
                self.k.get_or_insert(value.k);
                self.min_spacing.get_or_insert(value.min_spacing);
                self.sketches.push(value.into());
            }

            Ok(self)
        }
}

impl<'de> Deserialize<'de> for DatabaseVisitor {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where D: serde::de::Deserializer<'de>
    {
        let inspector = DatabaseVisitor::default();
        deserializer.deserialize_seq(inspector)
    }
}


pub fn inspect(args: InspectArgs){
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    let mut read_sketch_files = Vec::new();
    let mut genome_sketch_files = Vec::new();
    let mut two_stage_db_files = Vec::new();

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

        if file.ends_with(TWO_STAGE_DB_SUFFIX){
            two_stage_db_files.push(file);
        } else if genome_sketch_good_suffix{
            genome_sketch_files.push(file);
        } else if sample_sketch_good_suffix{
            read_sketch_files.push(file);
        } else {
            warn!(
                "{} file is not a .sylsp, .syldb, or .syl2db file. Skipping...",
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

    let mut syl2db_inspect = Vec::new();
    for file in two_stage_db_files.iter(){
        syl2db_inspect.push(get_syl2db_inspect(file));
    }
    let yaml = serde_yaml::to_string(&syl2db_inspect).unwrap();
    if !syl2db_inspect.is_empty(){
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

    let visitor: DatabaseVisitor = bincode::deserialize_from(genome_reader)
        .expect(&format!(
            "The database sketch `{}` is not a valid sketch. Perhaps it is an older, incompatible version ",
            &genome_sketch_file
        ));
    if visitor.sketches.is_empty() {
        warn!(
            "The database sketch `{}` is empty. Skipping...",
            &genome_sketch_file
        );
        return DatabaseSketch::default();
    }

    info!(
        "Database file {} processed with {} genomes",
        genome_sketch_file, visitor.sketches.len()
    );

    DatabaseSketch{
        database_file: genome_sketch_file.clone(),
        c: visitor.c.unwrap(),
        k: visitor.k.unwrap(),
        min_spacing_parameter: visitor.min_spacing.unwrap(),
        genome_files: visitor.sketches,
    }
}

fn get_syl2db_inspect(
    path: &String,
) -> Syl2DbInspect{
    let db = crate::twostage_db::open_file(path)
        .unwrap_or_else(|e| panic!("{} is not a valid two-stage database: {}", path, e));
    let file_size_bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);

    let genomes: Vec<Syl2DbGenomeInspect> = (0..db.len() as u32)
        .map(|g| {
            let meta = db.genome_meta(g);
            Syl2DbGenomeInspect{
                file_name: meta.file_name.clone(),
                first_contig_name: meta.first_contig_name.clone(),
                gn_size: meta.gn_size,
                min_spacing: meta.min_spacing,
                has_pseudotax: meta.has_pseudotax,
                screen_kmers_num: db.screen_index.screen_kmers(g) as usize,
            }
        })
        .collect();

    info!(
        "Two-stage database file {} processed with {} genomes",
        path, genomes.len()
    );

    Syl2DbInspect{
        database_file: path.clone(),
        file_size_bytes,
        c: db.c,
        k: db.k,
        screen_c: db.screen_c,
        effective_screen_c: db.effective_screen_c(),
        band_keys: db.screen_index.band.as_ref().map_or(0, |b| b.num_keys()),
        band_genomes: db.screen_index.band.as_ref().map_or(0, |b| b.num_genomes()),
        num_genomes: db.len(),
        genomes,
    }
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
    seq_sketch.into()
}
