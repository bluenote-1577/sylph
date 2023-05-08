use crate::cmdline::*;
use crate::constants::*;
use crate::seeding::*;
use crate::types::*;
use flate2::read::GzDecoder;
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use regex::Regex;
use seq_io::fasta::{Reader as ReaderA, Record as ReccordA};
use seq_io::fastq::{Reader as ReaderQ, Record as RecordQ};
use seq_io::parallel::read_parallel;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::{prelude::*, BufReader};
use std::path::Path;
use std::sync::Mutex;

pub fn is_fastq(file: &str) -> bool {
    if file.ends_with(".fq") || file.ends_with(".fnq") || file.ends_with(".fastq") 
    || file.ends_with(".fq.gz") || file.ends_with(".fnq.gz") || file.ends_with(".fastq.gz") {
        return true;
    } else {
        return false;
    }
}

pub fn is_fasta(file: &str) -> bool {
    if file.ends_with(".fa") || file.ends_with(".fna") || file.ends_with(".fasta") 
    || file.ends_with(".fa.gz") || file.ends_with(".fna.gz") || file.ends_with(".fasta.gz") {
        return true;
    } else {
        return false;
    }
}

//Can combine two paired sketches into one. Deprecated since paired end filenames have
//no standard format...
pub fn combine_sketches(mut sketches: Vec<SequencesSketch>) -> SequencesSketch {
    assert!(!sketches.is_empty());
    let re = Regex::new(PAIR_REGEX).unwrap();
    let num_sketches = sketches.len();
    let mut first_sketch = std::mem::take(&mut sketches[0]);
    for i in 1..num_sketches {
        let x = std::mem::take(&mut sketches[i]);
        for (key, val) in x.kmer_counts.into_iter() {
            let count = first_sketch.kmer_counts.entry(key).or_insert(0);
            *count += val;
        }
    }
    for cap in re.captures_iter(&first_sketch.file_name.clone()) {
        let new_file_name = format!("{}{}", &cap[1], &cap[3]);
        first_sketch.file_name = new_file_name;
    }
    return first_sketch;
}

//Collect paired-end reads names together.
pub fn collect_pairs_file_names<'a>(sequence_names: &Option<Vec<String>>) -> Vec<Vec<&str>> {
    let mut list_pair1: HashMap<_, _> = HashMap::default();
    let mut ret_pairs = vec![];
    let re = Regex::new(PAIR_REGEX).unwrap();
    if sequence_names.is_some() {
        for read_file in sequence_names.as_ref().unwrap().iter() {
            if re.is_match(read_file) {
                let front = re.captures(read_file).unwrap().get(1).unwrap().as_str();
                let pair = list_pair1.entry(front).or_insert(vec![]);
                pair.push(read_file.as_str());
            } else {
                ret_pairs.push(vec![read_file.as_str()]);
            }
        }
    }
    for (_, files) in list_pair1 {
        if files.len() == 1 || files.len() == 2 {
            ret_pairs.push(files);
        } else if files.len() > 2 {
            log::warn!("Something went wrong with paired-end read processing. Treating pairs in {:?} as separate.", &files);
            for file in files {
                ret_pairs.push(vec![file]);
            }
        }
    }

    return ret_pairs;
}

pub fn sketch(args: SketchArgs) {
    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else {
        level = log::LevelFilter::Info;
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();


    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    if args.files.is_empty() && args.list_sequence.is_none() {
        error!("No input sequences found. Exiting.");
        std::process::exit(1);
    }

    let mut all_files = vec![];

    if args.list_sequence.is_some() {
        let file_list = args.list_sequence.unwrap();
        let file = File::open(file_list).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            all_files.push(line.unwrap());
        }
    }

    all_files.extend(args.files);

    let mut read_inputs = vec![];
    let mut genome_inputs = vec![];
    for file in all_files {
        if args.read_force {
            read_inputs.push(file);
        } else if args.genome_force {
            genome_inputs.push(file);
        } else if is_fastq(&file){
            read_inputs.push(file);
        } else if is_fasta(&file){
            genome_inputs.push(file);
        } else {
            warn!("{} does not have a fasta/fastq/gzip type extension.", file);
        }
    }

    if !read_inputs.is_empty() {
        info!("Sketching sequences...");
    }

    let iter_vec: Vec<usize> = (0..read_inputs.len()).into_iter().collect();
    iter_vec.into_par_iter().for_each(|i| {
        let pref = Path::new(&args.read_prefix);

        let read_file = &read_inputs[i];
        let read_sketch_opt;
        if args.read_force{
            read_sketch_opt = sketch_sequences_needle(read_file, args.c, args.k)
        }
        else{
            read_sketch_opt = sketch_query(args.c, args.k, args.threads, read_file);
        }
        if read_sketch_opt.is_some() {
            let read_sketch = read_sketch_opt.unwrap();
            let read_file_path = Path::new(&read_sketch.file_name).file_name().unwrap();
            let file_path = pref.join(&read_file_path);

            let file_path_str = format!("{}.prs", file_path.to_str().unwrap());

            let mut read_sk_file = BufWriter::new(
                File::create(&file_path_str)
                    .expect(&format!("{} path not valid, exiting.", file_path_str)),
            );

            let enc = SequencesSketchEncode::new(read_sketch);
            bincode::serialize_into(&mut read_sk_file, &enc).unwrap();
            info!("Sketching {} complete.", file_path_str);
        }
    });

    if !genome_inputs.is_empty() {
        info!("Sketching genomes...");
        let iter_vec: Vec<usize> = (0..genome_inputs.len()).into_iter().collect();
        let counter: Mutex<usize> = Mutex::new(0);
        let pref = Path::new(&args.genome_prefix);
        let file_path_str = format!("{}.pgs", pref.to_str().unwrap());
        let mut genome_sk_file = BufWriter::new(
            File::create(&file_path_str).expect(&format!("{} not valid", file_path_str)),
        );
        let all_genome_sketches = Mutex::new(vec![]);

        iter_vec.into_par_iter().for_each(|i| {
            let genome_file = &genome_inputs[i];
            if args.individual{
                let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file);
                all_genome_sketches.lock().unwrap().extend(indiv_gn_sketches);
            }
            else{
                let genome_sketch = sketch_genome(args.c, args.k, genome_file);
                if genome_sketch.is_some() {
                    all_genome_sketches
                        .lock()
                        .unwrap()
                        .push(genome_sketch.unwrap());
                }
            }
            let mut c = counter.lock().unwrap();
            *c += 1;
            if *c % 100 == 0 && *c != 0 {
                info!("{} genomes processed.", *c);
            }
        });

        bincode::serialize_into(&mut genome_sk_file, &all_genome_sketches).unwrap();
        info!("Wrote all genome sketches to {}", file_path_str);
    }

    info!("Finished.");
}

pub fn sketch_genome_individual(c: usize, k: usize, ref_file: &str) -> Vec<GenomeSketch> {
    let reader = parse_fastx_file(&ref_file);
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        return vec![];
    } else {
        let mut reader = reader.unwrap();
        let mut return_vec = vec![];
        while let Some(record) = reader.next() {
            let mut return_genome_sketch = GenomeSketch::default();
            return_genome_sketch.c = c;
            return_genome_sketch.k = k;
            return_genome_sketch.file_name = ref_file.to_string();
            if record.is_ok() {
                let mut kmer_vec = vec![];
                let record = record.expect(&format!("Invalid record for file {}", ref_file));
                let contig_name = String::from_utf8_lossy(record.id()).to_string();
                return_genome_sketch.first_contig_name = contig_name;
                let seq = record.seq();
                unsafe {
                    extract_markers_avx2(&seq, &mut kmer_vec, c, k);
                }
                let mut kmer_set = MMHashSet::default();
                let mut duplicate_set = MMHashSet::default();
                let mut new_vec = Vec::with_capacity(kmer_vec.len());
                for km in kmer_vec.iter() {
                    if !kmer_set.contains(&km) {
                        kmer_set.insert(km);
                    } else {
                        duplicate_set.insert(km);
                    }
                }
                for km in kmer_vec.iter() {
                    if !duplicate_set.contains(&km) {
                        new_vec.push(*km);
                    }
                }
                return_genome_sketch.genome_kmers = new_vec;
                return_vec.push(return_genome_sketch);

            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return vec![];
            }
        }
        return return_vec
    }
}

pub fn sketch_genome(c: usize, k: usize, ref_file: &str) -> Option<GenomeSketch> {
    let reader = parse_fastx_file(&ref_file);
    let mut vec = vec![];
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        return None;
    } else {
        let mut reader = reader.unwrap();
        let mut first = true;
        let mut return_genome_sketch = GenomeSketch::default();
        return_genome_sketch.c = c;
        return_genome_sketch.k = k;
        return_genome_sketch.file_name = ref_file.to_string();
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {}", ref_file));
                if first {
                    let contig_name = String::from_utf8_lossy(record.id()).to_string();
                    return_genome_sketch.first_contig_name = contig_name;
                    first = false;
                }
                let seq = record.seq();
                unsafe {
                    extract_markers_avx2_positions(&seq, &mut vec, c, k);
                }
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return None;
            }
        }
        let mut kmer_set = MMHashSet::default();
        let mut duplicate_set = MMHashSet::default();
        let mut new_vec = Vec::with_capacity(vec.len());
        vec.sort();
        for (_,km) in vec.iter() {
            if !kmer_set.contains(&km) {
                kmer_set.insert(km);
            } else {
                duplicate_set.insert(km);
            }
        }
        let mut last_pos = 0;
        for (pos,km) in vec.iter() {
            if !duplicate_set.contains(&km) && (pos - last_pos > 75 || last_pos == 0){
            //if !duplicate_set.contains(&km){
                new_vec.push(*km);
                last_pos = *pos;
            }
        }
        return_genome_sketch.genome_kmers = new_vec;
        return Some(return_genome_sketch);
    }
}

pub fn sketch_query(
    c: usize,
    k: usize,
    threads: usize,
    query_file: &str,
) -> Option<SequencesSketch> {
    let read_file = query_file;
    let mut read_sketch = SequencesSketch::new(read_file.to_string(), c, k);
    let mut error_free = true;
    if is_fastq(read_file){
        let reader;
        if read_file.contains(".gz") || read_file.contains(".gzip") {
            let file = File::open(read_file).unwrap();
            let gz_decoder: Box<dyn Read + Send> = Box::new(BufReader::new(GzDecoder::new(file)));
            reader = ReaderQ::new(gz_decoder);
        } else {
            let file = File::open(read_file).unwrap();
            let decoder: Box<dyn Read + Send> = Box::new(BufReader::new(file));
            reader = ReaderQ::new(decoder);
        }
        read_parallel(
            reader,
            threads as u32,
            100,
            |record_set| {
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        //                        dbg!(String::from_utf8(record.seq().to_vec()));
                        extract_markers_avx2(record.seq(), &mut vec, c, k);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    if result.is_ok() {
                        let (_rec_set, vec) = result.unwrap();
                        for km in vec {
                            let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                            *c += 1;
                        }
                    } else {
                        log::warn!("{} was not a valid sequence file.", query_file);
                        error_free = false;
                        break;
                    }
                }
            },
        );
    } else {
        let reader;
        if read_file.contains(".gz") || read_file.contains(".gzip") {
            let file = File::open(read_file).unwrap();
            let gz_decoder: Box<dyn Read + Send> = Box::new(BufReader::new(GzDecoder::new(file)));
            reader = ReaderA::new(gz_decoder);
        } else {
            let file = File::open(read_file).unwrap();
            let decoder: Box<dyn Read + Send> = Box::new(BufReader::new(file));
            reader = ReaderA::new(decoder);
        }
        read_parallel(
            reader,
            threads as u32,
            100,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        extract_markers_avx2(record.seq(), &mut vec, c, k);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    if result.is_ok() {
                        let (_rec_set, vec) = result.as_ref().unwrap();
                        for km in vec {
                            let c = read_sketch.kmer_counts.entry(*km).or_insert(0);
                            *c += 1;
                        }
                    } else {
                        error_free = false;
                        log::warn!("{} was not a valid sequence file.", query_file);
                    }
                }
            },
        );
    }

    if error_free {
        return Some(read_sketch);
    } else {
        return None;
    }
}

pub fn sketch_sequences_needle(read_file: &str, c: usize, k: usize) -> Option<SequencesSketch> {
    let mut kmer_map = HashMap::default();
    let ref_file = &read_file;
    let reader = parse_fastx_file(&ref_file);
    let mut vec = vec![];
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
    } else {
        let mut reader = reader.unwrap();
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {}", ref_file));
                let seq = record.seq();
                unsafe {
                    extract_markers_avx2(&seq, &mut vec, c, k);
                }
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
            }
        }
        for km in vec {
            let c = kmer_map.entry(km).or_insert(0);
            *c += 1;
        }
    }

    return Some(SequencesSketch {
        kmer_counts: kmer_map,
        file_name: read_file.to_string(),
        c, 
        k
    });
}
