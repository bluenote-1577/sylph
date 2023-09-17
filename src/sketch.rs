use crate::cmdline::*;
use std::thread;
use std::time::Duration;

use memory_stats::memory_stats;

use crate::constants::*;
use crate::seeding::*;
use crate::types::*;
use flate2::read::GzDecoder;
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use regex::Regex;
use seq_io;
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

pub fn check_vram_and_block(max_ram: usize, file: &str){
    if let Some(usage) = memory_stats() {
        let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
        if (max_ram as f64) < gb_usage_curr{
            log::debug!("Max memory reached. Blocking sketch for {}. Curr memory {}, max mem {}", file, gb_usage_curr, max_ram);
        }
        while (max_ram as f64) < gb_usage_curr{
            let five_second = Duration::from_secs(3);
            thread::sleep(five_second);
            gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
            if (max_ram as f64) >= gb_usage_curr{
                log::debug!("Sketching for {} freed", file);
            }
        }

    }
}

pub fn extract_markers(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    #[cfg(any(target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2")  && false{
            use crate::avx2_seeding::*;
            unsafe {
                extract_markers_avx2(string, kmer_vec, c, k);
            }
        } else {
            fmh_seeds(string, kmer_vec, c, k);
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        fmh_seeds(string, kmer_vec, c, k);
    }
}

pub fn extract_markers_positions(string: &[u8], kmer_vec: &mut Vec<(usize, usize, u64)>, c: usize, k: usize, contig_number: usize) {
    #[cfg(any(target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2"){
            use crate::avx2_seeding::*;
            unsafe {
                extract_markers_avx2_positions(string, kmer_vec, c, k, contig_number);
            }
        } else {
            fmh_seeds_positions(string, kmer_vec, c, k, contig_number);
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        fmh_seeds_positions(string, kmer_vec, c, k, contig_number);
    }
}

pub fn is_fastq(file: &str) -> bool {
    if file.ends_with(".fq")
        || file.ends_with(".fnq")
        || file.ends_with(".fastq")
        || file.ends_with(".fq.gz")
        || file.ends_with(".fnq.gz")
        || file.ends_with(".fastq.gz")
    {
        return true;
    } else {
        return false;
    }
}

pub fn is_fasta(file: &str) -> bool {
    if file.ends_with(".fa")
        || file.ends_with(".fna")
        || file.ends_with(".fasta")
        || file.ends_with(".fa.gz")
        || file.ends_with(".fna.gz")
        || file.ends_with(".fasta.gz")
    {
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

    if args.files.is_empty()
        && args.list_sequence.is_none()
        && args.first_pair.is_empty()
        && args.second_pair.is_empty()
    {
        error!("No input sequences found; see sylph sketch -h for help. Exiting.");
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
        if args.sample_force{
            read_inputs.push(file);
        } else if args.query_force{
            genome_inputs.push(file);
        } else if is_fastq(&file) {
            read_inputs.push(file);
        } else if is_fasta(&file) {
            genome_inputs.push(file);
        } else {
            warn!("{} does not have a fasta/fastq/gzip type extension; skipping", file);
        }
    }

    let mut max_ram = usize::MAX;
    if args.max_ram.is_some(){
        max_ram = args.max_ram.unwrap();
        if max_ram < 10{
            log::error!("Max ram must be >= 10. Exiting.");
            std::process::exit(1);
        }
    }

    if !args.first_pair.is_empty() && !args.second_pair.is_empty() {
        info!("Sketching paired sequences...");
        if args.first_pair.len() != args.second_pair.len() {
            error!("Different number of paired sequences. Exiting.");
            std::process::exit(1);
        }
        let iter_vec: Vec<usize> = (0..args.first_pair.len()).into_iter().collect();
        iter_vec.into_par_iter().for_each(|i| {
            let read_file1 = &args.first_pair[i];
            let read_file2 = &args.second_pair[i];
            check_vram_and_block(max_ram, read_file1);

            let read_sketch_opt = sketch_pair_sequences(read_file1, read_file2, args.c, args.k);
            if read_sketch_opt.is_some() {
                let pref = Path::new(&args.sample_output_dir);
                let read_sketch = read_sketch_opt.unwrap();
                let read_file_path = Path::new(&read_sketch.file_name).file_name().unwrap();
                let file_path = pref.join(&read_file_path);

                let file_path_str = format!(
                    "{}.paired{}",
                    file_path.to_str().unwrap(),
                    SAMPLE_FILE_SUFFIX
                );

                let mut read_sk_file = BufWriter::new(
                    File::create(&file_path_str)
                        .expect(&format!("{} path not valid; exiting ", file_path_str)),
                );

                let enc = SequencesSketchEncode::new(read_sketch);
                bincode::serialize_into(&mut read_sk_file, &enc).unwrap();
                info!("Sketching {} complete.", file_path_str);
            }
        });
    }

    if !read_inputs.is_empty() {
        info!("Sketching non-paired sequences...");
    }

    let iter_vec: Vec<usize> = (0..read_inputs.len()).into_iter().collect();
    iter_vec.into_par_iter().for_each(|i| {
        let pref = Path::new(&args.sample_output_dir);
        let read_file = &read_inputs[i];

        check_vram_and_block(max_ram, read_file);

        let read_sketch_opt;
        if args.sample_force {
            read_sketch_opt = sketch_sequences_needle(read_file, args.c, args.k)
        } else {
            read_sketch_opt = sketch_sequences_needle(read_file, args.c, args.k)
            //read_sketch_opt = sketch_query(args.c, args.k, args.threads, read_file);
        }
        if read_sketch_opt.is_some() {
            let read_sketch = read_sketch_opt.unwrap();
            let read_file_path = Path::new(&read_sketch.file_name).file_name().unwrap();
            let file_path = pref.join(&read_file_path);

            let file_path_str = format!("{}{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);

            let mut read_sk_file = BufWriter::new(
                File::create(&file_path_str)
                    .expect(&format!("{} path not valid; exiting ", file_path_str)),
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
        let pref = Path::new(&args.query_out_name);
        let file_path_str = format!("{}{}", pref.to_str().unwrap(), QUERY_FILE_SUFFIX);
        let all_genome_sketches = Mutex::new(vec![]);

        iter_vec.into_par_iter().for_each(|i| {
            let genome_file = &genome_inputs[i];
            if args.individual {
                let indiv_gn_sketches =
                    sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                all_genome_sketches
                    .lock()
                    .unwrap()
                    .extend(indiv_gn_sketches);
            } else {
                let genome_sketch =
                    sketch_genome(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
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



        if all_genome_sketches.lock().unwrap().is_empty(){
            warn!("No valid queries (e.g. genomes or fasta files) to sketch; {} is not output", file_path_str);
        }
        else{
            let mut genome_sk_file = BufWriter::new(
                File::create(&file_path_str).expect(&format!("{} not valid ", file_path_str)),
            );
            info!("Wrote all genome sketches to {}", file_path_str);
            bincode::serialize_into(&mut genome_sk_file, &all_genome_sketches).unwrap();
        }
    }

    info!("Finished.");
}

pub fn sketch_genome_individual(
    c: usize,
    k: usize,
    ref_file: &str,
    min_spacing: usize,
    pseudotax: bool
) -> Vec<GenomeSketch> {
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
                let mut pseudotax_track_kmers = vec![];
                let mut kmer_vec = vec![];
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                let contig_name = String::from_utf8_lossy(record.id()).to_string();
                return_genome_sketch.first_contig_name = contig_name;
                let seq = record.seq();

                extract_markers_positions(&seq, &mut kmer_vec, c, k, 0);

                let mut kmer_set = MMHashSet::default();
                let mut duplicate_set = MMHashSet::default();
                let mut new_vec = Vec::with_capacity(kmer_vec.len());
                kmer_vec.sort();
                for (_, _pos, km) in kmer_vec.iter() {
                    if !kmer_set.contains(&km) {
                        kmer_set.insert(km);
                    } else {
                        duplicate_set.insert(km);
                    }
                }
                let mut last_pos = 0;
                for (_, pos, km) in kmer_vec.iter() {
                    if !duplicate_set.contains(&km)
                    {
                        if last_pos == 0 || pos - last_pos > min_spacing{
                            new_vec.push(*km);
                            last_pos = *pos;
                        }
                        else if pseudotax{
                            pseudotax_track_kmers.push(*km);
                        }
                    }
                }
                return_genome_sketch.genome_kmers = new_vec;
                if pseudotax{
                    return_genome_sketch.pseudotax_tracked_nonused_kmers = Some(pseudotax_track_kmers);
                }
                return_vec.push(return_genome_sketch);
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return vec![];
            }
        }
        return return_vec;
    }
}

pub fn sketch_genome(
    c: usize,
    k: usize,
    ref_file: &str,
    min_spacing: usize,
    pseudotax: bool
) -> Option<GenomeSketch> {
    let reader = parse_fastx_file(&ref_file);
    let mut vec = vec![];
    let mut pseudotax_track_kmers = vec![];
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
        let mut contig_number = 0;
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                if first {
                    let contig_name = String::from_utf8_lossy(record.id()).to_string();
                    return_genome_sketch.first_contig_name = contig_name;
                    first = false;
                }
                let seq = record.seq();

                extract_markers_positions(&seq, &mut vec, c, k, contig_number);
                
                contig_number += 1
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return None;
            }
        }
        let mut kmer_set = MMHashSet::default();
        let mut duplicate_set = MMHashSet::default();
        let mut new_vec = Vec::with_capacity(vec.len());
        vec.sort();
        for (_, _, km) in vec.iter() {
            if !kmer_set.contains(&km) {
                kmer_set.insert(km);
            } else {
                duplicate_set.insert(km);
            }
        }

        let mut last_pos = 0;
        let mut last_contig = 0;
        for (contig, pos, km) in vec.iter() {
            if !duplicate_set.contains(&km){
                if last_pos == 0 || last_contig != *contig || pos - last_pos > min_spacing
                {
                    new_vec.push(*km);
                    last_contig = *contig;
                    last_pos = *pos;
                }
                else if pseudotax{
                    pseudotax_track_kmers.push(*km);
                }
            }
        }
        return_genome_sketch.genome_kmers = new_vec;
        if pseudotax{
            return_genome_sketch.pseudotax_tracked_nonused_kmers = Some(pseudotax_track_kmers);
        }
        return Some(return_genome_sketch);
    }
}

pub fn sketch_pair_sequences(
    read_file1: &str,
    read_file2: &str,
    c: usize,
    k: usize,
) -> Option<SequencesSketch> {
    let r1o = parse_fastx_file(&read_file1);
    let r2o = parse_fastx_file(&read_file2);
    let mut read_sketch = SequencesSketch::new(read_file1.to_string(), c, k, true);
    if r1o.is_err() || r2o.is_err() {
        panic!("Paired end reading failed");
    }
    let mut reader1 = r1o.unwrap();
    let mut reader2 = r2o.unwrap();

    loop {
        let n1 = reader1.next();
        let n2 = reader2.next();
        if let Some(rec1_o) = n1 {
            if let Some(rec2_o) = n2 {
                if let Ok(rec1) = rec1_o {
                    if let Ok(rec2) = rec2_o {
                        let mut temp_vec1 = vec![];
                        let mut temp_vec2 = vec![];
                        extract_markers(&rec1.seq(), &mut temp_vec1, c, k);
                        extract_markers(&rec2.seq(), &mut temp_vec2, c, k);
                        for km in temp_vec1.iter() {
                            let c = read_sketch.kmer_counts.entry(*km).or_insert(0);
                            *c += 1;
                        }
                        for km in temp_vec2 {
                            if temp_vec1.contains(&km) {
                                continue;
                            }
                            let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                            *c += 1;
                        }
                    }
                } else {
                    return None;
                }
            }
        } else {
            break;
        }
    }
    return Some(read_sketch);
}

//This did not work
//pub fn sketch_pair_sequences(read_file1: &str, read_file2: &str) {
//    let mut queue: WorkQueue<(, _)> = WorkQueue::new();
//
//    let (results_tx, results_rx) = channel();
//
//    // Create a SyncFlag to share whether or not the worker threads should
//    // keep waiting on jobs.
//    let (mut more_jobs_tx, more_jobs_rx) = new_syncflag(true);
//
//    // This Vec is just for the controller to keep track of the worker threads.
//    let mut thread_handles = Vec::new();
//    for _ in 0..15 {
//        let mut t_queue = queue.clone();
//        let t_results_tx = results_tx.clone();
//        let t_more_jobs = more_jobs_rx.clone();
//        thread_handles.push(std::thread::spawn(move || {
//            while let Some(work_input) = t_queue.wait(&t_more_jobs) {
//                // Do some work. Totally contrived in this case.
//                let result = work_input;
//                // Send the results of the work to the main thread.
//                //
//                unsafe {
//                    extract_markers_avx2(&result.0, &mut vec![], 100, 21);
//                    extract_markers_avx2(&result.1, &mut vec![], 100, 21);
//                }
//                t_results_tx.send(("x", result)).unwrap();
//            }
//        }));
//    }
//
//    let mut reader1 = ReaderQ::new(BufReader::new(File::open(&read_file1).unwrap()));
//    let mut reader2 = ReaderQ::new(BufReader::new(File::open(&read_file2).unwrap()));
//    let mut rset1 = seq_io::fastq::RecordSet::default();
//    let mut rset2 = seq_io::fastq::RecordSet::default();
//
//    loop{
//        let ok1 = reader1.read_record_set(&mut rset1);
//        let ok2 = reader2.read_record_set(&mut rset2);
//
//        queue.push_work(rset1,rset2);
////        for (rec1,rec2) in rset1.into_iter().zip(rset2.into_iter()){
////            queue.push_work((rec1, rec2));
////        }
//    }
//
////    for (record1 , record2) in records1.into_iter().zip(records2.into_iter()){
////        if record1.is_ok() && record2.is_ok() {
////            let rec1 = record1.unwrap();
////            let rec2 = record2.unwrap();
////            let seq1 = rec1.seq;
////            let seq2 = rec2.seq;
////            queue.push_work((seq1, seq2));
////            //tx.send("x").unwrap();
////        } else {
////            println!("err");
////        }
////    }
//
//    more_jobs_tx.set(false);
//
//    // Join all the threads.
//    for thread_handle in thread_handles {
//        thread_handle.join().unwrap();
//    }
//}

pub fn sketch_query(
    c: usize,
    k: usize,
    _threads: usize,
    query_file: &str,
) -> Option<SequencesSketch> {
    let read_file = query_file;
    let mut read_sketch = SequencesSketch::new(read_file.to_string(), c, k, false);
    let mut error_free = true;
    if is_fastq(read_file) {
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
            10,
            100,
            |record_set| {
                let mut vec = vec![];
                for record in record_set.into_iter() {
                    //                        dbg!(String::from_utf8(record.seq().to_vec()));
                    extract_markers(record.seq(), &mut vec, c, k);
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
            10 as u32,
            2 * 100,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                for record in record_set.into_iter() {
                    extract_markers(record.seq(), &mut vec, c, k);
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
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                let seq = record.seq();
                extract_markers(&seq, &mut vec, c, k);
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
        k,
        paired: false,
    });
}
