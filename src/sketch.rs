use crate::cmdline::*;
use smallvec::SmallVec;
use smallvec::smallvec;
use std::fs;
use std::thread;
use std::time::Duration;
use memory_stats::memory_stats;
use fxhash::FxHashMap;

use crate::constants::*;
use crate::seeding::*;
use crate::types::*;
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::{prelude::*, BufReader};
use std::path::Path;
use std::sync::Mutex;
type Marker = u64;

pub fn check_vram_and_block(max_ram: usize, file: &str){
    if let Some(usage) = memory_stats() {
        let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
        if (max_ram as f64) < gb_usage_curr{
            log::debug!("Max memory reached. Blocking sketch for {}. Curr memory {}, max mem {}", file, gb_usage_curr, max_ram);
        }
        while (max_ram as f64) < gb_usage_curr{
            let five_second = Duration::from_secs(1);
            thread::sleep(five_second);
            if let Some(usage) = memory_stats() {
                gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
                if (max_ram as f64) >= gb_usage_curr{
                    log::debug!("Sketching for {} freed", file);
                }
            }
            else{
                break;
            }
        }

    }
}

pub fn extract_markers(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    #[cfg(any(target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2"){
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

pub fn sketch(args: SketchArgs) {
    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else if args.debug {
        level = log::LevelFilter::Debug;
    }
    else{
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
        && args.genomes.is_none()
        && args.reads.is_none()
        && args.list_genomes.is_none()
        && args.list_reads.is_none()
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
        if is_fastq(&file) {
            read_inputs.push(file);
        } else if is_fasta(&file) {
            genome_inputs.push(file);
        } else {
            warn!("{} does not have a fasta/fastq/gzip type extension; skipping", file);
        }
    }

    if let Some(genomes_syl_in) = args.genomes{
        for gn_file in genomes_syl_in{
            genome_inputs.push(gn_file);
        }
    }
    if let Some(reads_syl_in) = args.reads{
        for rd_file in reads_syl_in{
            read_inputs.push(rd_file);
        }
    }

    if args.list_reads.is_some() {
        let file_reads = args.list_reads.unwrap();
        let file = File::open(file_reads).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            read_inputs.push(line.unwrap());
        }
    }

    if args.list_genomes.is_some() {
        let file_genomes = args.list_genomes.unwrap();
        let file = File::open(file_genomes).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            genome_inputs.push(line.unwrap());
        }
    }

    let mut max_ram = usize::MAX;
    if args.max_ram.is_some(){
        max_ram = args.max_ram.unwrap();
        if max_ram < 7{
            log::error!("Max ram must be >= 7. Exiting.");
            std::process::exit(1);
        }
    }

    if args.first_pair.is_empty() && !args.second_pair.is_empty(){
        error!("Different number of paired sequences. Exiting.");
        std::process::exit(1);

    }
    if !args.first_pair.is_empty() && args.second_pair.is_empty(){
        error!("Different number of paired sequences. Exiting.");
        std::process::exit(1);
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
                let res = fs::create_dir_all(&args.sample_output_dir);
                if res.is_err(){
                    error!("Could not create directory at {}", args.sample_output_dir);
                    std::process::exit(1);
                }
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
        std::fs::create_dir_all(pref).expect("Could not create directory for output sample files (-d). Exiting...");

        let read_file = &read_inputs[i];

        check_vram_and_block(max_ram, read_file);

        let read_sketch_opt;
        read_sketch_opt = sketch_sequences_needle(read_file, args.c, args.k);

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
        let pref = Path::new(&args.db_out_name);
        let file_path_str = format!("{}{}", pref.to_str().unwrap(), QUERY_FILE_SUFFIX);
        let path = std::path::Path::new(&file_path_str);
        let prefix = path.parent().unwrap();
        std::fs::create_dir_all(prefix).expect("Could not create directory for output database file (-o). Exiting...");
        let all_genome_sketches = Mutex::new(vec![]);

        iter_vec.into_par_iter().for_each(|i| {
            let genome_file = &genome_inputs[i];
            if args.individual {
                let indiv_gn_sketches =
                    sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, !args.no_pseudotax);
                all_genome_sketches
                    .lock()
                    .unwrap()
                    .extend(indiv_gn_sketches);
            } else {
                let genome_sketch =
                    sketch_genome(args.c, args.k, genome_file, args.min_spacing_kmer, !args.no_pseudotax);
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
            warn!("No valid genomes to sketch; {} is not output", file_path_str);
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

                return_genome_sketch.gn_size = record.seq().len();
                return_genome_sketch.genome_kmers = new_vec;
                return_genome_sketch.min_spacing = min_spacing;
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

                return_genome_sketch.gn_size += seq.len();
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
        return_genome_sketch.min_spacing = min_spacing;
        if pseudotax{
            return_genome_sketch.pseudotax_tracked_nonused_kmers = Some(pseudotax_track_kmers);
        }
        return Some(return_genome_sketch);
    }
}

#[inline]
fn pair_kmer(s1: &[u8], s2: &[u8]) -> Option<([Marker;2], [Marker;2])>{
    let k = std::mem::size_of::<Marker>() * 4 ;
    if s1.len() < 2 * k + 1 || s2.len() < 2 * k + 1{
        return None
    }
    else{
        let mut kmer_f = 0;
        let mut kmer_g = 0;
        let mut kmer_r = 0;
        let mut kmer_t = 0;
        for i in 0..k{
            let nuc_1 = BYTE_TO_SEQ[s1[2*i] as usize] as Marker;
            let nuc_2 = BYTE_TO_SEQ[s2[2*i] as usize] as Marker;
            let nuc_3 = BYTE_TO_SEQ[s1[1+2*i] as usize] as Marker;
            let nuc_4 = BYTE_TO_SEQ[s2[1+2*i] as usize] as Marker;

            kmer_f <<= 2;
            kmer_f |= nuc_1;

            kmer_r <<= 2;
            kmer_r |= nuc_2;

            kmer_g <<= 2;
            kmer_g |= nuc_3;

            kmer_t <<= 2;
            kmer_t |= nuc_4;

        }
        return Some(([kmer_f, kmer_r], [kmer_g, kmer_t]));
    }
}

fn dup_removal_lsh(read_sketch: &mut SequencesSketch,
                   kmer_to_pair_table: &mut FxHashMap<u64, (SmallVec<[[Marker;2];1]>, SmallVec<[[Marker;2];1]>)>,
                   km: &u64,
                   kmer_pair: Option<([Marker;2],[Marker;2])>){
    let c = read_sketch.kmer_counts.entry(*km).or_insert(0);
    if *c < MAX_DEDUP_COUNT{ 
        if let Some(doublepairs) = kmer_pair{
            if kmer_to_pair_table.contains_key(km){ 
                let mut ret = false;
                let tables = kmer_to_pair_table.get_mut(km).unwrap();
                if tables.0.contains(&doublepairs.0){
                    ret = true;
                }
                else{
                    tables.0.push(doublepairs.0);
                }
                if tables.1.contains(&doublepairs.1){
                    ret = true;
                }
                else{
                    tables.1.push(doublepairs.1);
                }
                if ret{
                    return
                }
            }
            else{
                kmer_to_pair_table.insert(*km, (smallvec![doublepairs.0], smallvec![doublepairs.1]));
            }
        }
    }
    *c += 1;
    if *c == MAX_DEDUP_COUNT{
        kmer_to_pair_table.remove(km);
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
    let mut read_sketch = SequencesSketch::new(read_file1.to_string(), c, k, true, String::new(), 0.);
    if r1o.is_err() || r2o.is_err() {
        panic!("Paired end reading failed");
    }
    let mut reader1 = r1o.unwrap();
    let mut reader2 = r2o.unwrap();

    let mut kmer_to_pair_table : FxHashMap<u64,
    (SmallVec<[[Marker;2];1]>,
     SmallVec<[[Marker;2];1]>)> = FxHashMap::default();

    let mut mean_read_length:f64 = 0.;
    let mut counter:f64 = 0.;

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
                        let kmer_pair = pair_kmer(&rec1.seq(), &rec2.seq());

                        //moving average
                        counter += 1.;
                        mean_read_length = mean_read_length + 
                            ((rec1.seq().len() as f64) - mean_read_length) / counter;

                        for km in temp_vec1.iter() {
                            dup_removal_lsh(&mut read_sketch, &mut kmer_to_pair_table, km, kmer_pair); 
                            
                        }
                        for km in temp_vec2.iter() {
                            if temp_vec1.contains(km) {
                                continue;
                            }
                            dup_removal_lsh(&mut read_sketch, &mut kmer_to_pair_table, km, kmer_pair); 
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
    read_sketch.mean_read_length = mean_read_length;
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

pub fn sketch_sequences_needle(read_file: &str, c: usize, k: usize) -> Option<SequencesSketch> {
    let mut kmer_map = HashMap::default();
    let ref_file = &read_file;
    let reader = parse_fastx_file(&ref_file);
    let mut vec = vec![];
    let mut mean_read_length = 0.;
    let mut counter = 0.;
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
    } else {
        let mut reader = reader.unwrap();
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                let seq = record.seq();
                extract_markers(&seq, &mut vec, c, k);
                //moving average
                counter += 1.;
                mean_read_length = mean_read_length + 
                    ((seq.len() as f64) - mean_read_length) / counter;

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
        sample_name: String::new(),
        mean_read_length
    });
}
