use crate::seeding::*;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::types::*;
use log::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use seq_io::fasta::{Reader as ReaderA, Record as ReccordA};
use seq_io::fastq::{Reader as ReaderQ, Record as RecordQ};
use seq_io::parallel::read_parallel;
use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::Mutex;

pub fn sketch_sequences_needle(read_file: &str) -> SequencesSketch {
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
                    extract_markers_avx2(&seq, &mut vec);
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

    return SequencesSketch {
        kmer_counts: kmer_map,
        file_name: read_file.to_string(),
    };
}

pub fn sketch_sequences(read_file: &str) -> SequencesSketch {
    let mut read_sketch = SequencesSketch::new(read_file.to_string());
    if read_file.contains(".fq") || read_file.contains(".fastq") {
        let reader = ReaderQ::from_path(read_file).unwrap();
        read_parallel(
            reader,
            20,
            20,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        extract_markers_avx2(record.seq(), &mut vec);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    let (_rec_set, vec) = result.unwrap();
                    for km in vec {
                        let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                        *c += 1;
                    }
                }
            },
        );
    } else {
        let reader = ReaderA::from_path(read_file).unwrap();
        read_parallel(
            reader,
            1000,
            10,
            |record_set| {
                // this function does the heavy work
                let mut vec = vec![];
                unsafe {
                    for record in record_set.into_iter() {
                        extract_markers_avx2(record.seq(), &mut vec);
                    }
                }

                return vec;
            },
            |record_sets| {
                while let Some(result) = record_sets.next() {
                    let (_rec_set, vec) = result.unwrap();
                    for km in vec {
                        let c = read_sketch.kmer_counts.entry(km).or_insert(0);
                        *c += 1;
                    }
                }
            },
        );
    }

    return read_sketch;
}

pub fn sketch_genomes(files: Vec<String>, list_file: bool, derep: bool) -> GenomesSketch {
    let genome_files;
    if list_file {
        let file = File::open(&files[0]).unwrap();
        let reader = BufReader::new(file);
        let mut temp_vec = vec![];
        for line in reader.lines() {
            temp_vec.push(line.unwrap().trim().to_string());
        }
        genome_files = temp_vec;
    } else {
        genome_files = files;
    }
    let kmer_maps: Mutex<Vec<_>> = Mutex::new(vec![HashMap::default(); genome_files.len()]);
    let index_vec = (0..genome_files.len()).collect::<Vec<usize>>();
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &genome_files[i];
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
                        extract_markers_avx2(&seq, &mut vec);
                    }
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                }
            }
            let mut kmer_map = HashMap::default();
            for km in vec {
                let c = kmer_map.entry(km).or_insert(0);
                *c += 1;
            }
            let mut locked = kmer_maps.lock().unwrap();
            locked[i] = kmer_map;
        }
    });

    let gn_sketch = GenomesSketch {
        genome_kmer_counts: kmer_maps.into_inner().unwrap(),
        file_names: genome_files,
    };

    let ret;
    if derep {
        ret = derep_genome_sketch(gn_sketch);
    } else {
        ret = gn_sketch;
    }

    return ret;
}

fn derep_genome_sketch(genome_sketch: GenomesSketch) -> GenomesSketch {
    let mut used_genomes: HashSet<usize> = HashSet::default();
    let mut rep_genomes: HashSet<_> = HashSet::default();
    let mut kmer_to_genome_table: HashMap<u64, Vec<_>> = HashMap::default();
    let num_gn = genome_sketch.file_names.len();
    {
        let genome_kmers = &genome_sketch.genome_kmer_counts;

        for (i, kmer_map) in genome_kmers.iter().enumerate() {
            for kmer in kmer_map.keys() {
                let genomes = kmer_to_genome_table.entry(*kmer).or_insert(vec![]);
                genomes.push(i);
            }
        }

        for i in 0..genome_kmers.len() {
            if used_genomes.contains(&i) {
                continue;
            }
            rep_genomes.insert(i);
            let mut genome_containment_count: HashMap<_, _> = HashMap::default();
            let kmers = genome_kmers[i].keys();
            for kmer in kmers {
                let corresp_genomes = &kmer_to_genome_table[kmer];
                for j in corresp_genomes {
                    if used_genomes.contains(j) {
                        continue;
                    }
                    let c = genome_containment_count.entry(*j).or_insert(0);
                    *c += 1;
                }
            }

            for (j, count) in genome_containment_count {
                let ratio = count as f64
//                    / ((genome_kmers[i].len() as f64 + genome_kmers[j].len() as f64)/2.);
                    / f64::min(genome_kmers[i].len() as f64 , genome_kmers[j].len() as f64);
                let ani = f64::powf(ratio, 1. / 31.);
                if ani > 0.990 {
                    used_genomes.insert(j);
                }
            }
        }
    }

    let ret_kmer_map: Vec<HashMap<u64, u32>> = genome_sketch
        .genome_kmer_counts
        .into_iter()
        .enumerate()
        .filter(|(i, _)| rep_genomes.contains(i))
        .map(|(_, x)| x)
        .collect();
    let ret_file_names = genome_sketch
        .file_names
        .into_iter()
        .enumerate()
        .filter(|(i, _)| rep_genomes.contains(i))
        .map(|(_, x)| x)
        .collect();

    info!(
        "Genomes before dereplication: {}. Genomes after dereplicaction: {}",
        num_gn,
        ret_kmer_map.len()
    );

    return GenomesSketch {
        genome_kmer_counts: ret_kmer_map,
        file_names: ret_file_names,
    };
}
