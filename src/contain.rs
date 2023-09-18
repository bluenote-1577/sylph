use crate::cmdline::*;
use fxhash::FxHashMap;
use crate::constants::*;
use crate::inference::*;
use crate::sketch::*;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult, pseudotax: bool) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High {
        lambda_print = format!("HIGH");
    } else {
        lambda_print = format!("LOW");
    }
    let low_ani = ani_result.ani_ci.0;
    let high_ani = ani_result.ani_ci.1;
    let low_lambda = ani_result.lambda_ci.0;
    let high_lambda = ani_result.lambda_ci.1;

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none() {
        ci_ani = "NA-NA".to_string();
    } else {
        ci_ani = format!(
            "{:.2}-{:.2}",
            low_ani.unwrap() * 100.,
            high_ani.unwrap() * 100.
        );
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none() {
        ci_lambda = "NA-NA".to_string();
    } else {
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }

    if !pseudotax{
        println!(
            "{}\t{}\t{}\t{:.2}\t{}\t{:.3}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            print_final_ani,
            ani_result.naive_ani * 100.,
            ci_ani,
            ani_result.final_est_cov,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.contig_name,
        );
    }
    else{
    println!(
            "{}\t{}\t{:.4}\t{}\t{:.2}\t{}\t{:.3}\t{}\t{}\t{:.0}\t{:.3}\t{}/{}\t{}",
            ani_result.seq_name,
            ani_result.gn_name,
            ani_result.rel_abund.unwrap(),
            print_final_ani,
            ani_result.naive_ani * 100.,
            ci_ani,
            ani_result.final_est_cov,
            lambda_print,
            ci_lambda,
            ani_result.median_cov,
            ani_result.mean_cov,
            ani_result.containment_index.0,
            ani_result.containment_index.1,
            ani_result.contig_name,
        );

    }
}

pub fn contain(args: ContainArgs) {
    let level;
    if args.trace {
        level = log::LevelFilter::Trace;
    } else {
        level = log::LevelFilter::Info;
    }
    simple_logger::SimpleLogger::new()
        .with_level(level)
        .init()
        .unwrap();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    log::info!("Obtaining sketches...");
    let mut genome_sketch_files = vec![];
    let mut genome_files = vec![];
    let mut read_sketch_files = vec![];
    let mut read_files = vec![];

    for file in args.files.iter() {
        if file.ends_with(QUERY_FILE_SUFFIX) {
            genome_sketch_files.push(file);
        } else if file.ends_with(SAMPLE_FILE_SUFFIX) {
            read_sketch_files.push(file);
        } else if is_fasta(&file) {
            genome_files.push(file);
        } else if is_fastq(&file) {
            read_files.push(file);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    let genome_sketches = get_genome_sketches(&args, &genome_sketch_files, &genome_files);
    let genome_index_vec = (0..genome_sketches.len()).collect::<Vec<usize>>();
    log::info!("Finished obtaining genome sketches.");

    
    if genome_sketches.is_empty() {
        log::error!("No genome sketches found; see sylph contain -h for help. Exiting");
        std::process::exit(1);
    }

    if genome_sketches.first().unwrap().pseudotax_tracked_nonused_kmers.is_none() && args.pseudotax{
        log::error!("--pseudotax is enabled, but *.sylsample was not sketched with the --enable-pseudotax option. Exiting");
        std::process::exit(1);
    }

    read_files.extend(read_sketch_files.clone());
    let write_mutex = Mutex::new(true);
    let first_write = Mutex::new(true);
    let sequence_index_vec = (0..read_files.len()).collect::<Vec<usize>>();


    sequence_index_vec.into_iter().for_each(|j| {
        let is_sketch = j >= read_files.len() - read_sketch_files.len();
        let sequence_sketch = get_seq_sketch(&args, read_files[j], is_sketch, genome_sketches[0].c, genome_sketches[0].k);
        if sequence_sketch.is_some(){
            let sequence_sketch = sequence_sketch.unwrap();
            let stats_vec_seq: Mutex<Vec<AniResult>> = Mutex::new(vec![]);
            genome_index_vec.par_iter().for_each(|i| {
                let genome_sketch = &genome_sketches[*i];
                let mut res = get_stats(&args, &genome_sketch, &sequence_sketch, None);
                if res.is_some() {
                    res.as_mut().unwrap().genome_sketch_index = *i;
                    stats_vec_seq.lock().unwrap().push(res.unwrap());
                    
                }
            });
            let mut stats_vec_seq = stats_vec_seq.into_inner().unwrap();

            if args.pseudotax{
                log::info!("Pseudotax enabled. Reassigning k-mers for {} genomes...", stats_vec_seq.len());
                let winner_map = winner_table(&stats_vec_seq, &genome_sketches);
            //If pseudotax, get k-mer to genome table: table <k_mer, &genome_sketch> = table(results)
                let remaining_genomes = stats_vec_seq.iter().map(|x| &genome_sketches[x.genome_sketch_index]).collect::<Vec<&GenomeSketch>>();
                let stats_vec_seq_2 = Mutex::new(vec![]);
                remaining_genomes.into_par_iter().for_each(|genome_sketch|{
                    let res = get_stats(&args, &genome_sketch, &sequence_sketch, Some(&winner_map));
                    if res.is_some() {
                        stats_vec_seq_2.lock().unwrap().push(res.unwrap());
                    }
                });
                stats_vec_seq = stats_vec_seq_2.into_inner().unwrap();
                log::info!("{} genomes passing reassigned k-mer threshold. ", stats_vec_seq.len());

                let total_cov = stats_vec_seq.iter().map(|x| x.final_est_cov).sum::<f64>();
                for thing in stats_vec_seq.iter_mut(){
                    thing.rel_abund = Some(thing.final_est_cov/total_cov * 100.);
                }
            
            //for loop over genomes in results
            //Reassign k-mers to genomes: get_stats(table)
            }

            if args.pseudotax{
                stats_vec_seq.sort_by(|x,y| y.rel_abund.unwrap().partial_cmp(&x.rel_abund.unwrap()).unwrap());
            }
            else{
                stats_vec_seq.sort_by(|x,y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
            }

            if *first_write.lock().unwrap(){
                *first_write.lock().unwrap() = false;
                print_header(args.pseudotax);

            }
            let _tmp = write_mutex.lock().unwrap();
            for res in stats_vec_seq{
                print_ani_result(&res, args.pseudotax);
            }
        }
        log::info!("Finished contain for {}.", &read_files[j]);
    });

    log::info!("Finished contain.");
}

fn winner_table<'a>(results : &Vec<AniResult>, genome_sketches: &'a Vec<GenomeSketch>) -> FxHashMap<Kmer, &'a GenomeSketch> {
    let mut kmer_to_genome_map = FxHashMap::default();
    let mut return_map = FxHashMap::default();
    for res in results.iter(){
        let gn_sketch = &genome_sketches[res.genome_sketch_index];
        for kmer in gn_sketch.genome_kmers.iter(){
            let v = kmer_to_genome_map.entry(*kmer).or_insert(vec![]);
            v.push((res.final_est_ani, gn_sketch));
        }
        
        if gn_sketch.pseudotax_tracked_nonused_kmers.is_some(){
            for kmer in gn_sketch.pseudotax_tracked_nonused_kmers.as_ref().unwrap().iter(){
                let v = kmer_to_genome_map.entry(*kmer).or_insert(vec![]);
                v.push((res.final_est_ani, gn_sketch));
            }
        }
    }

    for (kmer, list) in kmer_to_genome_map.iter_mut(){
        list.sort_by(|a,b| b.partial_cmp(&a).unwrap());
        return_map.insert(*kmer, list.first().unwrap().1);
    }

    return return_map;
}

fn print_header(pseudotax: bool) {
    if !pseudotax{
    println!(
            "Sample_file\tQuery_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            );
    }
    else{
    println!(
            "Sample_file\tQuery_file\tRelative_abundance\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov_geq1\tContainment_ind\tContig_name",
            );
    }
}

fn get_genome_sketches(
    args: &ContainArgs,
    genome_sketch_files: &Vec<&String>,
    genome_files: &Vec<&String>,
) -> Vec<GenomeSketch> {
    let mut lowest_genome_c = None;
    let mut current_k = None;

    let genome_sketches = Mutex::new(vec![]);

    for genome_sketch_file in genome_sketch_files {
        let file = File::open(genome_sketch_file).expect(&format!("The sketch `{}` not a valid file ", genome_sketch_file));
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch ",
                &genome_sketch_file
            ));
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    return genome_sketches.into_inner().unwrap();
}

fn get_seq_sketch(
    args: &ContainArgs,
    read_file: &str,
    is_sketch_file: bool,
    genome_c: usize,
    genome_k: usize,
) -> Option<SequencesSketch> {
    if is_sketch_file {
        let read_sketch_file = read_file;
        let file = File::open(read_sketch_file.clone()).expect(&format!(
            "The sketch `{}` not a valid file ",
            &read_sketch_file
        ));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch_enc: SequencesSketchEncode = bincode::deserialize_from(read_reader).expect(
            &format!("The sketch `{}` is not a valid sketch ", read_sketch_file),
        );
        let read_sketch = SequencesSketch::from_enc(read_sketch_enc);
        if read_sketch.c > genome_c {
            error!("{} value of -c for {} is {} -- greater than the smallest value of -c for a genome sketch {}. Exiting.", read_file, read_sketch.c, read_sketch_file, genome_c);
            return None;
        }

        return Some(read_sketch);
    } else {
        if genome_c < args.c {
            error!("{} error: value of -c for contain = {} -- greater than the smallest value of -c for a genome sketch = {}. Continuing without sketching.", read_file, args.c, genome_c);
            return None;
        } else if genome_k != args.k {
            error!(
                "{} -k {} is not equal to -k {} found in sketches. Continuing without sketching.",
                read_file, args.k, genome_k
            );
            return None;
        } else {
            let read_sketch_opt = sketch_query(args.c, args.k, args.threads, &read_file);
            return read_sketch_opt;
        }
    }
}

fn _get_sketches_rewrite(args: &ContainArgs) -> (Vec<SequencesSketch>, Vec<GenomeSketch>) {
    let mut read_sketch_files = vec![];
    let mut genome_sketch_files = vec![];
    let mut read_files = vec![];
    let mut genome_files = vec![];
    for file in args.files.iter() {
        if file.ends_with(QUERY_FILE_SUFFIX) {
            genome_sketch_files.push(file);
        } else if file.ends_with(SAMPLE_FILE_SUFFIX) {
            read_sketch_files.push(file);
        } else if is_fasta(&file) {
            genome_files.push(file);
        } else if is_fastq(&file) {
            read_files.push(file);
        } else {
            warn!(
                "{} file extension is not a sketch or a fasta/fastq file.",
                &file
            );
        }
    }

    let genome_sketches = Mutex::new(vec![]);
    let read_sketches = Mutex::new(vec![]);
    //read c can be lower than lowest genome c.
    let mut lowest_genome_c = None;
    let mut current_k = None;

    read_sketch_files.into_par_iter().for_each(|read_sketch_file|{
        let file = File::open(read_sketch_file.clone()).expect(&format!("The sketch `{}` not a valid file ", &read_sketch_file));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch_enc: SequencesSketchEncode = bincode::deserialize_from(read_reader).expect(&format!(
            "The sketch `{}` is not a valid sketch ",
            read_sketch_file
        ));
        let read_sketch = SequencesSketch::from_enc(read_sketch_enc);
        if lowest_genome_c.is_some() && read_sketch.c > lowest_genome_c.unwrap(){
            error!("Value of -c for {} is {} -- greater than the smallest value of -c for a genome sketch {}. Exiting.", read_sketch.c, read_sketch_file, lowest_genome_c.unwrap());
            std::process::exit(1);
        }
        read_sketches.lock().unwrap().push(read_sketch);
    });

    for genome_sketch_file in genome_sketch_files {
        let file =
            File::open(genome_sketch_file.clone()).expect(&format!("The sketch `{}` not a valid file ", genome_sketch_file));
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader)
            .expect(&format!(
                "The sketch `{}` is not a valid sketch ",
                &genome_sketch_file
            ));
        if genome_sketches_vec.is_empty() {
            continue;
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none() {
            lowest_genome_c = Some(c);
        } else if lowest_genome_c.unwrap() < c {
            lowest_genome_c = Some(c);
        }
        if current_k.is_none() {
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        } else if current_k.unwrap() != k {
            error!("Query sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer, args.pseudotax);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer, args.pseudotax);
                if genome_sketch_opt.is_some() {
                    genome_sketches.lock().unwrap().push(genome_sketch_opt.unwrap());
                }
            }
        }
    });

    read_files.into_par_iter().for_each(|read_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            let read_sketch_opt = sketch_query(args.c, args.k, args.threads, &read_file);
            if read_sketch_opt.is_some() {
                read_sketches.lock().unwrap().push(read_sketch_opt.unwrap());
            }
        }
    });

    return (
        read_sketches.into_inner().unwrap(),
        genome_sketches.into_inner().unwrap(),
    );
}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &SequencesSketch,
    winner_map: Option<&FxHashMap<Kmer, & GenomeSketch>>
) -> Option<AniResult<'a>> {
    if genome_sketch.k != sequence_sketch.k {
        log::error!(
            "k parameter for reads {} != k parameter for genome {}",
            sequence_sketch.k,
            genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c < sequence_sketch.c {
        log::error!(
            "c parameter for reads {} > c parameter for genome {}",
            sequence_sketch.c,
            genome_sketch.c
        );
        std::process::exit(1);
    }
    let mut contain_count = 0;
    let mut covs = vec![];
    let gn_kmers = &genome_sketch.genome_kmers;

    //let start_t_initial = Instant::now();
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            if winner_map.is_some(){
                if winner_map.unwrap()[kmer] != genome_sketch{
                    continue
                }
            }
            contain_count += 1;
            covs.push(sequence_sketch.kmer_counts[kmer]);
        }
    }
    //log::trace!("Hashing time {:?}", Instant::now() - start_t_initial);
    //let start_t_initial = Instant::now();
    if covs.is_empty() {
        return None;
    }
    let naive_ani = f64::powf(
        contain_count as f64 / gn_kmers.len() as f64,
        1. / genome_sketch.k as f64,
    );
    covs.sort();
    //let covs = &covs[0..covs.len() * 99 / 100];
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov).unwrap();
    let mut max_cov = f64::MAX;
    if median_cov < 30.{
        for i in covs.len() / 2..covs.len(){
            let cov = covs[i];
            if pois.cdf(cov.into()) < CUTOFF_PVALUE {
                max_cov = cov as f64;
            } else {
                break;
            }
        }
    }
    log::trace!("COV VECTOR for {}/{}: {:?}, {}", sequence_sketch.file_name, genome_sketch.file_name ,covs, max_cov);

    let mut full_covs = vec![0; gn_kmers.len() - contain_count];
    for cov in covs.iter() {
        if (*cov as f64) < max_cov {
            full_covs.push(*cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    let geq1_mean_cov = full_covs.iter().sum::<u32>() as f64 / covs.len() as f64;

    let use_lambda;
    if median_cov > MEDIAN_ANI_THRESHOLD {
        use_lambda = AdjustStatus::High
    } else {
        let test_lambda;
        if args.ratio {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        } else if args.mme {
            test_lambda = mme_lambda(&full_covs)
        } else if args.nb {
            test_lambda = binary_search_lambda(&full_covs)
        } else if args.mle {
            test_lambda = mle_zip(&full_covs, sequence_sketch.k as f64)
        } else {
            test_lambda = ratio_lambda(&full_covs, args.min_count_correct)
        };
        if test_lambda.is_none() {
            use_lambda = AdjustStatus::Low
        } else {
            use_lambda = AdjustStatus::Lambda(test_lambda.unwrap());
        }
    }

    let final_est_cov;

    if let AdjustStatus::Lambda(lam) = use_lambda {
        final_est_cov = lam
    } else {
        final_est_cov = median_cov;
    }

    let opt_lambda;
    if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {
        opt_lambda = None
    } else {
        opt_lambda = Some(final_est_cov)
    };

    let opt_est_ani = ani_from_lambda(opt_lambda, mean_cov, sequence_sketch.k as f64, &full_covs);
    
    let final_est_ani;
    if opt_lambda.is_none() || opt_est_ani.is_none() || args.no_adj {
        final_est_ani = naive_ani;
    } else {
        final_est_ani = opt_est_ani.unwrap();
    }

    let min_ani = if args.minimum_ani.is_some() {args.minimum_ani.unwrap()/100. }
        else if args.pseudotax { MIN_ANI_P_DEF } 
        else { MIN_ANI_DEF };
    if final_est_ani < min_ani {
        return None;
    }

    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if !args.no_ci && opt_lambda.is_some() {
        let bootstrap = bootstrap_interval(&full_covs, sequence_sketch.k as f64, &args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    }

    

    let ani_result = AniResult {
        naive_ani,
        final_est_ani,
        final_est_cov,
        seq_name: sequence_sketch.file_name.clone(),
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov: geq1_mean_cov,
        median_cov,
        containment_index: (contain_count, gn_kmers.len()),
        lambda: use_lambda,
        ani_ci: (low_ani, high_ani),
        lambda_ci: (low_lambda, high_lambda),
        genome_sketch_index: usize::MAX,
        rel_abund: None
    };
    //log::trace!("Other time {:?}", Instant::now() - start_t_initial);

    return Some(ani_result);
}

fn _ani_from_lambda_moment(lambda: Option<f64>, mean: f64, k: f64) -> Option<f64> {
    if lambda.is_none() {
        return None;
    }
    let lambda = lambda.unwrap();
    let pi = ((lambda + 1.) * mean - mean * mean - mean) / ((lambda + 1.) * mean - mean);
    let ret_ani;
    let ani = f64::powf(1. - pi, 1. / k);
    if ani < 0. || ani.is_nan() {
        ret_ani = None;
    } else {
        if ani > 1. {
            ret_ani = Some(1.)
        } else {
            ret_ani = Some(ani);
        }
    }
    return ret_ani;
}

fn ani_from_lambda(lambda: Option<f64>, _mean: f64, k: f64, full_cov: &[u32]) -> Option<f64> {
    if lambda.is_none() {
        return None;
    }
    let mut contain_count = 0;
    let mut _zero_count = 0;
    for x in full_cov {
        if *x != 0 {
            contain_count += 1;
        } else {
            _zero_count += 1;
        }
    }

    let lambda = lambda.unwrap();
    let adj_index =
        contain_count as f64 / (1. - f64::powf(2.78281828, -lambda)) / full_cov.len() as f64;
    let ret_ani;
    //let ani = f64::powf(1. - pi, 1./k);
    let ani = f64::powf(adj_index, 1. / k);
    if ani < 0. || ani.is_nan() {
        ret_ani = None;
    } else {
        if ani > 1. {
            ret_ani = Some(ani)
        } else {
            ret_ani = Some(ani);
        }
    }
    return ret_ani;
}

fn mle_zip(full_covs: &[u32], _k: f64) -> Option<f64> {
    let mut num_zero = 0;
    let mut count_set: HashSet<_> = HashSet::default();

    for x in full_covs {
        if *x == 0 {
            num_zero += 1;
        } else {
            count_set.insert(x);
        }
    }

    //Lack of information for inference, retun None.
    if count_set.len() == 1 {
        return None;
    }

    if full_covs.len() - num_zero < SAMPLE_SIZE_CUTOFF {
        return None;
    }

    let mean = mean(&full_covs).unwrap();
    let lambda = newton_raphson(
        (num_zero as f32 / full_covs.len() as f32).into(),
        mean.into(),
    );
    //    log::trace!("lambda,pi {} {} {}", lambda,pi, num_zero as f64 / full_covs.len() as f64);
    let ret_lambda;
    if lambda < 0. || lambda.is_nan() {
        ret_lambda = None
    } else {
        ret_lambda = Some(lambda);
    }

    return ret_lambda;
}

fn newton_raphson(rat: f64, mean: f64) -> f64 {
    let mut curr = mean / (1. - rat);
    //    dbg!(1. - mean,rat);
    for _ in 0..1000 {
        let t1 = (1. - rat) * curr;
        let t2 = mean * (1. - f64::powf(2.78281828, -curr));
        let t3 = 1. - rat;
        let t4 = mean * (f64::powf(2.78281828, -curr));
        curr = curr - (t1 - t2) / (t3 - t4);
    }
    return curr;
}

fn var(data: &[u32]) -> Option<f32> {
    if data.is_empty() {
        return None;
    }
    let mean = mean(data).unwrap();
    let mut var = 0.;
    for x in data {
        var += (*x as f32 - mean) * (*x as f32 - mean)
    }
    return Some(var / data.len() as f32);
}

fn mean(data: &[u32]) -> Option<f32> {
    let sum = data.iter().sum::<u32>() as f32;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f32),
        _ => None,
    }
}

fn bootstrap_interval(
    covs_full: &Vec<u32>,
    k: f64,
    args: &ContainArgs,
) -> (Option<f64>, Option<f64>, Option<f64>, Option<f64>) {
    fastrand::seed(7);
    let num_samp = covs_full.len();
    let iters = 100;
    let mut res_ani = vec![];
    let mut res_lambda = vec![];

    for _ in 0..iters {
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp {
            rand_vec.push(covs_full[fastrand::usize(..covs_full.len())]);
        }
        let lambda;
        if args.ratio {
            lambda = ratio_lambda(&rand_vec, args.min_count_correct);
        } else if args.mme {
            lambda = mme_lambda(&rand_vec);
        } else if args.nb {
            lambda = binary_search_lambda(&rand_vec);
        } else if args.mle {
            lambda = mle_zip(&rand_vec, k);
        } else {
            lambda = ratio_lambda(&rand_vec,args.min_count_correct);
        }
        let ani = ani_from_lambda(lambda, mean(&rand_vec).unwrap().into(), k, &rand_vec);
        if ani.is_some() && lambda.is_some() {
            if !ani.unwrap().is_nan() && !lambda.unwrap().is_nan() {
                res_ani.push(ani);
                res_lambda.push(lambda);
            }
        }
    }
    res_ani.sort_by(|x, y| x.partial_cmp(y).unwrap());
    res_lambda.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if res_ani.len() < 50 {
        return (None, None, None, None);
    }
    let suc = res_ani.len();
    let low_ani = res_ani[suc * 5 / 100 - 1];
    let high_ani = res_ani[suc * 95 / 100 - 1];
    let low_lambda = res_lambda[suc * 5 / 100 - 1];
    let high_lambda = res_lambda[suc * 95 / 100 - 1];

    return (low_ani, high_ani, low_lambda, high_lambda);
}

fn ratio_lambda(full_covs: &Vec<u32>, min_count_correct: f64) -> Option<f64> {
    let mut num_zero = 0;
    let mut count_map: FxHashMap<_, _> = FxHashMap::default();

    for x in full_covs {
        if *x == 0 {
            num_zero += 1;
        } else {
            let c = count_map.entry(*x as usize).or_insert(0);
            *c += 1;
        }
    }

    //Lack of information for inference, retun None.
    if count_map.len() == 1 {
        return None;
    }

    if full_covs.len() - num_zero < SAMPLE_SIZE_CUTOFF {
        return None;
    } else {
        let mut sort_vec: Vec<(_, _)> = count_map.iter().map(|x| (x.1, x.0)).collect();
        sort_vec.sort_by(|x, y| y.cmp(&x));
        let most_ind = sort_vec[0].1;
        if !count_map.contains_key(&(most_ind + 1)) {
            return None;
        }
        let count_p1 = count_map[&(most_ind + 1)] as f64;
        let count = count_map[&most_ind] as f64;
        if count_p1 < min_count_correct || count < min_count_correct{
            return None;
        }
        let lambda = Some(count_p1 / count * ((most_ind + 1) as f64));
        return lambda;
    }
}

fn mme_lambda(full_covs: &[u32]) -> Option<f64> {
    let mut num_zero = 0;
    let mut count_set: HashSet<_> = HashSet::default();

    for x in full_covs {
        if *x == 0 {
            num_zero += 1;
        } else {
            count_set.insert(x);
        }
    }

    //Lack of information for inference, retun None.
    if count_set.len() == 1 {
        return None;
    }

    if full_covs.len() - num_zero < SAMPLE_SIZE_CUTOFF {
        return None;
    }

    let mean = mean(&full_covs).unwrap();
    let var = var(&full_covs).unwrap();
    let lambda = var / mean + mean - 1.;
    if lambda < 0. {
        return None;
    } else {
        return Some(lambda as f64);
    }
}
