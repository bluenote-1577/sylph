use crate::cmdline::*;
use crate::constants::*;
use crate::inference::*;
use crate::sketch::*;
use crate::types::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::fs;
use std::io::BufReader;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult) {
    let cov_cor_print;
    if let Some(cor_ani) = ani_result.cov_cor_ani {
        cov_cor_print = format!("{:.2}", cor_ani * 100.);
    } else {
        cov_cor_print = format!("NA");
    }

    let lambda_print;
    if let Some(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else {
        lambda_print = format!("NA");
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

    println!(
        "{}\t{}\t{:.2}\t{}\t{}\t{}\t\t{}\t{:.2}\t{}\t{}/{}",
        ani_result.seq_name,
        ani_result.gn_name,
        ani_result.ani * 100.,
        cov_cor_print,
        ci_ani,
        lambda_print,
        ci_lambda,
        ani_result.median_cov,
        ani_result.contig_name,
        ani_result.containment_index.0,
        ani_result.containment_index.1,
    );
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

    let (sequence_sketches, genome_sketches) = get_sketches(&args);
    let genome_index_vec = (0..genome_sketches.len()).collect::<Vec<usize>>();
    let stats_vec: Mutex<Vec<AniResult>> = Mutex::new(vec![]);

    genome_index_vec.into_par_iter().for_each(|i| {
        let genome_sketch = &genome_sketches[i];
        let sequence_index_vec = (0..sequence_sketches.len()).collect::<Vec<usize>>();
        sequence_index_vec.into_par_iter().for_each(|j| {
            let sequence_sketch = &sequence_sketches[j];
            let res = get_stats(&args, &genome_sketch, &sequence_sketch);
            if res.is_some() {
                stats_vec.lock().unwrap().push(res.unwrap());
            }
        });
    });

    let mut result_vec = stats_vec.into_inner().unwrap();
    result_vec.sort_by(|x, y| y.cov_cor_ani.partial_cmp(&x.cov_cor_ani).unwrap());
    for res in result_vec {
        print_ani_result(&res);
    }

    log::info!("Finished contain.");
}

fn get_sketches(args: &ContainArgs) -> (Vec<SequencesSketch>, Vec<GenomeSketch>) {
    let genome_sketches: Mutex<Vec<GenomeSketch>> = Mutex::new(vec![]);
    let sequence_sketches: Mutex<Vec<SequencesSketch>> = Mutex::new(vec![]);
    let sequence_sketches_clone;
    if args.sequence_sketches.is_some(){
        sequence_sketches_clone = args.sequence_sketches.as_ref().unwrap().clone();
    }
    else if args.sequence_folder.is_some(){
        let paths = fs::read_dir(args.sequence_folder.as_ref().unwrap()).unwrap();
        sequence_sketches_clone = paths.into_iter().map(|x| x.unwrap().path().display().to_string()).collect();
    }
    else{
        panic!("Input sequence sketch argument not valid");
    }
    //    let paired_sequences = collect_pairs_file_names(sequence_sketches_clone);
    //    for pair in paired_sequences.iter(){
    //        if pair.len() > 1{
    //            if pair[0].contains(".prs"){
    //                log::warn!("{:?} appear to be paired sequence sketches and will be merged.", pair)
    //            }
    //        }
    //    }

    let sketch_c = Mutex::new(None);
    let sketch_k = Mutex::new(None);
    let sketch_name_track = Mutex::new(None);

    let sequence_index_vec =
        (0..sequence_sketches_clone.len()).collect::<Vec<usize>>();
    sequence_index_vec.into_par_iter().for_each(|i| {
        let sequence_file_name = &sequence_sketches_clone[i];
        let sequence_sketch;
        if sequence_file_name.contains(".prs") {
            let sequence_file = File::open(sequence_file_name).expect("Sequence sketch not a valid file");
            let seq_reader = BufReader::with_capacity(10_000_000, sequence_file);
            let sequence_sketch_enc: SequencesSketchEncode =
                bincode::deserialize_from(seq_reader).unwrap();
            {
                let mut s_c = sketch_c.lock().unwrap();
                let mut s_k = sketch_k.lock().unwrap();
                let mut s_n = sketch_name_track.lock().unwrap();
                if s_c.is_none() {
                    *s_c = Some(sequence_sketch_enc.c);
                    *s_n = Some(sequence_sketch_enc.file_name.clone());
                } else {
                    if s_c.unwrap() != sequence_sketch_enc.c {
                        log::error!(
                            "Value of 'c' for {} is not equal to value of 'c' for {}. Exiting.",
                            sequence_file_name,
                            s_n.as_ref().unwrap()
                        );
                    }
                }
                if s_k.is_none() {
                    *s_k = Some(sequence_sketch_enc.k);
                } else {
                    if s_k.unwrap() != sequence_sketch_enc.k {
                        log::error!(
                            "Value of 'k' for {} is not equal to value of 'k' for {}. Exiting.",
                            sequence_file_name,
                            s_n.as_ref().unwrap()
                        );
                    }
                }
            }
            sequence_sketch = SequencesSketch::from_enc(sequence_sketch_enc);
        } else {
            {
                let s_c = sketch_c.lock().unwrap();
                let s_k = sketch_k.lock().unwrap();
                let s_n = sketch_name_track.lock().unwrap();
                if s_c.is_some() {
                    if s_c.unwrap() != args.c {
                        log::error!(
                            "Value of 'c' for {} is not equal to value of 'c' for {}. Exiting.",
                            sequence_file_name,
                            s_n.as_ref().unwrap()
                        );
                    }
                    if s_k.unwrap() != args.k {
                        log::error!(
                            "Value of 'k' for {} is not equal to value of 'k' for {}. Exiting.",
                            sequence_file_name,
                            s_n.as_ref().unwrap()
                        );
                    }
                }
            }
            sequence_sketch = sketch_query(args.c, args.k, args.threads, sequence_file_name);
        }
        let mut locked = sequence_sketches.lock().unwrap();
        locked.push(sequence_sketch);
    });

    log::info!("Sequence sketch loading complete.");

    let sequence_sketches = sequence_sketches.into_inner().unwrap();
    if sequence_sketches.is_empty() {
        log::error!("No sequence sketches. Exiting");
        std::process::exit(1);
    }

    let c = sequence_sketches[0].c;
    let k = sequence_sketches[0].k;

    let genome_sketches_clone;
    if args.genome_sketches.is_some(){
        genome_sketches_clone = args.genome_sketches.as_ref().unwrap().clone();
    }
    else if args.genome_folder.is_some(){
        let paths = fs::read_dir(args.genome_folder.as_ref().unwrap()).unwrap();
        genome_sketches_clone = paths.into_iter().map(|x| x.unwrap().path().display().to_string()).collect();
        
    }
    else{
        panic!("Input sequence sketch argument not valid");
    }

    let genome_index_vec = (0..genome_sketches_clone.len()).collect::<Vec<usize>>();
    genome_index_vec.into_par_iter().for_each(|i| {
        let genome_name = &genome_sketches_clone[i];
        let genome_sketch;
        if genome_name.contains(".prg") {
            let genome_file =
                File::open(&genome_sketches_clone[i]).expect("Genome sketch {} not a valid file");
            let genome_reader = BufReader::with_capacity(10_000_000, genome_file);
            let genome_sketch_m: GenomeSketch =
                bincode::deserialize_from(genome_reader).expect(&format!(
                    "Genome sketch {} is not a valid sketch.",
                    &genome_sketches_clone[i]
                ));
            genome_sketch = genome_sketch_m;
        } else {
            let mut sketch_args = SketchArgs::default();
            sketch_args.c = c;
            sketch_args.k = k;
            genome_sketch = sketch_genome(c, k, &genome_sketches_clone[i]).unwrap();
        }
        let mut locked = genome_sketches.lock().unwrap();
        locked.push(genome_sketch);
    });
    let genome_sketches = genome_sketches.into_inner().unwrap();
    log::info!("Genome sketch loading complete.");

    return (sequence_sketches, genome_sketches);
}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &'a SequencesSketch,
) -> Option<AniResult<'a>> {
    if genome_sketch.k != sequence_sketch.k {
        log::error!(
            "k parameter for query {} != k parameter for reference {}",
            sequence_sketch.k,
            genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c != sequence_sketch.c {
        log::error!(
            "c parameter for query {} != c parameter for reference {}",
            sequence_sketch.c,
            genome_sketch.c
        );
        std::process::exit(1);
    }
    let mut contain_count = 0;
    let mut covs = vec![];
    let gn_kmers = &genome_sketch.genome_kmers;
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            contain_count += 1;
            covs.push(sequence_sketch.kmer_counts[kmer]);
        }
    }
    if covs.is_empty() {
        return None;
    }
    let ani = f64::powf(
        contain_count as f64 / gn_kmers.len() as f64,
        1. / genome_sketch.k as f64,
    );
    covs.sort();
    //let covs = &covs[0..covs.len() * 99 / 100];
    let median_cov = covs[covs.len() / 2] as f64;
    let pois = Poisson::new(median_cov).unwrap();
    let mut max_cov = median_cov;
    for i in covs[covs.len() / 2]..1000000 {
        if pois.cdf(i.into()) < CUTOFF_PVALUE {
            max_cov = i as f64;
        } else {
            break;
        }
    }
    log::trace!("{:?}, {}", covs, max_cov);

    let mut full_covs = vec![0; gn_kmers.len() - contain_count];
    for cov in covs.into_iter() {
        if (cov as f64) < max_cov {
            full_covs.push(cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    log::trace!("mean cov {}", mean_cov);

    let lambda = mle_zip(&full_covs, sequence_sketch.k as f64);
    let skmer_lambda = ratio_lambda(&full_covs);
    let mme_lambda = mme_lambda(&full_covs);
    //    let nb_lambda = binary_search_lambda(&full_covs);
    //    if skmer_lambda.is_none() && lambda.is_some(){
    //eprintln!("l {:?} skl {:?}", lambda, skmer_lambda);
    //    }

    let use_lambda = if args.ratio {
        skmer_lambda
    } else if args.mme {
        mme_lambda
    } else if args.nb {
        binary_search_lambda(&full_covs)
    } else {
        lambda
    };
    let cov_cor_ani = ani_from_lambda(use_lambda, mean_cov, sequence_sketch.k as f64, &full_covs);
    if ani < 0.75 {
        return None;
    }
    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if args.ci {
        let bootstrap = bootstrap_interval(&full_covs, sequence_sketch.k as f64, &args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    } else {
    }

    let ani_result = AniResult {
        ani,
        cov_cor_ani,
        seq_name: sequence_sketch.file_name.as_str(),
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov,
        median_cov,
        containment_index: (contain_count, gn_kmers.len()),
        lambda: use_lambda,
        ani_ci: (low_ani, high_ani),
        lambda_ci: (low_lambda, high_lambda),
    };

    return Some(ani_result);
}

fn ani_from_lambda_moment(lambda: Option<f64>, mean: f64, k: f64) -> Option<f64> {
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

fn ani_from_lambda(lambda: Option<f64>, mean: f64, k: f64, full_cov: &[u32]) -> Option<f64> {
    if lambda.is_none() {
        return None;
    }
    let mut contain_count = 0;
    let mut zero_count = 0;
    for x in full_cov {
        if *x != 0 {
            contain_count += 1;
        } else {
            zero_count += 1;
        }
    }

    let lambda = lambda.unwrap();
    let pi = 1. - mean / lambda;
    let adj_index =
        contain_count as f64 / (1. - f64::powf(2.78281828, -lambda)) / full_cov.len() as f64;
    let ret_ani;
    //let ani = f64::powf(1. - pi, 1./k);
    let ani = f64::powf(adj_index, 1. / k);
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

fn mle_zip(full_covs: &[u32], k: f64) -> Option<f64> {
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

    return (ret_lambda);
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
            lambda = ratio_lambda(&rand_vec);
        } else if args.mme {
            lambda = mme_lambda(&rand_vec);
        } else if args.nb {
            lambda = binary_search_lambda(&rand_vec);
        } else {
            lambda = mle_zip(&rand_vec, k);
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

fn ratio_lambda(full_covs: &Vec<u32>) -> Option<f64> {
    let mut num_zero = 0;
    let mut count_map: HashMap<_, _> = HashMap::default();

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
        let lambda = Some(
            (count_map[&(most_ind + 1)] as f64 / count_map[&most_ind] as f64)
                * ((most_ind + 1) as f64),
        );
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
