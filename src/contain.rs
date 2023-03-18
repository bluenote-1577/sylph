use crate::cmdline::*;
use statrs::distribution::{Poisson, Discrete, DiscreteCDF};
use crate::sketch::*;
use crate::types::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

pub fn contain(args: ContainArgs) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let (sequence_sketches, genome_sketches) = get_sketches(&args);
    let genome_index_vec = (0..genome_sketches.len()).collect::<Vec<usize>>();

    genome_index_vec.into_par_iter().for_each(|i| {
        let genome_sketch = &genome_sketches[i];
        let sequence_index_vec = (0..sequence_sketches.len()).collect::<Vec<usize>>();
        sequence_index_vec.into_par_iter().for_each(|j| {
            let sequence_sketch = &sequence_sketches[j];
            get_stats(&args, &genome_sketch, &sequence_sketch);
        });
    });

    log::info!("Finished contain.");
}

fn get_sketches(args: &ContainArgs) -> (Vec<SequencesSketch>, Vec<GenomeSketch>) {
    let genome_sketches: Mutex<Vec<GenomeSketch>> = Mutex::new(vec![]);
    let sequence_sketches: Mutex<Vec<SequencesSketch>> = Mutex::new(vec![]);

    let genome_index_vec = (0..args.genome_sketches.len()).collect::<Vec<usize>>();
    let sequence_index_vec = (0..args.sequence_sketches.len()).collect::<Vec<usize>>();

    sequence_index_vec.into_par_iter().for_each(|i| {
        let sequence_name = &args.sequence_sketches[i];
        let sequence_file = File::open(sequence_name).expect("Sequence sketch not a valid file");
        let seq_reader = BufReader::with_capacity(10_000_000, sequence_file);
        let sequence_sketch_enc: SequencesSketchEncode =
            bincode::deserialize_from(seq_reader).unwrap();
        let sequence_sketch = SequencesSketch::from_enc(sequence_sketch_enc);
        let mut locked = sequence_sketches.lock().unwrap();
        locked.push(sequence_sketch);
    });

    let sequence_sketches = sequence_sketches.into_inner().unwrap();
    if sequence_sketches.is_empty() {
        log::error!("No sequence sketches. Exiting");
        std::process::exit(1);
    }

    let c = sequence_sketches[0].c;
    let k = sequence_sketches[0].k;

    log::info!("Sequence sketch loading complete.");

    genome_index_vec.into_par_iter().for_each(|i| {
        let genome_name = &args.genome_sketches[i];
        let genome_sketch;
        if genome_name.contains(".prg") {
            let genome_file =
                File::open(&args.genome_sketches[i]).expect("Genome sketch {} not a valid file");
            let genome_reader = BufReader::with_capacity(10_000_000, genome_file);
            let genome_sketch_m: GenomeSketch = bincode::deserialize_from(genome_reader).unwrap();
            genome_sketch = genome_sketch_m;
        } else {
            let mut sketch_args = SketchArgs::default();
            sketch_args.c = c;
            sketch_args.k = k;
            genome_sketch = sketch_genome(c, k, &args.genome_sketches[i]).unwrap();
        }
        let mut locked = genome_sketches.lock().unwrap();
        locked.push(genome_sketch);
    });
    let genome_sketches = genome_sketches.into_inner().unwrap();
    log::info!("Genome sketch loading complete.");

    return (sequence_sketches, genome_sketches);
}

fn get_stats(args: &ContainArgs, genome_sketch: &GenomeSketch, sequence_sketch: &SequencesSketch) {
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
        return;
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
    for i in covs[covs.len()/2]..200{
        if pois.cdf(i.into()) < 0.99999{
            max_cov = i as f64;
        }
    }
    log::trace!("{:?}, {}", covs, max_cov);

    let mut full_covs = vec![0; gn_kmers.len() - contain_count];
    for cov in covs.into_iter(){
        if (cov as f64) < max_cov{
            full_covs.push(cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    log::trace!("mean cov {}", mean_cov);


    let (cov_cor_ani,lambda) = correct_ani_abund(&full_covs, sequence_sketch.k as f64);
    let ani_result = AniResult {
        ani,
        cov_cor_ani,
        seq_name: sequence_sketch.file_name.as_str(),
        gn_name: genome_sketch.file_name.as_str(),
        contig_name: genome_sketch.first_contig_name.as_str(),
        mean_cov,
        median_cov,
        containment_index: (contain_count, gn_kmers.len()),
        lambda
    };

    if ani_result.ani < 0.75 {
        return;
    }
    let (low_ani,high_ani,low_lambda,high_lambda) = bootstrap_interval(&full_covs, sequence_sketch.k as f64);
    let cov_cor_print;
    if let Some(cor_ani) = ani_result.cov_cor_ani{
        cov_cor_print = cor_ani * 100.;
    }
    else{
        cov_cor_print = -1.;
    }

    let lambda_print;
    if let Some(lambda) = ani_result.lambda{
        lambda_print = lambda;
    }
    else{
        lambda_print= -1.;
    }

    let ci_ani;
    if low_ani.is_none() || high_ani.is_none(){
        ci_ani = "NA-NA".to_string();
    }
    else{
        ci_ani = format!("{:.2}-{:.2}", low_ani.unwrap() * 100., high_ani.unwrap() * 100.);
    }

    let ci_lambda;
    if low_lambda.is_none() || high_lambda.is_none(){
        ci_lambda = "NA-NA".to_string();
    }
    else{
        ci_lambda = format!("{:.2}-{:.2}", low_lambda.unwrap(), high_lambda.unwrap());
    }

    println!(
        "{}\t{}\t{:.2}\t{:.2}\t{}\t{:.3}\t\t{}\t{:.2}\t{}\t{}/{}",
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
    )
}

fn correct_ani_abund(full_covs: &[u32], k: f64) -> (Option<f64>,Option<f64>) {
//    let thresh = 5;
//    let mut leq_threshold = 0;
//    for cov in covs.iter() {
//        if *cov < thresh {
//            leq_threshold += 1;
//        }
//    }

    let (mle_ani, mle_lambda) = mle_zip(full_covs, k);
//    if (leq_threshold as f64) < 0.9 * covs.len() as f64 {
//        return (mle_ani,mle_lambda);
//    }
    
    return (mle_ani, mle_lambda);
}

fn mle_zip(full_covs: &[u32], k: f64) -> (Option<f64>,Option<f64>) {

    let mut num_zero = 0;

    for x in full_covs{
        if *x == 0{
            num_zero += 1;
        }
    }

    let mean = mean(&full_covs).unwrap();
    let (lambda, pi) = newton_raphson((num_zero as f32 / full_covs.len() as f32).into(), mean.into());
//    log::trace!("lambda,pi {} {} {}", lambda,pi, num_zero as f64 / full_covs.len() as f64);
    let ani = f64::powf(1. - pi, 1./k);
    let ret_lambda;
    let ret_ani;
    if lambda < 0.{
        ret_lambda = None
    }
    else{
        ret_lambda = Some(lambda);
    }
    if ani < 0. || ani > 1.{
        ret_ani = None;
    }
    else{
        ret_ani = Some(ani);
    }
    return (ret_ani, ret_lambda);

}

fn newton_raphson(rat: f64, mean: f64) -> (f64,f64){
    let mut curr = mean / (1. - rat);
//    dbg!(1. - mean,rat);
    for _ in 0..1000{
        let t1 = (1. - rat) * curr;
        let t2 = mean * (1. - f64::powf(2.78281828, -curr));
        let t3 = 1. - rat;
        let t4 = mean * (f64::powf(2.78281828,-curr));
        curr = curr - (t1 - t2) / (t3 - t4);
    }
    return (curr, 1. - mean / curr);
}

fn mean(data: &[u32]) -> Option<f32> {
    let sum = data.iter().sum::<u32>() as f32;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f32),
        _ => None,
    }
}

fn bootstrap_interval(covs_full: &Vec<u32>, k: f64) -> (Option<f64>,Option<f64>,Option<f64>,Option<f64>){
    fastrand::seed(7);
    let num_samp = covs_full.len();
    let iters = 100;
    let mut res_ani = vec![];
    let mut res_lambda = vec![];

    for _ in 0..iters{
        let mut rand_vec = vec![];
        rand_vec.reserve(num_samp);
        for _ in 0..num_samp{
            rand_vec.push(covs_full[fastrand::usize(..covs_full.len())]);
        }
        let (ani,lambda) = mle_zip(&rand_vec, k);
        if ani.is_some() && lambda.is_some(){
            res_ani.push(ani);
            res_lambda.push(lambda);
        }
    }
    res_ani.sort_by(|x,y| x.partial_cmp(y).unwrap());
    res_lambda.sort_by(|x,y| x.partial_cmp(y).unwrap());
    if res_ani.len() < 20{
        return (None, None, None, None);
    }
    let suc = res_ani.len();
    let low_ani = res_ani[suc * 5 / 100 - 1];
    let high_ani = res_ani[suc * 95 / 100 - 1];
    let low_lambda = res_lambda[suc * 5 / 100 - 1];
    let high_lambda = res_lambda[suc * 95 / 100 - 1];

    return (low_ani, high_ani, low_lambda, high_lambda);
    

}
