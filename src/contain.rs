use crate::cmdline::*;
use log::*;
use crate::constants::*;
use crate::inference::*;
use crate::sketch::*;
use crate::types::*;
use rayon::prelude::*;
use statrs::distribution::{DiscreteCDF, Poisson};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;

fn print_ani_result(ani_result: &AniResult) {
    let print_final_ani = format!("{:.2}", f64::min(ani_result.final_est_ani * 100., 100.));
    let lambda_print;
    if let AdjustStatus::Lambda(lambda) = ani_result.lambda {
        lambda_print = format!("{:.3}", lambda);
    } else if ani_result.lambda == AdjustStatus::High{
        lambda_print = format!("HIGH");
    }
    else{
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

    let (sequence_sketches, genome_sketches) = get_sketches_rewrite(&args);
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
    result_vec.sort_by(|x, y| y.final_est_ani.partial_cmp(&x.final_est_ani).unwrap());
    println!(
        "Read_file\tGenome_file\tAdjusted_ANI\tNaive_ANI\tANI_5-95_percentile\tEff_cov\tEff_lambda\tLambda_5-95_percentile\tMedian_cov\tMean_cov\tContainment_ind\tContig_name",
        );
    for res in result_vec {
        print_ani_result(&res);
    }

    log::info!("Finished contain.");
}


fn get_sketches_rewrite(args: &ContainArgs) -> (Vec<SequencesSketch>, Vec<GenomeSketch>) {
    let mut read_sketch_files = vec![];
    let mut genome_sketch_files = vec![];
    let mut read_files = vec![];
    let mut genome_files = vec![];
    for file in args.files.iter(){
        if file.ends_with(".pgs"){
            genome_sketch_files.push(file);
        }
        else if file.ends_with(".prs"){
            read_sketch_files.push(file);
        }
        else if is_fasta(&file){
            genome_files.push(file);
        }
        else if is_fastq(&file){
            read_files.push(file);
        }
        else{
            warn!("{} file extension is not a sketch or a fasta/fastq file.", &file);
        }
    }

    let genome_sketches = Mutex::new(vec![]);
    let read_sketches = Mutex::new(vec![]);
    //read c can be lower than lowest genome c. 
    let mut lowest_genome_c = None;
    let mut current_k = None;

    for genome_sketch_file in genome_sketch_files{
        let file = 
                File::open(genome_sketch_file.clone()).expect("Genome sketch {} not a valid file");
        let genome_reader = BufReader::with_capacity(10_000_000, file);
        let genome_sketches_vec: Vec<GenomeSketch> = bincode::deserialize_from(genome_reader).expect(&format!(
            "Genome sketch {} is not a valid sketch.",
            &genome_sketch_file
        ));
        if genome_sketches_vec.is_empty(){
            continue
        }
        let c = genome_sketches_vec.first().unwrap().c;
        let k = genome_sketches_vec.first().unwrap().k;
        if lowest_genome_c.is_none(){
            lowest_genome_c = Some(c);
        }
        else if lowest_genome_c.unwrap() < c{
            lowest_genome_c = Some(c);
        }
        if current_k.is_none(){
            current_k = Some(genome_sketches_vec.first().unwrap().k);
        }
        else if current_k.unwrap() != k{
            error!("Genome sketches have inconsistent -k. Exiting.");
            std::process::exit(1);
        }
        genome_sketches.lock().unwrap().extend(genome_sketches_vec);
    }

    read_sketch_files.into_par_iter().for_each(|read_sketch_file|{
        let file = 
                File::open(read_sketch_file.clone()).expect(&format!("Read sketch {} not a valid file", &read_sketch_file));
        let read_reader = BufReader::with_capacity(10_000_000, file);
        let read_sketch_enc: SequencesSketchEncode = bincode::deserialize_from(read_reader).expect(&format!(
            "Read sketch {} is not a valid sketch.",
            read_sketch_file
        ));
        let read_sketch = SequencesSketch::from_enc(read_sketch_enc);
        if lowest_genome_c.is_some() && read_sketch.c > lowest_genome_c.unwrap(){
            error!("Value of -c for {} is {} -- greater than the smallest value of -c for a genome sketch {}. Exiting.", read_sketch.c, read_sketch_file, lowest_genome_c.unwrap());
            std::process::exit(1);
        }
        read_sketches.lock().unwrap().push(read_sketch);
    });

    genome_files.into_par_iter().for_each(|genome_file|{
        if lowest_genome_c.is_some() && lowest_genome_c.unwrap() < args.c{
            error!("Value of -c for contain is {} -- greater than the smallest value of -c for a genome sketch {}. Continuing without sketching.", args.c, lowest_genome_c.unwrap());
        }
        else if current_k.is_some() && current_k.unwrap() != args.k{
            error!("-k {} is not equal to -k {} found in sketches. Continuing without sketching.", args.k, current_k.unwrap());
        }
        else {
            if args.individual{
            let indiv_gn_sketches = sketch_genome_individual(args.c, args.k, genome_file, args.min_spacing_kmer);
                genome_sketches.lock().unwrap().extend(indiv_gn_sketches);

            }
            else{
                let genome_sketch_opt = sketch_genome(args.c, args.k, &genome_file, args.min_spacing_kmer);
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

    return (read_sketches.into_inner().unwrap(), genome_sketches.into_inner().unwrap()) ;

}

fn get_stats<'a>(
    args: &ContainArgs,
    genome_sketch: &'a GenomeSketch,
    sequence_sketch: &'a SequencesSketch,
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
    for kmer in gn_kmers.iter() {
        if sequence_sketch.kmer_counts.contains_key(kmer) {
            contain_count += 1;
            covs.push(sequence_sketch.kmer_counts[kmer]);
        }
    }
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
    for cov in covs.iter() {
        if (*cov as f64) < max_cov {
            full_covs.push(*cov);
        }
    }
    let mean_cov = full_covs.iter().sum::<u32>() as f64 / full_covs.len() as f64;
    let geq1_mean_cov = covs.iter().sum::<u32>() as f64 / covs.len() as f64;

    let use_lambda;
    if median_cov > MEDIAN_ANI_THRESHOLD {
        use_lambda = AdjustStatus::High
    } else {
        let test_lambda;
        if args.ratio {
            test_lambda = ratio_lambda(&full_covs)
        } else if args.mme {
            test_lambda = mme_lambda(&full_covs)
        } else if args.nb {
            test_lambda = binary_search_lambda(&full_covs)
        } else if args.mle{
            test_lambda = mle_zip(&full_covs, sequence_sketch.k as f64)
        } else{
            test_lambda = ratio_lambda(&full_covs)
        };
        if test_lambda.is_none(){
            use_lambda = AdjustStatus::Low
        }
        else{
            use_lambda = AdjustStatus::Lambda(test_lambda.unwrap());
        }
    }

    let final_est_cov;


    if let AdjustStatus::Lambda(lam) = use_lambda{
        final_est_cov = lam
    }
    else{
        if median_cov < MEDIAN_ANI_THRESHOLD{
            final_est_cov = geq1_mean_cov;
        }
        else{
            final_est_cov = geq1_mean_cov;
        }
    }

    let opt_lambda;
    if use_lambda == AdjustStatus::Low || use_lambda == AdjustStatus::High {opt_lambda = None} else { opt_lambda = Some(final_est_cov) };

    let opt_est_ani =
        ani_from_lambda(opt_lambda, mean_cov, sequence_sketch.k as f64, &full_covs);
    if naive_ani < args.minimum_ani {
        return None;
    }
    let (mut low_ani, mut high_ani, mut low_lambda, mut high_lambda) = (None, None, None, None);
    if args.ci && opt_lambda.is_some() {
        let bootstrap = bootstrap_interval(&full_covs, sequence_sketch.k as f64, &args);
        low_ani = bootstrap.0;
        high_ani = bootstrap.1;
        low_lambda = bootstrap.2;
        high_lambda = bootstrap.3;
    }

    let final_est_ani;
    if opt_lambda.is_none() || opt_est_ani.is_none() || args.no_adj{
        final_est_ani = naive_ani;
    }
    else{
        final_est_ani = opt_est_ani.unwrap();
    }

    let ani_result = AniResult {
        naive_ani,
        final_est_ani,
        final_est_cov,
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
        } else if args.mle {
            lambda = mle_zip(&rand_vec, k);
        }
        else{
            lambda = ratio_lambda(&rand_vec);
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
        let count_p1 = count_map[&(most_ind + 1)] as f64;
        let count = count_map[&most_ind] as f64;
        if count_p1 < MINIMUM_COUNT_RATIO || count < MINIMUM_COUNT_RATIO{
            return None;
        }
        let lambda = Some(
            count_p1 / count
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
