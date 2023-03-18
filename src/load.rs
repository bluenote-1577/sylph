use crate::cmdline::*;
use rand::distributions::WeightedIndex;

use itertools::*;
use puruspe::ln_gamma;
use rand::distributions::Distribution;


use crate::constants;
use crate::io_utils::*;
use crate::types::*;
use log::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::sync::Mutex;
pub const N_STRAIN:usize = 1;

pub fn pois_approx(lambda: f64, k: f64) -> f64 {
    return 1.;
    //N ( μ=λ, σ=√λ)
    let sigma = lambda.powf(0.5);
    return 1. / sigma / (2. * 3.1415_f64).powf(0.5)
        * 2.713_f64.powf(-0.5 * ((k - lambda) / sigma).powf(2.))
        + 0.000001;
}

pub fn load(args: Load) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let genome_file =
        File::open(args.references_sketch.clone()).expect("Genome sketch not a valid file");
    let genome_reader = BufReader::with_capacity(10_000_000, genome_file);
    let genome_sketch: GenomesSketch = bincode::deserialize_from(genome_reader).unwrap();
    info!("Genome sketch loading complete.");

    let mut results = vec![];
    for query_file in args.query_sketch.iter() {
        let res_vec = work(&args, query_file, &genome_sketch);
        info!("Finished computation for {}. Writing...", query_file);
        output_metaphlan_format(&res_vec, genome_sketch.tax_tree.as_ref(), query_file, &args);
        results.push(res_vec);
    }

    print_sample_abund_matrix(&results, &args);
}

fn work<'a>(
    args: &Load,
    query_file: &str,
    genome_sketch: &'a GenomesSketch,
) -> Vec<StrainResult<'a>> {
    let sequence_file = File::open(query_file).expect("Sequence sketch not a valid file");
    let seq_reader = BufReader::with_capacity(10_000_000, sequence_file);
    let sequence_sketch_enc: SequencesSketchEncode = bincode::deserialize_from(seq_reader).unwrap();
    let sequence_sketch = SequencesSketch::from_enc(sequence_sketch_enc);
    info!("Sequence sketch loading complete.");

    if genome_sketch.k != sequence_sketch.k {
        error!(
            "k parameter for query {} != k parameter for reference {}",
            sequence_sketch.k, genome_sketch.k
        );
        std::process::exit(1);
    }
    if genome_sketch.c != sequence_sketch.c {
        error!(
            "c parameter for query {} != c parameter for reference {}",
            sequence_sketch.c, genome_sketch.c
        );
        std::process::exit(1);
    }
    let ani_desc_debug: Mutex<Vec<_>> = Mutex::new(vec![]);
    let anis: Mutex<MMHashMap<_, _>> = Mutex::new(MMHashMap::default());
    let useable_refs: Mutex<Vec<_>> = Mutex::new(vec![]);
    let used_kmers: Mutex<MMHashSet<_>> = Mutex::new(MMHashSet::default());
    let index_vec = (0..genome_sketch.file_names.len()).collect::<Vec<usize>>();
    let progress: Mutex<usize> = Mutex::new(0);
    index_vec.into_par_iter().for_each(|i| {
        let mut count = 0;
        let mut covs = vec![];
        let gn_kmers = &genome_sketch.genome_kmer_counts[i];
        for kmer in gn_kmers.iter() {
            if sequence_sketch.kmer_counts.contains_key(kmer) {
                used_kmers.lock().unwrap().insert(kmer);
                count += 1;
                covs.push(sequence_sketch.kmer_counts[kmer]);
            }
        }
        //        dbg!(count,gn_kmers.len(),sequence_sketch.kmer_counts.len());
        let ani = f64::powf(
            count as f64 / gn_kmers.len() as f64,
            1. / genome_sketch.k as f64,
        );
        covs.sort();
        println!(">{}", genome_sketch.file_names[i]);
        for cov in covs.iter() {
            println!("{}", cov);
        }

        {
            let mut num_proc = progress.lock().unwrap();
            *num_proc += 1;
            if *num_proc % 1000 == 0 && *num_proc > 0 {
                info!("{} references processed", num_proc);
            }
        }
        if ani > args.min_ani {
            let ani = correct_ani_abund(ani, &covs, sequence_sketch.k as f64);
            let median_cov = covs[covs.len() / 2];
            let mean_cov =
                covs.iter().enumerate().map(|(_, x)| *x).sum::<u32>() as f32 / (covs.len() as f32);
            anis.lock().unwrap().insert(
                i,
                (
                    ani,
                    mean_cov,
                    &genome_sketch.file_names[i],
                    genome_sketch.relative_taxid[i],
                    covs.len(),
                ),
            );
            useable_refs.lock().unwrap().push(i);
            ani_desc_debug.lock().unwrap().push((
                ani,
                median_cov,
                mean_cov,
                genome_sketch.file_names[i].clone(),
            ));
        }
    });

    let mut used_kmer_c = 0;
    for kmer in used_kmers.into_inner().unwrap() {
        used_kmer_c += sequence_sketch.kmer_counts[kmer];
    }
    let mut all_kmers_sequence = 0;
    for val in sequence_sketch.kmer_counts.values() {
        all_kmers_sequence += (*val) as usize;
    }
    info!("Number of kmers found in references passing filter: {}. Number of total kmers in query: {}. Fraction: {}", used_kmer_c, all_kmers_sequence, used_kmer_c as f64 / all_kmers_sequence as f64);

    let mut ani_desc_debug = ani_desc_debug.into_inner().unwrap();
    let useable_refs = useable_refs.into_inner().unwrap();
    let anis = anis.into_inner().unwrap();
    ani_desc_debug.sort_by(|x, y| x.partial_cmp(&y).unwrap());
    for ani in ani_desc_debug {
        debug!("{}\t{}\t{}\t{}", ani.0, ani.1, ani.2, ani.3);
    }

    let mut ret_vec = vec![];
    let mut clusters = cluster_genomes(&genome_sketch, &useable_refs, args.cluster);
    for (i, cluster) in clusters.iter_mut().enumerate() {
        cluster.sort();
        debug!("Cluster {:?}", &cluster);
        if cluster.len() > 1 {
            let anis_cluster: Vec<f64> = cluster.iter().map(|x| anis[x].0).collect();
            //        dbg!(&cluster);
            let mut res = em(&sequence_sketch, &genome_sketch, &cluster, &anis_cluster);
            for result in res.iter_mut() {
                result.total_number_kmers = all_kmers_sequence;
                result.em_abund.1 = i;
            }
            ret_vec.extend(res);
        } else {
            let ani_res = &anis[&cluster[0]];
            let res = StrainResult {
                ani: ani_res.0 as f32,
                cov: ani_res.1 as f32,
                length: ani_res.4 as f32,
                em_abund: (1., i),
                kmer_fraction: 1.,
                file_name_ref_gn: ani_res.2,
                tax_id: ani_res.3,
                total_number_kmers: all_kmers_sequence as usize,
            };
            ret_vec.push(res)
        }
    }
    return ret_vec;
}

fn cluster_genomes(
    genome_sketch: &GenomesSketch,
    useable_refs: &Vec<usize>,
    cluster_val: f64,
) -> Vec<Vec<usize>> {
    let mut clusters = vec![];
    let mut used_genomes: MMHashSet<_> = MMHashSet::default();
    let mut kmer_to_genome_table: MMHashMap<u64, Vec<_>> = MMHashMap::default();
    let genome_kmers = &genome_sketch.genome_kmer_counts;

    for i in useable_refs.iter() {
        let kmer_map = &genome_kmers[*i];
        //This may occur for eukaryotic contaminiation
        if kmer_map.len() > 10_000_000 {
            continue;
        }
        for kmer in kmer_map.iter() {
            let genomes = kmer_to_genome_table.entry(*kmer).or_insert(vec![]);
            genomes.push(*i);
        }
    }

    for i in useable_refs.iter() {
        if used_genomes.contains(i) {
            continue;
        }
        used_genomes.insert(*i);

        let mut cluster = vec![];
        cluster.push(*i);

        let mut genome_containment_count: MMHashMap<_, _> = MMHashMap::default();

        let kmers = genome_kmers[*i].iter();
        if kmers.len() > 10_000_00 {
            continue;
        }
        for kmer in kmers {
            let corresp_genomes = &kmer_to_genome_table[kmer];
            for j in corresp_genomes {
                if !used_genomes.contains(j) {
                    let c = genome_containment_count.entry(*j).or_insert(0);
                    *c += 1;
                }
            }
        }

        for (j, count) in genome_containment_count {
            let ratio = count as f64
                / f64::min(genome_kmers[*i].len() as f64, genome_kmers[j].len() as f64);
            let ani = f64::powf(ratio, 1. / genome_sketch.k as f64);
            if ani > cluster_val {
                cluster.push(j);
                used_genomes.insert(j);
            }
        }

        clusters.push(cluster);
    }
    return clusters;
}

fn em<'a>(
    sequence_sketch: &SequencesSketch,
    genome_sketch: &'a GenomesSketch,
    useable_refs: &Vec<usize>,
    anis: &Vec<f64>,
) -> Vec<StrainResult<'a>> {
    let thresh = 200;
    let anis_k: Vec<f64> = anis
        .iter()
        .map(|x| f64::powf(*x, sequence_sketch.k as f64))
        .collect();
    let mut genome_kmer_counts_in_sequence = vec![0.; useable_refs.len()];
    let full_lengths: Vec<f64> = genome_sketch
        .genome_kmer_counts
        .iter()
        .enumerate()
        .filter(|(i, _)| useable_refs.contains(i))
        .map(|(_, x)| x.len() as f64)
        .collect();
    let mut kmer_union: MMHashSet<_> = MMHashSet::default();
    let mut equiv_classes: HashMap<Vec<usize>, (f64, f64)> = HashMap::default();
    let mut kmer_to_refs: HashMap<u64, Vec<usize>> = HashMap::default();
    for (iter_num, i) in useable_refs.iter().enumerate() {
        let kmer_map = &genome_sketch.genome_kmer_counts[*i];
        let mut useable_count = 0;
        for kmer in kmer_map.iter() {
            if sequence_sketch.kmer_counts.contains_key(kmer) {
                useable_count += 1;
                if sequence_sketch.kmer_counts[kmer] > thresh {
                    continue;
                }
                kmer_union.insert(kmer);
                let vec = kmer_to_refs.entry(*kmer).or_insert(vec![]);
                vec.push(iter_num);
            }
        }
        genome_kmer_counts_in_sequence[iter_num] = useable_count as f64;
    }
    let lengths = genome_kmer_counts_in_sequence;

    for (kmer, class) in kmer_to_refs {
        let c = equiv_classes.entry(class).or_insert((0., 0.));
        c.0 += sequence_sketch.kmer_counts[&kmer] as f64;
        c.1 += 1.;
    }

    let mut total_kmers = 0.;
    println!(">all");
    for kmer in kmer_union.iter() {
        total_kmers += sequence_sketch.kmer_counts[kmer] as f64;
        println!("{}", sequence_sketch.kmer_counts[kmer]);
    }
    //    dbg!(useable_refs.len());
    //    dbg!(kmer_union.len(), genome_sketch.genome_kmer_counts[useable_refs[0]].len());
    let test = false;
    let epsilon = 0.000;
    let anis_k_sum: f64 = anis_k.iter().sum();
    let mut thetas: Vec<f64> = anis_k.iter().map(|x| x / anis_k_sum).collect();
    //    let mut thetas = vec![1. / useable_refs.len() as f64; useable_refs.len()];
    for _ in 0..200 {
        let mut new_thetas = vec![0.; useable_refs.len()];
        let avg_cov_kmer = (0..thetas.len())
            .map(|x| thetas[x] * total_kmers / lengths[x])
            .collect::<Vec<_>>();
        for (class, count) in equiv_classes.iter() {
            let mut denom = 0.;
            if test {
                for j in 0..thetas.len() {
                    if class.contains(&j) {
                        denom += (1. - epsilon) / full_lengths[j]
                            * thetas[j]
                            * anis_k[j]
                            * pois_approx(avg_cov_kmer[j], thetas[j] * (count.0 / count.1)) as f64;
                    } else {
                        denom += epsilon / full_lengths[j]
                            * thetas[j]
                            * anis_k[j]
                            * pois_approx(avg_cov_kmer[j], thetas[j] * (count.0 / count.1)) as f64;
                    }
                    //                denom += 1. / lengths[*j] * thetas[*j];
                }
                for i in 0..thetas.len() {
                    let num;
                    if class.contains(&i) {
                        num = (1. - epsilon) * (count.0) as f64 / full_lengths[i]
                            * thetas[i]
                            * anis_k[i]
                            * pois_approx(avg_cov_kmer[i], thetas[i] * (count.0 / count.1)) as f64;
                    } else {
                        num = epsilon * (count.0) as f64 / full_lengths[i]
                            * thetas[i]
                            * anis_k[i]
                            * pois_approx(avg_cov_kmer[i], thetas[i] * (count.0 / count.1)) as f64;
                    }
                    //                let num = (*count)  as f64 / lengths[*i] * thetas[*i];
                    new_thetas[i] += num / denom;
                }
            } else {
                for j in class.iter() {
                    denom += 1. / full_lengths[*j] * thetas[*j] * anis_k[*j];
                    //                denom += 1. / lengths[*j] * thetas[*j];
                }
                for i in class.iter() {
                    let num = (count.0) as f64 / full_lengths[*i] * thetas[*i] * anis_k[*i];
                    //                let num = (*count)  as f64 / lengths[*i] * thetas[*i];
                    new_thetas[*i] += num / denom;
                }
            }
        }
        for theta in new_thetas.iter_mut() {
            *theta /= total_kmers;
        }
        //dbg!(&thetas);
        thetas = new_thetas;
    }
    let mut res = vec![];
    let mut temp = (0..useable_refs.len())
        .into_iter()
        .map(|i| (anis[i], i))
        .collect::<Vec<_>>();

    temp.sort_by(|x, y| y.partial_cmp(&x).unwrap());
    let sorted_use_ref_clone = temp.into_iter().map(|x| x.1).collect::<Vec<_>>();

    let nearest_neighbour_num_kmers =
        genome_sketch.genome_kmer_counts[useable_refs[sorted_use_ref_clone[0]]].len();
    opt_strain(
        &thetas,
        sequence_sketch,
        genome_sketch,
        total_kmers / lengths[0],
        useable_refs,
        &kmer_union,
    );

    //fn opt_strain<'a>(
    //    initial_thetas: Vec<f64>,
    //    sequence_sketch: &SequencesSketch,
    //    genome_sketch: &'a GenomesSketch,
    //    initial_total_cov: f64,
    //    useable_refs: &Vec<usize>,
    //    kmer_union: &MMHashSet<u64>,
    //) {

    //Checkingthe number of k-mers explained by each genome, similar
    //to variance epxlained by PCs. Doesn't work well though...
    let mut used_kmers = MMHashSet::default();
    for i in sorted_use_ref_clone {
        let mut non_used_kmer_count = 0;
        let kmer_map = &genome_sketch.genome_kmer_counts[i];
        for kmer in kmer_map.iter() {
            if sequence_sketch.kmer_counts.contains_key(kmer) {
                if !used_kmers.contains(kmer) {
                    non_used_kmer_count += 1;
                    used_kmers.insert(kmer);
                }
            }
        }
        if thetas[i] > constants::EM_ABUND_CUTOFF {
            let load_res = StrainResult {
                ani: anis[i] as f32,
                cov: (thetas[i] * total_kmers / lengths[i]) as f32,
                length: lengths[i] as f32,
                em_abund: (thetas[i] as f32, usize::MAX),
                file_name_ref_gn: &genome_sketch.file_names[useable_refs[i]],
                tax_id: genome_sketch.relative_taxid[useable_refs[i]],
                kmer_fraction: non_used_kmer_count as f32 / nearest_neighbour_num_kmers as f32,
                total_number_kmers: 0,
            };
            res.push(load_res);
        }
    }
    return res;
}

fn correct_ani_abund(ani: f64, covs: &Vec<u32>, k: f64) -> f64 {
    let thresh = 5;
    let mut leq_threshold = 0;
    for cov in covs.iter() {
        if *cov < thresh {
            leq_threshold += 1;
        }
    }

    if (leq_threshold as f64) < 0.9 * covs.len() as f64 {
        return ani;
    }

    let mut counts: HashMap<_, _> = HashMap::default();
    for cov in covs.iter() {
        let c = counts.entry(cov).or_insert(0);
        *c += 1;
    }

    let mut lambda_ests = vec![];
    let mut weights = vec![];
    for i in 1..10 {
        if counts.contains_key(&i) && counts.contains_key(&(i + 1)) {
            lambda_ests.push(counts[&(i + 1)] as f64 / counts[&i] as f64 * (i + 1) as f64);
            weights.push(f64::min(counts[&(i + 1)] as f64, counts[&i] as f64));
        }
    }
    //    dbg!(&lambda_ests, &counts, &weights);
    let sum = weights.iter().sum::<f64>();
    if covs.len() < 30 || weights.len() == 0 {
        return ani;
    }

    let mut lambda = 0.;
    for i in 0..weights.len() {
        lambda += weights[i] * lambda_ests[i];
    }
    lambda /= sum;
    //    lambda = lambda_ests[0];
    let ret = f64::min(
        1.0,
        ani * f64::powf(1. / (1. - f64::powf(2.718, -lambda)), 1. / k),
    );
    trace!("{},{},{}", lambda, ret, ani);
    return ret;
}

fn opt_strain<'a>(
    initial_thetas: &Vec<f64>,
    sequence_sketch: &SequencesSketch,
    genome_sketch: &'a GenomesSketch,
    initial_total_cov: f64,
    useable_refs: &Vec<usize>,
    kmer_union: &MMHashSet<&u64>,
) {
    let opt_solver = false;
    let mixture_add_model = true;
    let K = useable_refs.len();
    let J = N_STRAIN;
    //Get kmer count matrix
    let mut kmer_count_matrix = vec![vec![0.; useable_refs.len()]; kmer_union.len()];
    let mut kmer_union_vec = kmer_union.iter().map(|x| **x).collect::<Vec<u64>>();
    kmer_union_vec.sort();
    for (i, kmer) in kmer_union_vec.iter().enumerate() {
        for (j, index) in useable_refs.iter().enumerate() {
            if genome_sketch.genome_kmer_counts[*index].contains(kmer) {
                kmer_count_matrix[i][j] = sequence_sketch.kmer_counts[kmer] as f64;
            }
        }
    }

    let mut thetas = (0..J)
        .map(|x| (initial_thetas[x], x))
        .collect::<Vec<(f64, usize)>>();

    dbg!(useable_refs.len());
    if mixture_add_model{
        thetas.sort_by(|x, y| y.partial_cmp(&x).unwrap());
        let mut A = vec![vec![0;J];kmer_union_vec.len()];
        let mut covs = vec![];
        let mut r = 10.;
        for i in 0..J{
            covs.push(thetas[i].0 * initial_total_cov);
        }
        let mut rng = rand::thread_rng();
        let mut kmer_count_vec = vec![];
        for kmer in kmer_union_vec.iter(){
            kmer_count_vec.push(sequence_sketch.kmer_counts[kmer] as f64);
        }
        let mut cov_theta = vec![0.; J * J];
        for i in 0..J {
            cov_theta[i + i * J] = initial_total_cov/20.;
        }
        let mvn_theta = MultivariateNormal::new(vec![0.; J], cov_theta).unwrap();
        let mvn_r = Normal::new(0., 0.5).unwrap();
        let mut best_ll = f64::MIN;
        for _ in 0..10000{
            let theta_samp = mvn_theta.sample(&mut rng);
            let r_samp = mvn_r.sample(&mut rng);
            let mut new_covs = covs.clone();
            for i in 0..J{
                new_covs[i] += theta_samp[i];
                if new_covs[i] < 0.{
                    new_covs[i] = 0.;
                }
            }
            let mut new_r = r.clone();
            new_r += r_samp;
            if new_r < 1.{
                new_r = 1.
            }
            for i in 0..kmer_count_vec.len(){
                let mut tup_probs = vec![];
                let mut tups = vec![];
                for (j,k,l) in iproduct!(0..2,0..2,0..2){
                    if j == 0 && k == 0 && l ==0{
                        continue
                    }
                    let tup = vec![j,k,l];
                    let mut total_cov = 0.;
                    for j in 0..J{
                        total_cov += covs[j] * tup[j] as f64;
                    }
                    if total_cov == 0.{
                        let prob = 0.;
                        tup_probs.push(prob);
                        tups.push(tup);
                        continue;

                    }
                    let mut prob = nb_pmf(r / (total_cov + r), r, kmer_count_vec[i]);
                    if prob.is_nan(){
                       prob = 0.; 
                    }
                    tup_probs.push(prob);
                    tups.push(tup)
                }

                if tup_probs == vec![0.;tup_probs.len()]{
                    dbg!(&tup_probs, &tups, &A[i], &covs, kmer_count_vec[i]);
                    panic!();
                }
                normalize_theta(&mut tup_probs);
                let dist2 = WeightedIndex::new(&tup_probs).unwrap();
                let tup = tups[dist2.sample(&mut rng)].clone();
                A[i] = tup;
            }
            let mut old_ll = 0.;
            let mut new_ll = 0.;
            for i in 0..kmer_count_vec.len(){
                let mut total_cov_old = 0.;
                for j in 0..J{
                    total_cov_old += covs[j] * A[i][j] as f64;
                }
                let prob_old;
                if total_cov_old == 0.{
                    prob_old = 0.0000000000001;
                }
                else{
                    prob_old = nb_pmf(r / (total_cov_old  + r), r, kmer_count_vec[i]);
                }
                old_ll += prob_old.ln();
            }
            for _ in 0..10{
                for i in 0..kmer_count_vec.len(){
                    let mut total_cov_new = 0.;
                    for j in 0..J{
                        total_cov_new += new_covs[j] * A[i][j] as f64;
                    }
                    let prob_new;
                    if total_cov_new == 0.{
                        prob_new = 0.0000000000001;
                    }
                    else{
                        prob_new = nb_pmf(new_r / (total_cov_new + new_r), new_r,  kmer_count_vec[i]);
                    }
                    new_ll += prob_new.ln();
                }
                if new_ll > old_ll{
                    r = new_r;
                    covs = new_covs;
                    break;
                }
                else{
                    let n = Uniform::new(0.0, 1.0).unwrap();
                    let p = (new_ll - old_ll).exp();
                    if n.sample(&mut rng) < p {
                        r = new_r;
                        covs = new_covs;
                        break;
                    }
                }
                if best_ll < new_ll{
                    let aic = (2 * N_STRAIN * kmer_union_vec.len() + N_STRAIN + 1) as f64 - 2. * best_ll;

//                    dbg!(&covs, r, best_ll, &A[0..5], &kmer_count_vec[0..5], aic);
                    best_ll = new_ll;
                }
            }
        }
    }

    else if opt_solver {
        //set up parameters
        let num_strain = J;
        let mut thetas = (0..num_strain)
            .map(|x| initial_thetas[x])
            .collect::<Vec<f64>>();
        let normal = thetas.iter().sum::<f64>();
        thetas = thetas.into_iter().map(|x| x / normal).collect();
        dbg!(&thetas);

        //[thetas]_j,[tau]_{jk},[sigma_mult]_j, c
        let num_params = num_strain * K + num_strain + num_strain + 1;
        let mut tau = vec![vec![0.; num_strain]; num_strain];
        for i in 0..num_strain {
            tau[i][i] = 1.;
        }
        let mut initial_guess = vec![0.; num_params];
        for i in 0..num_strain {
            initial_guess[i] = thetas[i];
            initial_guess[num_strain + thetas.len() * i + i] = 1.;
            initial_guess[num_strain + thetas.len() * num_strain + i] = 0.2;
        }
        initial_guess[num_params - 1] = initial_total_cov;

        let tol = vec![0.001; num_strain + 1];
        let mut opt = nlopt::Nlopt::<_, Vec<Vec<f64>>>::new(
            nlopt::Algorithm::Isres,
            num_params,
            obj_fn,
            nlopt::Target::Maximize,
            kmer_count_matrix,
        );
        let res = opt.set_lower_bounds(&vec![0.; num_params]);
        dbg!(res);
        let mut bounds = vec![1.; num_params];
        bounds[num_params - 1] = 5000.;
        let res = opt.set_upper_bounds(&bounds);
        dbg!(res);
        //    let res =
        //        opt.add_inequality_mconstraint(num_strain + 1, simplex_constraint, vec![vec![0.; J]], &tol);
        let res = opt.add_equality_constraint(theta_constraint, vec![vec![0.; J]], 0.01);
        dbg!(res);

        let res = opt.optimize(&mut initial_guess);
        dbg!(res);
        dbg!(initial_guess);
    }
    //DO MCMC
    else {

        let mut mu_scale_std = vec![0.2; J];
        let mut theta_param = vec![1.; J];
        let mut taus = vec![vec![0.; K]; J];
        let mut cov = initial_total_cov;

        thetas.sort_by(|x, y| y.partial_cmp(&x).unwrap());
//        for i in 0..J {
//            theta_param[thetas[i].1] = thetas[i].0;
//            taus[i][thetas[i].1] = 1.;
//        }
        for i in 0..J {
            if taus[i] == vec![0.; K] {
                taus[i] = vec![1. / (K as f64); K];
            }
        }
        normalize_theta(&mut theta_param);

        let mut old_ll = f64::MIN;
        let mut cov_theta = vec![0.; J * J];
        let mut cov_mu = vec![0.; J * J];
        let mut cov_taus = vec![0.; K * K];

        for i in 0..J {
            cov_theta[i + i * J] = 0.01;
            cov_mu[i + i * J] = 0.01;
        }
        for j in 0..K {
            cov_taus[j + j * K] = 0.01;
        }
        let mvn_theta = MultivariateNormal::new(vec![0.; J], cov_theta).unwrap();
        let mvn_mu = MultivariateNormal::new(vec![0.; J], cov_mu).unwrap();
        let mvn_taus = MultivariateNormal::new(vec![0.; K], cov_taus).unwrap();
        let mvn_cov = Normal::new(0., 1.).unwrap();

        let mut max_ll = f64::MIN;
        for x in 0..100000 {
//            dbg!(&theta_param, &mu_scale_std, &taus, &cov, &old_ll);
            if old_ll > max_ll{
                dbg!("NEW");
                max_ll = old_ll;
                dbg!(&theta_param, &old_ll, &mu_scale_std, &cov, &taus);
            }

            let mut proposed_theta = theta_param.clone();
            let mut proposed_mu_scale_std = mu_scale_std.clone();
            let mut proposed_taus = taus.clone();
            let mut proposed_cov = cov.clone();
            let mut r = rand::thread_rng();
            let theta_samp = mvn_theta.sample(&mut r);
            let mu_samp = mvn_mu.sample(&mut r);
            let mut tau_samps = vec![];
            for _ in 0..J {
                let tau_samp = mvn_taus.sample(&mut r);
                tau_samps.push(tau_samp);
            }
            let cov_samp = mvn_cov.sample(&mut r);
            for j in 0..theta_samp.len() {
                proposed_theta[j] += theta_samp[j];
                if proposed_theta[j] < 0. {
                    proposed_theta[j] = 0.;
                }
                //Worried about mu_scale becoming too small, invalidates distribution.
                proposed_mu_scale_std[j] += mu_samp[j];
                if proposed_mu_scale_std[j] < 0.0 {
                    proposed_mu_scale_std[j] = 0.0;
                }
                for k in 0..K {
                    proposed_taus[j][k] += tau_samps[j][k];
                    if proposed_taus[j][k] < 0. {
                        proposed_taus[j][k] = 0.;
                    }
                }
                normalize_theta(&mut proposed_taus[j]);
            }
            normalize_theta(&mut proposed_theta);
            proposed_cov += cov_samp;
            if proposed_cov < 1. {
                proposed_cov = 1.;
            }
            (0..proposed_theta.len())
                .map(|i| proposed_theta[i] + theta_samp[i])
                .collect::<Vec<f64>>();
            let mut concat_vec = vec![];
            concat_vec.append(&mut proposed_theta.clone());
            for j in 0..J {
                concat_vec.append(&mut proposed_taus[j].clone());
            }
            concat_vec.append(&mut proposed_mu_scale_std.clone());
            concat_vec.push(cov);

            let ll = obj_fn(&concat_vec, None, &mut kmer_count_matrix);
            if ll > old_ll {
                old_ll = ll;
                theta_param = proposed_theta;
                mu_scale_std = proposed_mu_scale_std;
                taus = proposed_taus;
                cov = proposed_cov;
            } else {
                let n = Uniform::new(0.0, 1.0).unwrap();
                let p = (ll - old_ll).exp();
                if n.sample(&mut r) < p {
                    old_ll = ll;
                    theta_param = proposed_theta;
                    mu_scale_std = proposed_mu_scale_std;
                    taus = proposed_taus;
                    cov = proposed_cov;
                }
            }
//            dbg!(&theta_param);
        }
    }
}

fn normalize_theta(thetas: &mut Vec<f64>) {
    let norm = thetas.iter().sum::<f64>();
    
    let len = thetas.len() ;
    for theta in thetas.iter_mut() {
        if norm == 0.{
            *theta = 1./len as f64;
        }
        else{
            *theta = *theta / norm;
        }
    }
}

fn pr_from_musigma(mu: f64, sigma: f64) -> (f64, f64) {
    let p = (sigma * sigma - mu) / (sigma * sigma);
    let r = (mu * mu) / (sigma * sigma - mu);
    return (p, r);
}

fn nb_pmf(p: f64, r: f64, k: f64) -> f64 {
    let top = ln_gamma(r + k);
    if top.is_infinite(){
        dbg!(p,r,k);
        panic!();
    }
    let bot = ln_gamma(k + 1.) + ln_gamma(r);
    if bot.is_infinite(){
        dbg!(p,r,k);
        panic!();
    }
    let log_int = top - bot + r * p.ln() + k * (1.-p).ln();
    if log_int.is_infinite() || log_int.is_nan(){
        dbg!(top,bot, p, r);
        panic!();
    }
//    dbg!(top,bot,log_int, log_int.exp());
    return log_int.exp();

//    return top / bot * p.powf(k) * (1. - p).powf(r);
}

fn obj_fn(x: &[f64], y: Option<&mut [f64]>, z: &mut Vec<Vec<f64>>) -> f64 {
    //    let mus = (0..5)
    //        .map(|a| x[5 * z[0].len() + 5 + a])
    //        .collect::<Vec<f64>>();
    //    let std = (0..5)
    //        .map(|a| x[5 * z[0].len() + 2 * 5 + a])
    //        .collect::<Vec<f64>>();
    //    let prs = (0..5)
    //        .map(|a| pr_from_musigma(mus[a], std[a]))
    //        .collect::<Vec<(f64, f64)>>();
    //    let nbs = (0..5)
    //        .map(|a| NegativeBinomial::new(prs[a].0, prs[a].1).unwrap())
    //        .collect::<Vec<_>>();
//    dbg!(&x);
    let mut ret = 0.;
    for i in 0..z.len() {
        let mut sum: f64 = 0.;
        for j in 0..N_STRAIN {
            if x[j] == 0.{
                continue
            }
            for k in 0..z[0].len() {
                let mut mu = 0.;
                if z[i][k] != 0.{
                    mu = 1. * x[j] * x.iter().last().unwrap();
                }
                else{
                    continue;
                }
                let sigma_mult = x[N_STRAIN + N_STRAIN * z[0].len() + j];
                let sigma = sigma_mult * mu;
                if sigma*sigma < mu{
                    continue
                }
                let (p, r) = pr_from_musigma(mu, sigma);
                let nbpmf = nb_pmf(p, r, z[i][k] );

                if nbpmf.is_nan(){
//                    dbg!(nb_pmf(p, r, z[i][k]), x[j * 5 + k], x[j]);
//                    dbg!(i,j,k);
                    //panic!("{},{},{},{:?},{},{}",sigma_mult, sigma, mu, z[j], p, r);
                    continue;
                }
                else{
                    sum += nbpmf * x[j * N_STRAIN + k] * x[j];
                    if sum > 1.{
                        dbg!(x[j], x[j * N_STRAIN + k], nbpmf);
                    }
                }
            }
        }
        ret += sum.ln();
    }
    return ret;
}

fn mu_sigma_from_pr(p: f64, r:f64) -> (f64,f64){
    let mu = r * (1. - p) / p;
    let sigma = (r * (1. - p)/ p / p).sqrt();
    return (mu,sigma)

}

fn simplex_constraint(
    result: &mut [f64],
    x: &[f64],
    grad: Option<&mut [f64]>,
    z: &mut Vec<Vec<f64>>,
) {
    let mut sum = 1.;
    for i in 0..N_STRAIN {
        sum += x[i];
    }
    result[0] = sum;

    for i in 1..6 {
        let mut sum = -1.;
        for j in (i * z[0].len())..(i * z[0].len() + z[0].len()) {
            sum -= x[j];
        }
        result[i] = sum;
    }
    dbg!(&x, &result);
}

fn theta_constraint(x: &[f64], y: Option<&mut [f64]>, z: &mut Vec<Vec<f64>>) -> f64 {
    return -1. + x[0] + x[1] + x[2] + x[3] + x[4];
}
