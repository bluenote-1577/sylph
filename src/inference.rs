use statrs::function::gamma::*;
use fxhash::FxHashMap;
use crate::constants::*;
use std::collections::HashSet;

pub fn r_from_moments_lambda(m: f64, v: f64, lambda: f64) -> f64{
    //return (v / m - 1. - lambda + m) / lambda;
    //return 1000.;
    return lambda / (v - 1. + lambda + m)
}

pub fn ratio_formula(val: f64, r: f64, lambda: f64) -> f64{
    if r < 100.{
       return gamma(r + val + 1.) / (val + 1.) / gamma(r + val) * lambda / (r + lambda)
    }
    else{
       return (r + val + 1.) / (val + 1.) * lambda / (r + lambda)
    }
}

fn ratio_from_moments_lambda(val: f64, lambda: f64, m: f64, v: f64) -> Option<f64>{
    let r = r_from_moments_lambda(m, v, lambda);
    if r < 0.{
        return None;
    }
    return Some(ratio_formula(val, r, lambda));
}

pub fn binary_search_lambda(full_covs: &[u32]) -> Option<f64>{
    if full_covs.len() == 0{
        return None
    }
    let m = mean(full_covs).unwrap();
    let v = var(full_covs).unwrap();
    let mut _nonzero = 0;
    let mut ones = 0;
    let mut twos = 0;

    for x in full_covs{
        if *x != 0{
            _nonzero += 1;
        }
        if *x == 1{
            ones += 1;
        }
        else if *x == 2{
            twos += 1;
        }
    }



    let ratio_est = twos as f64 / ones as f64;

    let left = f64::max(0.003, m - 2.);
    let right = m + 5.;
    let endpoints = (left,right);
    let mut best = None;
    let mut best_val = 10000.;
    for i in 0..10000{
        let test = (endpoints.1 - endpoints.0)/10000. * i as f64 + endpoints.0;
        let proposed = ratio_from_moments_lambda(1.,test , m, v);
        if proposed.is_some(){
            let p = proposed.unwrap() - ratio_est;
            if p.abs() < best_val{
                best_val = p.abs();
                best = Some(test);
            }
        }
    }
    if best.is_none(){
        return None
    }
    let best = best.unwrap();
    let r = r_from_moments_lambda(m,v,best);
    dbg!(m,v, ratio_est, r, best);
//    let ratio_adj = 1. - pi;
//    dbg!(best,best_val, r, nonzero, full_covs.len(), f64::powf(nonzero as f64 / full_covs.len() as f64, 1./k));
//    dbg!(1. - m / best, f64::powf(1. - m/best, 1. / k as f64));
//    dbg!(zeros_from_nb, ratio_adj, f64::powf(ratio_adj, 1. / k));
//    let val = 1.;
//    let mut endpoints_output = (ratio_from_moments_lambda(1., endpoints.0, m, v) - ratio_est, ratio_from_moments_lambda(1., endpoints.1, m, v) - ratio_est);
//    if endpoints_output.0 < 0. || endpoints_output.1 > 0.{
//        return None;
//    }
//    dbg!(endpoints_output);
//    for _ in 0..100{
//        let proposed = ratio_from_moments_lambda(1., (endpoints.1 + endpoints.0)/2., m, v) - ratio_est;
//        if proposed > 0.{
//            endpoints.0 = (endpoints.1 + endpoints.0)/2.;
//            endpoints_output.0 = proposed;
//        }
//        else{
//            endpoints.1 = (endpoints.1 + endpoints.0)/2.;
//            endpoints_output.1 = proposed;
//        }
//        curr = endpoints_output.0;
//        dbg!(endpoints, endpoints_output, proposed);
//    }
//
    return Some(best);
}

pub fn var(data: &[u32]) -> Option<f64> {
    if data.is_empty(){
        return None
    }
    let mean = mean(data).unwrap();
    let mut var = 0.;
    for x in data{
        var += (*x as f64 - mean) * (*x as f64 - mean)
    }
    return Some(var / data.len() as f64);
}

pub fn mean(data: &[u32]) -> Option<f64> {
    let sum = data.iter().sum::<u32>() as f64;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f64),
        _ => None,
    }
}

pub fn mme_lambda(full_covs: &[u32]) -> Option<f64> {
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

pub fn mle_zip(full_covs: &[u32], _k: f64) -> Option<f64> {
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
        let t2 = mean * (1. - f64::exp(-curr));
        let t3 = 1. - rat;
        let t4 = mean * (f64::exp(-curr));
        curr = curr - (t1 - t2) / (t3 - t4);
    }
    return curr;
}

pub fn ratio_lambda(full_covs: &Vec<u32>, min_count_correct: f64) -> Option<f64> {
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
