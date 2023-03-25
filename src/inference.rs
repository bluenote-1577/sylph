use crate::cmdline::*;
use statrs::function::gamma::*;

fn r_from_moments_lambda(m: f64, v: f64, lambda: f64) -> f64{
    //return (v / m - 1. - lambda + m) / lambda;
    //return 1000.;
    return lambda / (v - 1. + lambda + m)
}

fn ratio_formula(val: f64, r: f64, lambda: f64) -> f64{
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
    let mut nonzero = 0;
    let mut ones = 0;
    let mut twos = 0;

    for x in full_covs{
        if *x != 0{
            nonzero += 1;
        }
        if *x == 1{
            ones += 1;
        }
        else if *x == 2{
            twos += 1;
        }
    }

    let zeros = full_covs.len() as f64 - nonzero as f64;


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
    let zeros_from_nb = f64::powf(r/(r + best), r);
    let pi = zeros / full_covs.len() as f64 - zeros_from_nb;
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
