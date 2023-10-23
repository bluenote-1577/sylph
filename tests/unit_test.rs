use assert_cmd::prelude::*; // Add methods on commands
use sylph::seeding;

#[test]
fn test_hash(){

    let key = 19238239812933123;
    println!("{}", format!("{key:b}"));
    let h = seeding::mm_hash64(key);
    println!("{}", format!("{h:b}"));
    let rev = seeding::rev_hash_64(h);
    println!("{}", format!("{rev:b}"));
    assert!(rev == key);

    if is_x86_feature_detected!("avx2"){
        unsafe{
            let key = key as i64;
            println!("{}", format!("{key:b}"));
            use std::arch::x86_64::*;
            use sylph::avx2_seeding::*;
            let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, key);
                let hash_256 = mm_hash256(rolling_kmer_f_marker);
                let v1 = _mm256_extract_epi64(hash_256, 0);
                println!("{}", format!("{v1:b}"));
                assert!(v1 == h as i64);
        }

    }
}
