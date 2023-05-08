use crate::types::*;
//use smallvec::SmallVec;
use std::arch::x86_64::*;

pub fn decode(byte: u64) -> u8 {
    if byte == 0 {
        return b'A';
    } else if byte == 1 {
        return b'C';
    } else if byte == 2 {
        return b'G';
    } else if byte == 3 {
        return b'T';
    } else {
        panic!("decoding failed")
    }
}
pub fn print_string(kmer: u64, k: usize) {
    let mut bytes = vec![];
    let mask = 3;
    for i in 0..k {
        let val = kmer >> 2 * i;
        let val = val & mask;
        bytes.push(decode(val));
    }
    dbg!(std::str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
}
#[inline]
fn _position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

//pub fn fmh_seeds(
//    string: &[u8],
//    sketch_params: &SketchParams,
//    contig_index: ContigIndex,
//    new_sketch: &mut Sketch,
//    seed: bool,
//) {
//    if seed && new_sketch.kmer_seeds_k.is_none() {
//        new_sketch.kmer_seeds_k = Some(KmerSeeds::default());
//    }
//    let marker_k = K_MARKER_DNA;
//    let kmer_seeds_k = &mut new_sketch.kmer_seeds_k;
//    let marker_seeds = &mut new_sketch.marker_seeds;
//    let k = sketch_params.k;
//    let c = sketch_params.c;
//    if k > 16 {
//        panic!("Value of k > {} for DNA; not allowed.", marker_k);
//    }
//    if string.len() < 2 * marker_k {
//        return;
//    }
//    let mut rolling_kmer_f_marker: MarkerBits = 0;
//    let mut rolling_kmer_f_seed: MarkerBits;
//    let mut rolling_kmer_r_marker: MarkerBits = 0;
//    let mut rolling_kmer_r_seed: MarkerBits;
//    let seed_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * k);
//
//    let marker_reverse_shift_dist = 2 * (marker_k - 1);
//    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
//    let marker_rev_mask = !(0 | (3 << 2 * marker_k - 2));
//    let len = string.len();
//    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
//    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;
//
//    let threshold = u64::MAX / (c as u64);
//    let threshold_marker = u64::MAX / (sketch_params.marker_c as u64);
//    for i in 0..marker_k - 1 {
//        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
//        //        let nuc_f = KmerEnc::encode(string[i]
//        let nuc_r = 3 - nuc_f;
//        rolling_kmer_f_marker <<= 2;
//        rolling_kmer_f_marker |= nuc_f;
//        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
//        rolling_kmer_r_marker >>= 2;
//        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
//    }
//    for i in marker_k..len {
//        let nuc_byte = string[i] as usize;
//        let nuc_f = BYTE_TO_SEQ[nuc_byte];
//        let nuc_r = 3 - nuc_f;
//        rolling_kmer_f_marker <<= 2;
//        rolling_kmer_f_marker |= nuc_f;
//        rolling_kmer_f_marker &= marker_mask;
//        rolling_kmer_r_marker >>= 2;
//        rolling_kmer_r_marker &= marker_rev_mask;
//        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
//        //        rolling_kmer_r &= max_mask;
//        //        KmerEnc::print_string(rolling_kmer_f, k);
//        //        KmerEnc::print_string(rolling_kmer_r, k);
//        //
//        rolling_kmer_f_seed = rolling_kmer_f_marker & seed_mask;
//        rolling_kmer_r_seed = rolling_kmer_r_marker & seed_mask;
//        let canonical_seed = rolling_kmer_f_seed < rolling_kmer_r_seed;
//
//        let canonical_kmer_seed = if canonical_seed {
//            rolling_kmer_f_seed
//        } else {
//            rolling_kmer_r_seed
//        };
//
////        let hash_seed = mm_hashi64(canonical_kmer_seed as i64);
//        let hash_seed = mm_hash64(canonical_kmer_seed);
//        if hash_seed < threshold {
//            if seed {
//                let kmer_seeds = &mut kmer_seeds_k.as_mut().unwrap();
//                let kmer_positions = kmer_seeds
//                    //Since we fix k = 15,can represent seeds as 32bits
//                    .entry(canonical_kmer_seed as SeedBits)
//                    .or_insert(SmallVec::<[SeedPosition; SMALL_VEC_SIZE]>::new());
//                //                    .or_insert(vec![]);
//                kmer_positions.push(SeedPosition {
//                    pos: i as GnPosition,
//                    canonical: canonical_seed,
//                    contig_index,
//                    phase: 0,
//                });
//            }
//            let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
//            let canonical_kmer_marker = if canonical_marker {
//                rolling_kmer_f_marker
//            } else {
//                rolling_kmer_r_marker
//            };
//
//            if hash_seed < threshold_marker {
//                marker_seeds.insert(canonical_kmer_marker);
//            }
//        }
//    }
//}

#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    let mut key = kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    return key;
}
//use bio::data_structures::interval_tree::IntervalTree;
//use fxhash::{hash, FxHashMap, FxHashSet};
//
#[target_feature(enable = "avx2")]
pub unsafe fn extract_markers_avx2(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    if string.len() < k {
        return;
    }
    let len = (string.len() - k + 1) / 4;
    let string1 = &string[0..len + k - 1];
    let string2 = &string[len..2 * len + k - 1];
    let string3 = &string[2 * len..3 * len + k - 1];
    let string4 = &string[3 * len..4 * len + k - 1];
    if string.len() < 2 * k {
        return;
    }

    let use_40 = if 2 * (k - 1) == 40 {
        true
    } else if 2 * (k - 1) == 60 {
        false
    } else {
        panic!()
    };
    const TWO_K_MINUS_2_60: i32 = 60;
    const TWO_K_MINUS_2_40: i32 = 40;
    let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_kmer_r_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);
    for i in 0..k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);

        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);

        let shift_nuc_r;
        if use_40 {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
        } else {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
        }
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);
    }

    let marker_mask = (Kmer::MAX >> (std::mem::size_of::<Kmer>() * 8 - 2 * k)) as i64;
    let rev_marker_mask: i64 = !(0 | (3 << 2 * k - 2));
    //    let rev_marker_mask = i64::from_le_bytes(rev_marker_mask.to_le_bytes());
    //    dbg!(u64::MAX / (c as u64));
    //    dbg!((u64::MAX / (c as u64)) as i64);
    let threshold_marker = u64::MAX / c as u64;

    let mm256_marker_mask = _mm256_set_epi64x(marker_mask, marker_mask, marker_mask, marker_mask);
    let mm256_rev_marker_mask = _mm256_set_epi64x(
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
    );

    for i in k - 1..(len + k - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);
        rolling_kmer_f_marker = _mm256_and_si256(rolling_kmer_f_marker, mm256_marker_mask);

        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);
        let shift_nuc_r;
        if use_40 {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
        } else {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
        }
        rolling_kmer_r_marker = _mm256_and_si256(rolling_kmer_r_marker, mm256_rev_marker_mask);
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);

        let compare_marker = _mm256_cmpgt_epi64(rolling_kmer_r_marker, rolling_kmer_f_marker);

        let canonical_markers_256 =
            _mm256_blendv_epi8(rolling_kmer_r_marker, rolling_kmer_f_marker, compare_marker);

        //        dbg!(rolling_kmer_f_marker,rolling_kmer_r_marker);
        //        dbg!(print_string(u64::from_ne_bytes(_mm256_extract_epi64(rolling_kmer_f_marker,1).to_ne_bytes()), 31));
        let hash_256 = mm_hash256(canonical_markers_256);
        let v1 = _mm256_extract_epi64(hash_256, 0) as u64;
        let v2 = _mm256_extract_epi64(hash_256, 1) as u64;
        let v3 = _mm256_extract_epi64(hash_256, 2) as u64;
        let v4 = _mm256_extract_epi64(hash_256, 3) as u64;
        //        let threshold_256 = _mm256_cmpgt_epi64(cmp_thresh, hash_256);
        //        let m1 = _mm256_extract_epi64(threshold_256, 0);
        //        let m2 = _mm256_extract_epi64(threshold_256, 1);
        //        let m3 = _mm256_extract_epi64(threshold_256, 2);
        //        let m4 = _mm256_extract_epi64(threshold_256, 3);

        if v1 < threshold_marker {
            kmer_vec.push(v1 as u64);
        }
        if v2 < threshold_marker {
            kmer_vec.push(v2 as u64);
        }
        if v3 < threshold_marker {
            kmer_vec.push(v3 as u64);
        }
        if v4 < threshold_marker {
            kmer_vec.push(v4 as u64);
        }
    }
}

#[target_feature(enable = "avx2")]
pub unsafe fn extract_markers_avx2_positions(string: &[u8], kmer_vec: &mut Vec<(usize,u64)>, c: usize, k: usize) {
    if string.len() < k {
        return;
    }
    let len = (string.len() - k + 1) / 4;
    let string1 = &string[0..len + k - 1];
    let string2 = &string[len..2 * len + k - 1];
    let string3 = &string[2 * len..3 * len + k - 1];
    let string4 = &string[3 * len..4 * len + k - 1];
    if string.len() < 2 * k {
        return;
    }

    let use_40 = if 2 * (k - 1) == 40 {
        true
    } else if 2 * (k - 1) == 60 {
        false
    } else {
        panic!()
    };
    const TWO_K_MINUS_2_60: i32 = 60;
    const TWO_K_MINUS_2_40: i32 = 40;
    let mut rolling_kmer_f_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let mut rolling_kmer_r_marker = _mm256_set_epi64x(0, 0, 0, 0);
    let rev_sub = _mm256_set_epi64x(3, 3, 3, 3);
    for i in 0..k - 1 {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);

        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);

        let shift_nuc_r;
        if use_40 {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
        } else {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
        }
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);
    }

    let marker_mask = (Kmer::MAX >> (std::mem::size_of::<Kmer>() * 8 - 2 * k)) as i64;
    let rev_marker_mask: i64 = !(0 | (3 << 2 * k - 2));
    //    let rev_marker_mask = i64::from_le_bytes(rev_marker_mask.to_le_bytes());
    //    dbg!(u64::MAX / (c as u64));
    //    dbg!((u64::MAX / (c as u64)) as i64);
    let threshold_marker = u64::MAX / c as u64;

    let mm256_marker_mask = _mm256_set_epi64x(marker_mask, marker_mask, marker_mask, marker_mask);
    let mm256_rev_marker_mask = _mm256_set_epi64x(
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
        rev_marker_mask,
    );

    for i in k - 1..(len + k - 1) {
        let nuc_f1 = BYTE_TO_SEQ[string1[i] as usize] as i64;
        let nuc_f2 = BYTE_TO_SEQ[string2[i] as usize] as i64;
        let nuc_f3 = BYTE_TO_SEQ[string3[i] as usize] as i64;
        let nuc_f4 = BYTE_TO_SEQ[string4[i] as usize] as i64;
        let f_nucs = _mm256_set_epi64x(nuc_f4, nuc_f3, nuc_f2, nuc_f1);
        let r_nucs = _mm256_sub_epi64(rev_sub, f_nucs);

        rolling_kmer_f_marker = _mm256_slli_epi64(rolling_kmer_f_marker, 2);
        rolling_kmer_f_marker = _mm256_or_si256(rolling_kmer_f_marker, f_nucs);
        rolling_kmer_f_marker = _mm256_and_si256(rolling_kmer_f_marker, mm256_marker_mask);

        rolling_kmer_r_marker = _mm256_srli_epi64(rolling_kmer_r_marker, 2);
        let shift_nuc_r;
        if use_40 {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_40);
        } else {
            shift_nuc_r = _mm256_slli_epi64(r_nucs, TWO_K_MINUS_2_60);
        }
        rolling_kmer_r_marker = _mm256_and_si256(rolling_kmer_r_marker, mm256_rev_marker_mask);
        rolling_kmer_r_marker = _mm256_or_si256(rolling_kmer_r_marker, shift_nuc_r);

        let compare_marker = _mm256_cmpgt_epi64(rolling_kmer_r_marker, rolling_kmer_f_marker);

        let canonical_markers_256 =
            _mm256_blendv_epi8(rolling_kmer_r_marker, rolling_kmer_f_marker, compare_marker);

        //        dbg!(rolling_kmer_f_marker,rolling_kmer_r_marker);
        //        dbg!(print_string(u64::from_ne_bytes(_mm256_extract_epi64(rolling_kmer_f_marker,1).to_ne_bytes()), 31));
        let hash_256 = mm_hash256(canonical_markers_256);
        let v1 = _mm256_extract_epi64(hash_256, 0) as u64;
        let v2 = _mm256_extract_epi64(hash_256, 1) as u64;
        let v3 = _mm256_extract_epi64(hash_256, 2) as u64;
        let v4 = _mm256_extract_epi64(hash_256, 3) as u64;
        //        let threshold_256 = _mm256_cmpgt_epi64(cmp_thresh, hash_256);
        //        let m1 = _mm256_extract_epi64(threshold_256, 0);
        //        let m2 = _mm256_extract_epi64(threshold_256, 1);
        //        let m3 = _mm256_extract_epi64(threshold_256, 2);
        //        let m4 = _mm256_extract_epi64(threshold_256, 3);

        if v1 < threshold_marker {
            kmer_vec.push((i, v1 as u64));
        }
        if v2 < threshold_marker {
            kmer_vec.push((len + i, v2 as u64));
        }
        if v3 < threshold_marker {
            kmer_vec.push((2*len + i, v3 as u64));
        }
        if v4 < threshold_marker {
            kmer_vec.push((3*len + i, v4 as u64));
        }
    }
}
