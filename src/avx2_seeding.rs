use std::arch::x86_64::*;
use crate::types::*;

#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    let mut key = kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    //TODO this is bugged. Fix after release...
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
pub unsafe fn extract_markers_avx2_positions(string: &[u8], kmer_vec: &mut Vec<(usize, usize,u64)>, c: usize, k: usize, contig_number: usize) {
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
            kmer_vec.push((contig_number, i, v1 as u64));
        }
        if v2 < threshold_marker {
            kmer_vec.push((contig_number, len + i, v2 as u64));
        }
        if v3 < threshold_marker {
            kmer_vec.push((contig_number, 2*len + i, v3 as u64));
        }
        if v4 < threshold_marker {
            kmer_vec.push((contig_number, 3*len + i, v4 as u64));
        }
    }
}
