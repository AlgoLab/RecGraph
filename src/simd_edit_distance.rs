use std::time::*;
use std::{arch::x86_64::*, intrinsics::transmute};
pub fn exec_no_simd(read: &Vec<u8>, reference: &Vec<u8>) -> i32 {
    let mut m: Vec<Vec<i32>> = vec![vec![0; read.len()]; reference.len()];
    for i in 0..reference.len() {
        m[i][0] = i as i32
    }
    for j in 0..read.len() {
        m[0][j] = j as i32
    }
    for i in 1..reference.len() {
        for j in 1..read.len() {
            let l = m[i][j - 1] + 1;
            let u = m[i - 1][j] + 1;
            let d = m[i - 1][j - 1] + if read[j] == reference[i] { 0 } else { 1 };

            m[i][j] = *[l, u, d].iter().min().unwrap();
        }
    }
    //m.iter().for_each(|line| println!("{:?}", line));
    //println!();
    m[reference.len() - 1][read.len() - 1]
}
//TODO: l value not implemented yet
#[target_feature(enable = "sse2")]
pub unsafe fn exec_simd(read: &Vec<u8>, reference: &Vec<u8>) -> i32 {
    let mut m: Vec<Vec<i32>> = vec![vec![0; read.len()]; reference.len()];
    for i in 0..reference.len() {
        m[i][0] = i as i32
    }
    for j in 0..read.len() {
        m[0][j] = j as i32
    }

    let one_simd: __m128i = _mm_set1_epi32(1);
    for i in 1..reference.len() {
        for j in (1..read.len()).step_by(4) {
            let us = _mm_add_epi32(_mm_loadu_epi32(m[i - 1].get_unchecked(j)), one_simd);
            let i_simd: __m128i = _mm_set1_epi32(reference[i] as i32);
            let read_simd = set_simd_read_slice(read, j);
            //let read_simd = _mm_set1_epi32(read[j] as i32);
            let eq_diff_mask = _mm_cmpeq_epi32(i_simd, read_simd);
            let eq_ds = _mm_load_epi32(m[i - 1].get_unchecked(j - 1));
            let neq_ds = _mm_add_epi32(_mm_loadu_epi32(m[i - 1].get_unchecked(j - 1)), one_simd);
            let ds = _mm_or_epi32(
                _mm_and_si128(eq_diff_mask, eq_ds),
                _mm_andnot_si128(eq_diff_mask, neq_ds),
            );
            let best_choice = _mm_cmplt_epi32(ds, us);
            let result = _mm_or_epi32(
                _mm_and_si128(best_choice, ds),
                _mm_andnot_si128(best_choice, us),
            );

            _mm_storeu_epi32(m[i].get_unchecked_mut(j), result);
        }
    }
    //m.iter().for_each(|line| println!("{:?}", line));
    //println!();

    m[reference.len() - 1][read.len() - 1]
}
#[cfg(target_arch = "x86_64")]
unsafe fn set_simd_read_slice(read: &Vec<u8>, j: usize) -> __m128i {
    let mut idx = j;
    let mut read_as_i32 = vec![];
    while idx < read.len() {
        read_as_i32.push(read[idx] as i32);
        idx += 1;
    }
    while read_as_i32.len() < 4 {
        read_as_i32.push(0);
    }
    let simd_read_slice =
        transmute::<[i32; 4], __m128i>(read_as_i32[0..4].try_into().expect("err"));
    simd_read_slice
}

#[cfg(target_arch = "x86_64")]
pub unsafe fn prova() {
    let s1 = "TGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATATGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATAA".chars().map(|c| c as u8).collect::<Vec<u8>>();
    let s2 = "TGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATATGATATAAAGAAATGAGATTTATTGCCTTGTGGGGGGAAGGGATGTGGTTGTGATAA".chars().map(|c| c as u8).collect::<Vec<u8>>();

    let now1 = Instant::now();
    exec_no_simd(&s1, &s2);
    let t1 = now1.elapsed();

    let now2 = Instant::now();
    exec_simd(&s1, &s2);
    let t2 = now2.elapsed();

    println!("No Simd: {}, Simd: {}", t1.as_micros(), t2.as_micros());
}
#[cfg(test)]
mod tests {

    #[test]
    fn test_no_simd_version() {
        let s1 = "AAA".chars().map(|c| c as u8).collect::<Vec<u8>>();
        let s2 = "AAAA".chars().map(|c| c as u8).collect::<Vec<u8>>();

        let s3 = "ACAAA".chars().map(|c| c as u8).collect::<Vec<u8>>();
        let s4 = "AAAAC".chars().map(|c| c as u8).collect::<Vec<u8>>();

        let ed1 = super::exec_no_simd(&s1, &s2);
        let ed2 = super::exec_no_simd(&s3, &s4);

        assert_eq!(ed1, 1);
        assert_eq!(ed2, 2);
    }
}
