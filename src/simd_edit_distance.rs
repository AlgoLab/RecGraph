use std::arch::x86_64::*;
use std::cmp;
//use std::ptr;

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

#[target_feature(enable = "avx2")]
pub unsafe fn exec_simd_avx2(read: &Vec<u8>, reference: &Vec<u8>) -> f32 {
    let mut m: Vec<Vec<f32>> = vec![vec![0f32; read.len()]; reference.len()];
    for i in 0..reference.len() {
        m[i][0] = i as f32
    }
    for j in 0..read.len() {
        m[0][j] = j as f32
    }
    let one_simd = _mm256_set1_ps(1.0);
    let max_multiple = (read.len() / 8) * 8;
    let read_f32 = read[0..max_multiple + 1]
        .iter()
        .map(|c| *c as f32)
        .collect::<Vec<f32>>();
    for i in 1..reference.len() {
        for j in (1..max_multiple + 1).step_by(8) {
            let us = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j)), one_simd);

            let eq_char = _mm256_cmp_ps(
                _mm256_loadu_ps(read_f32.get_unchecked(j)), //read chars simd
                _mm256_set1_ps(reference[i] as f32),        // reference char simd
                _CMP_EQ_OS,
            );
            let neq_ds = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)), one_simd);
            let eq_ds = _mm256_loadu_ps(m[i - 1].get_unchecked(j - 1));
            let ds = _mm256_blendv_ps(neq_ds, eq_ds, eq_char);

            let best_choice = _mm256_cmp_ps(ds, us, _CMP_LT_OS);
            let result = _mm256_blendv_ps(us, ds, best_choice);

            _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);
            for idx in j..cmp::min(j + 8, read.len()) {
                if m[i][idx - 1] + 1f32 < m[i][idx] {
                    m[i][idx] = m[i][idx - 1] + 1f32
                }
            }
        }
        for j in max_multiple + 1..read.len() {
            let l = m[i][j - 1] + 1f32;
            let u = m[i - 1][j] + 1f32;
            let d = m[i - 1][j - 1] + if read[j] == reference[i] { 0f32 } else { 1f32 };

            m[i][j] = [l, u, d].into_iter().reduce(f32::min).unwrap();
        }
    }

    m[reference.len() - 1][read.len() - 1]
}
/*
#[target_feature(enable = "sse2")]
pub unsafe fn exec_simd(read: &Vec<u8>, reference: &Vec<u8>) -> i32 {
    let mut m: Vec<Vec<i32>> = vec![vec![0; read.len()]; reference.len()];
    for i in 0..reference.len() {
        m[i][0] = i as i32
    }
    for j in 0..read.len() {
        m[0][j] = j as i32
    }
    let max_multiple = (read.len() / 4) * 4;
    let one_simd: __m128i = _mm_set1_epi32(1);
    for i in 1..reference.len() {
        for j in (1..max_multiple + 1).step_by(4) {
            let us = _mm_add_epi32(_mm_loadu_epi32(m[i - 1].get_unchecked(j)), one_simd);

            let read_simd = _mm_setr_epi32(
                read[j] as i32,
                read[j + 1] as i32,
                read[j + 2] as i32,
                read[j + 3] as i32,
            );
            let eq_diff_mask = _mm_cmpeq_epi32(_mm_set1_epi32(reference[i] as i32), read_simd);
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
            for idx in j..cmp::min(j + 4, read.len()) {
                if m[i][idx - 1] + 1 < m[i][idx] {
                    m[i][idx] = m[i][idx - 1] + 1
                }
            }
        }
        for j in max_multiple + 1..read.len() {
            let l = m[i][j - 1] + 1;
            let u = m[i - 1][j] + 1;
            let d = m[i - 1][j - 1] + if read[j] == reference[i] { 0 } else { 1 };

            m[i][j] = *[l, u, d].iter().min().unwrap();
        }
    }

    m[reference.len() - 1][read.len() - 1]
}
*/
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
    #[test]
    fn test_avx2_simd_version() {
        let s1 = "AAAAAAAAAAAAAAA"
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>();
        let s2 = "AAAAAAAAAAAAAAAA"
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>();

        let s3 = "ACAAAAAAAAAAAAAAA"
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>();
        let s4 = "AAAACAAAAAAAAAAAA"
            .chars()
            .map(|c| c as u8)
            .collect::<Vec<u8>>();
        unsafe {
            let ed1 = super::exec_simd_avx2(&s1, &s2);
            let ed2 = super::exec_simd_avx2(&s3, &s4);
            assert_eq!(ed1, 1f32);
            assert_eq!(ed2, 2f32);
        }
    }
}
