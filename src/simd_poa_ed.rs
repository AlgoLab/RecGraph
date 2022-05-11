use std::{arch::x86_64::*, cmp};

use crate::graph::LnzGraph;

pub fn exec_no_simd(read: &Vec<u8>, graph: &LnzGraph) -> f32 {
    let mut m: Vec<Vec<f32>> = vec![vec![0f32; read.len()]; graph.lnz.len()];
    for i in 1..graph.lnz.len() - 1 {
        if !graph.nwp[i] {
            m[i][0] = m[i - 1][0] + 1f32;
        } else {
            let pred = graph.pred_hash.get(&i).unwrap();
            let best_p = pred.iter().min().unwrap();
            m[i][0] = m[*best_p][0] + 1f32;
        }
    }
    for j in 0..read.len() {
        m[0][j] = j as f32
    }
    for i in 1..graph.lnz.len() - 1 {
        for j in 1..read.len() {
            if !graph.nwp[i] {
                let l = m[i][j - 1] + 1f32;
                let u = m[i - 1][j] + 1f32;
                let d = m[i - 1][j - 1]
                    + if read[j] == graph.lnz[i] as u8 {
                        0f32
                    } else {
                        1f32
                    };

                m[i][j] = [l, u, d].into_iter().reduce(f32::min).unwrap();
            } else {
                let mut u = 0f32;
                let mut d = 0f32;
                let mut first = true;
                for p in graph.pred_hash.get(&i).unwrap() {
                    if first {
                        u = m[*p][j];
                        d = m[*p][j - 1];
                        first = false
                    }
                    if m[*p][j] < u {
                        u = m[*p][j];
                    }
                    if m[*p][j - 1] < d {
                        d = m[*p][j - 1];
                    }
                }
                u += 1f32;
                d += if read[j] == graph.lnz[i] as u8 {
                    0f32
                } else {
                    1f32
                };
                let l = m[i][j - 1] + 1f32;

                m[i][j] = [l, u, d].into_iter().reduce(f32::min).unwrap();
            }
        }
    }
    let mut best_result = 0f32;
    let mut first = true;
    for p in graph.pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        if first {
            best_result = m[*p][read.len() - 1];
            first = false;
        }
        if m[*p][read.len() - 1] < best_result {
            best_result = m[*p][read.len() - 1];
        }
    }
    best_result
}

#[target_feature(enable = "avx2")]
pub unsafe fn exec(read: &Vec<u8>, graph: &LnzGraph) -> f32 {
    let mut m: Vec<Vec<f32>> = vec![vec![0f32; read.len()]; graph.lnz.len()];

    for i in 1..graph.lnz.len() - 1 {
        if !graph.nwp[i] {
            m[i][0] = m[i - 1][0] + 1f32;
        } else {
            let pred = graph.pred_hash.get(&i).unwrap();
            let best_p = pred.iter().min().unwrap();
            m[i][0] = m[*best_p][0] + 1f32;
        }
    }
    for j in 0..read.len() {
        m[0][j] = j as f32
    }

    let max_multiple = (read.len() / 8) * 8;
    let read_f32 = &read[0..max_multiple + 1]
        .iter()
        .map(|c| *c as f32)
        .collect::<Vec<f32>>();
    let one_simd = _mm256_set1_ps(1.0);

    for i in 1..graph.lnz.len() - 1 {
        for j in (1..max_multiple + 1).step_by(8) {
            if !graph.nwp[i] {
                let us = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j)), one_simd);

                let eq_char = _mm256_cmp_ps(
                    _mm256_loadu_ps(read_f32.get_unchecked(j)), //read chars simd
                    _mm256_set1_ps(graph.lnz[i] as u8 as f32),  // reference char simd
                    _CMP_EQ_OS,
                );
                let neq_ds =
                    _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)), one_simd);
                let eq_ds = _mm256_loadu_ps(m[i - 1].get_unchecked(j - 1));
                let ds = _mm256_blendv_ps(neq_ds, eq_ds, eq_char);

                let best_choice = _mm256_cmp_ps(ds, us, _CMP_LT_OS);
                let result = _mm256_blendv_ps(us, ds, best_choice);

                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);
            } else {
                let preds = graph.pred_hash.get(&i).unwrap();
                let mut result = _mm256_set1_ps(0f32);
                let mut first = true;
                for p in preds {
                    let us = _mm256_add_ps(_mm256_loadu_ps(m[*p].get_unchecked(j)), one_simd);

                    let eq_char = _mm256_cmp_ps(
                        _mm256_loadu_ps(read_f32.get_unchecked(j)), //read chars simd
                        _mm256_set1_ps(graph.lnz[i] as u8 as f32),  // reference char simd
                        _CMP_EQ_OS,
                    );
                    let neq_ds =
                        _mm256_add_ps(_mm256_loadu_ps(m[*p].get_unchecked(j - 1)), one_simd);
                    let eq_ds = _mm256_loadu_ps(m[*p].get_unchecked(j - 1));
                    let ds = _mm256_blendv_ps(neq_ds, eq_ds, eq_char);

                    let best_choice = _mm256_cmp_ps(ds, us, _CMP_LT_OS);
                    let tmp_result = _mm256_blendv_ps(us, ds, best_choice);
                    if first {
                        result = tmp_result;
                        first = false;
                    } else {
                        let best_choice = _mm256_cmp_ps(tmp_result, result, _CMP_LT_OS);
                        let update_result = _mm256_blendv_ps(result, tmp_result, best_choice);
                        result = update_result
                    }
                }
                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);
            }
            // update with l for each one
            for idx in j..cmp::min(j + 8, read.len()) {
                if m[i][idx - 1] + 1f32 < m[i][idx] {
                    m[i][idx] = m[i][idx - 1] + 1f32
                }
            }
        }
        for j in max_multiple + 1..read.len() {
            if !graph.nwp[i] {
                let l = m[i][j - 1] + 1f32;
                let u = m[i - 1][j] + 1f32;
                let d = m[i - 1][j - 1]
                    + if read[j] == graph.lnz[i] as u8 {
                        0f32
                    } else {
                        1f32
                    };

                m[i][j] = [l, u, d].into_iter().reduce(f32::min).unwrap();
            } else {
                let mut u = 0f32;
                let mut d = 0f32;
                let mut first = false;
                for p in graph.pred_hash.get(&i).unwrap() {
                    if first {
                        u = m[*p][j];
                        d = m[*p][j - 1];
                        first = false
                    }
                    if m[*p][j] < u {
                        u = m[*p][j];
                    }
                    if m[*p][j - 1] < d {
                        d = m[*p][j - 1];
                    }
                }
                u += 1f32;
                d += if read[j] == graph.lnz[i] as u8 {
                    0f32
                } else {
                    1f32
                };
                let l = m[i][j - 1] + 1f32;

                m[i][j] = [l, u, d].into_iter().reduce(f32::min).unwrap();
            }
        }
    }
    let mut best_result = 0f32;
    let mut first = true;
    for p in graph.pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        if first {
            best_result = m[*p][read.len() - 1];
            first = false;
        }
        if m[*p][read.len() - 1] < best_result {
            best_result = m[*p][read.len() - 1];
        }
    }
    best_result
}
