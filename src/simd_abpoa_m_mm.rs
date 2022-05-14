use std::{arch::x86_64::*, cmp};

use crate::{graph::LnzGraph, utils};
//TODO: no match(i, j) set i = 0 and j = 0 before;
#[target_feature(enable = "avx2")]
pub unsafe fn exec(
    read: &Vec<u8>,
    graph: &LnzGraph,
    score_match: f32,
    score_mis: f32,
    bta: usize,
    r_values: &Vec<usize>,
) -> f32 {
    let min_score = 2.0 * read.len() as f32 * score_mis;
    let mut m: Vec<Vec<f32>> = vec![vec![min_score; read.len()]; graph.lnz.len()];
    let mut path: Vec<Vec<f32>> = vec![vec![-1f32; read.len()]; graph.lnz.len()];

    let read_f32 = &read.iter().map(|c| *c as f32).collect::<Vec<f32>>();

    let mismatch_simd = _mm256_set1_ps(score_mis);
    let match_simd = _mm256_set1_ps(score_match);
    let d_move_simd = _mm256_set1_ps(0.1);
    let u_move_simd = _mm256_set1_ps(0.2);

    let mut best_scoring_pos = vec![0; graph.lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); graph.lnz.len()];

    m[0][0] = 0f32;
    path[0][0] = 0f32;

    for i in 1..graph.lnz.len() - 1 {
        if !graph.nwp[i] {
            m[i][0] = m[i - 1][0] + score_mis;
            path[i][0] = (i - 1) as f32 + 0.2;
        } else {
            let pred = graph.pred_hash.get(&i).unwrap();
            let best_p = pred.iter().min().unwrap();
            m[i][0] = m[*best_p][0] + score_mis;
            path[i][0] = *best_p as f32 + 0.2;
        }
    }
    let p_arr: Vec<usize> = vec![];
    let (left, right) = utils::set_ampl_for_row(
        0,
        &p_arr,
        r_values[0],
        &best_scoring_pos,
        read.len(),
        bta,
        true,
    );
    for j in 1..right {
        m[0][j] = m[0][j - 1] + score_mis;
        path[0][j] = 0.3;
    }
    ampl_for_row[0] = (left, right);

    for i in 1..graph.lnz.len() - 1 {
        let mut p_arr = &vec![];
        if graph.nwp[i] {
            p_arr = graph.pred_hash.get(&i).unwrap()
        }
        let (left, right) = utils::set_ampl_for_row(
            i,
            p_arr,
            r_values[i],
            &best_scoring_pos,
            read.len(),
            bta,
            true,
        );
        ampl_for_row[i] = (left, right);
        let mut best_col = left;
        let start = if left == 0 { 1 } else { left };
        let end = if right == read.len() {
            ((right - start) / 8) * 8 + start
        } else {
            right
        };
        for j in (start..end).step_by(8) {
            if !graph.nwp[i] {
                let us = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j)), mismatch_simd);

                let eq_char = _mm256_cmp_ps(
                    _mm256_loadu_ps(read_f32.get_unchecked(j)), //read chars simd
                    _mm256_set1_ps(graph.lnz[i] as u8 as f32),  // reference char simd
                    _CMP_EQ_OS,
                );
                let neq_ds = _mm256_add_ps(
                    _mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)),
                    mismatch_simd,
                );
                let eq_ds =
                    _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)), match_simd);
                let ds = _mm256_blendv_ps(neq_ds, eq_ds, eq_char);

                let best_choice = _mm256_cmp_ps(ds, us, _CMP_GT_OS);
                let result = _mm256_blendv_ps(us, ds, best_choice);

                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);

                let dir_result = _mm256_blendv_ps(u_move_simd, d_move_simd, best_choice);
                let path_update = _mm256_add_ps(_mm256_set1_ps((i - 1) as f32), dir_result);
                _mm256_storeu_ps(path[i].get_unchecked_mut(j), path_update);
            } else {
                let preds = graph.pred_hash.get(&i).unwrap();
                let mut best_ds = _mm256_set1_ps(0f32);
                let mut best_us = _mm256_set1_ps(0f32);
                let mut pred_best_us = _mm256_set1_ps(0f32);
                let mut pred_best_ds = _mm256_set1_ps(0f32);
                let mut first = true;
                for p in preds {
                    let us = _mm256_loadu_ps(m[*p].get_unchecked(j));

                    let ds = _mm256_loadu_ps(m[*p].get_unchecked(j - 1));
                    let pred_simd = _mm256_set1_ps(*p as f32);
                    if first {
                        first = false;
                        best_us = us;
                        best_ds = ds;
                        pred_best_us = pred_simd;
                        pred_best_ds = pred_simd;
                    } else {
                        let best_us_choices = _mm256_cmp_ps(us, best_us, _CMP_GT_OS);
                        best_us = _mm256_blendv_ps(best_us, us, best_us_choices);
                        pred_best_us = _mm256_blendv_ps(pred_best_us, pred_simd, best_us_choices);

                        let best_ds_choices = _mm256_cmp_ps(ds, best_ds, _CMP_GT_OS);
                        best_ds = _mm256_blendv_ps(best_ds, ds, best_ds_choices);
                        pred_best_ds = _mm256_blendv_ps(pred_best_ds, pred_simd, best_ds_choices);
                    }
                }
                best_us = _mm256_add_ps(best_us, mismatch_simd);

                let eq_char = _mm256_cmp_ps(
                    _mm256_loadu_ps(read_f32.get_unchecked(j)), //read chars simd
                    _mm256_set1_ps(graph.lnz[i] as u8 as f32),  // reference char simd
                    _CMP_EQ_OS,
                );
                let neq_ds = _mm256_add_ps(best_ds, mismatch_simd);
                best_ds = _mm256_add_ps(best_ds, match_simd);
                best_ds = _mm256_blendv_ps(neq_ds, best_ds, eq_char);

                let best_choice = _mm256_cmp_ps(best_ds, best_us, _CMP_GT_OS);
                let result = _mm256_blendv_ps(best_us, best_ds, best_choice);
                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);

                pred_best_ds = _mm256_add_ps(pred_best_ds, d_move_simd);
                pred_best_us = _mm256_add_ps(pred_best_us, u_move_simd);

                let dir_result = _mm256_blendv_ps(pred_best_us, pred_best_ds, best_choice);
                _mm256_storeu_ps(path[i].get_unchecked_mut(j), dir_result);
            }
            for idx in j..cmp::min(j + 8, read.len()) {
                if m[i][idx - 1] + score_mis > m[i][idx] {
                    m[i][idx] = m[i][idx - 1] + score_mis;
                    path[i][idx] = i as f32 + 0.3;
                }
            }

            if m[i][j] >= m[i][best_col] {
                best_col = j
            }
        }
        // set last position without simd
        if end != right {
            for j in end..read.len() {
                if !graph.nwp[i] {
                    let l = m[i][j - 1] + score_mis;
                    let u = m[i - 1][j] + score_mis;
                    let d = m[i - 1][j - 1]
                        + if read[j] == graph.lnz[i] as u8 {
                            score_match
                        } else {
                            score_mis
                        };

                    m[i][j] = [l, u, d].into_iter().reduce(f32::max).unwrap();
                    if m[i][j] == d {
                        path[i][j] = (i - 1) as f32 + 0.1;
                    } else if m[i][j] == u {
                        path[i][j] = (i - 1) as f32 + 0.2;
                    } else {
                        path[i][j] = i as f32 + 0.3;
                    }
                } else {
                    let mut u = 0f32;
                    let mut u_pred = 0;
                    let mut d = 0f32;
                    let mut d_pred = 0;
                    let mut first = true;
                    for p in graph.pred_hash.get(&i).unwrap() {
                        if first {
                            u = m[*p][j];
                            d = m[*p][j - 1];
                            u_pred = *p;
                            d_pred = *p;
                            first = false
                        }
                        if m[*p][j] > u {
                            u = m[*p][j];
                            u_pred = *p;
                        }
                        if m[*p][j - 1] > d {
                            d = m[*p][j - 1];
                            d_pred = *p;
                        }
                    }
                    u += score_mis;
                    d += if read[j] == graph.lnz[i] as u8 {
                        score_match
                    } else {
                        score_mis
                    };
                    let l = m[i][j - 1] + score_mis;

                    m[i][j] = [l, u, d].into_iter().reduce(f32::max).unwrap();

                    if m[i][j] == d {
                        path[i][j] = d_pred as f32 + 0.1;
                    } else if m[i][j] == u {
                        path[i][j] = u_pred as f32 + 0.2;
                    } else {
                        path[i][j] = i as f32 + 0.3;
                    }
                }
                if m[i][j] >= m[i][best_col] {
                    best_col = j
                }
            }
        }
        best_scoring_pos[i] = best_col;
    }
    let mut best_result = 0f32;
    let mut first = true;
    for p in graph.pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        if first {
            best_result = m[*p][read.len() - 1];
            first = false;
        }
        if m[*p][read.len() - 1] > best_result {
            best_result = m[*p][read.len() - 1];
        }
    }
    //rebuild_path(&path);
    best_result
}

fn rebuild_path(path: &Vec<Vec<f32>>) {
    let mut row = path.len() - 2;
    let mut col = path[row].len() - 1;
    while path[row][col] != 0.0 {
        let val = path[row][col];
        let pred = val as usize;
        let dir = val - (val as i32) as f32;
        row = pred;
        col = match (dir * 11f32) as i32 {
            1 => col - 1,
            2 => col,
            3 => col - 1,
            _ => {
                panic!();
            }
        }
    }
}
