use crate::gaf_output::GAFStruct;
use crate::graph::LnzGraph;
use crate::{bitfield_path as bf, gaf_output, utils};
use bitvec::prelude::*;
use std::arch::x86_64::*;
use std::{collections::HashMap, vec};

/// if avx2 features are detected this is the chosen algorithm
#[target_feature(enable = "avx2")]
pub unsafe fn exec_simd(
    read: &[char],
    seq_name: (&str, usize),
    graph: &LnzGraph,
    score_matrix: &HashMap<(char, char), f32>,
    bta: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
    r_values: &Vec<usize>,
) -> (i32, Option<GAFStruct>) {
    let min_score = 2.0 * read.len() as f32 * score_matrix.get(&(read[1], '-')).unwrap();
    let mut m: Vec<Vec<f32>> = vec![vec![min_score; read.len()]; graph.lnz.len()];
    let mut path: Vec<Vec<f32>> = vec![vec![-1f32; read.len()]; graph.lnz.len()];
    // value in path cells use integer part for predecessor, decimal for the choice
    // diagonal -> 0.1
    // upper -> 0.2
    // left -> 0.3
    let d_move_simd = _mm256_set1_ps(0.1);
    let u_move_simd = _mm256_set1_ps(0.2);

    let mut best_scoring_pos = vec![0; graph.lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); graph.lnz.len()];

    //base cases
    m[0][0] = 0f32;
    path[0][0] = 0.00;
    for i in 1..graph.lnz.len() - 1 {
        if !graph.nwp[i] {
            m[i][0] = m[i - 1][0] + score_matrix.get(&(graph.lnz[i], '-')).unwrap();
            path[i][0] = (i - 1) as f32 + 0.2;
        } else {
            let pred = graph.pred_hash.get(&i).unwrap();
            let best_p = pred.iter().min().unwrap();
            m[i][0] = m[*best_p][0] + score_matrix.get(&(graph.lnz[i], '-')).unwrap();
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
        m[0][j] = m[0][j - 1] + score_matrix.get(&(read[j], '-')).unwrap();
        path[0][j] = 0.3;
    }
    ampl_for_row[0] = (left, right);

    for i in 1..graph.lnz.len() - 1 {
        // set ampl for current row
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

        // us_update is set with the mismatch score of lnz[i], indel
        let us_update = _mm256_set1_ps(*score_matrix.get(&(graph.lnz[i], '-')).unwrap());
        for j in (start..end).step_by(8) {
            // ds_update is set with the match/mismatch score of lnz[i], read[k] for k in j..j+8
            let ds_update = _mm256_set_ps(
                *score_matrix.get(&(graph.lnz[i], read[j + 7])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 6])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 5])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 4])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 3])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 2])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j + 1])).unwrap(),
                *score_matrix.get(&(graph.lnz[i], read[j])).unwrap(),
            );
            if !graph.nwp[i] {
                let us = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j)), us_update);

                let ds = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)), ds_update);

                // best_choice is a mask with true if ds[a] > us[a], false otherwise
                let best_choice = _mm256_cmp_ps(ds, us, _CMP_GT_OS);

                // in result ds[a] is saved if best_choice[a] is true, otherwise us[a]
                let result = _mm256_blendv_ps(us, ds, best_choice);

                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);

                //path is updated with the same choice taken before
                let dir_result = _mm256_blendv_ps(u_move_simd, d_move_simd, best_choice);
                let path_update = _mm256_add_ps(_mm256_set1_ps((i - 1) as f32), dir_result);
                _mm256_storeu_ps(path[i].get_unchecked_mut(j), path_update);
            } else {
                // similar to !graph.nwp[i], but in this case iteration on every predecessor is needed

                let preds = graph.pred_hash.get(&i).unwrap();
                let mut best_us = _mm256_loadu_ps(m[preds[0]].get_unchecked(j));
                let mut best_ds = _mm256_loadu_ps(m[preds[0]].get_unchecked(j - 1));
                let mut pred_best_us = _mm256_set1_ps(preds[0] as f32);
                let mut pred_best_ds = _mm256_set1_ps(preds[0] as f32);

                for p in preds[1..].iter() {
                    let us = _mm256_loadu_ps(m[*p].get_unchecked(j));
                    let ds = _mm256_loadu_ps(m[*p].get_unchecked(j - 1));
                    let pred_simd = _mm256_set1_ps(*p as f32);

                    let best_us_choices = _mm256_cmp_ps(us, best_us, _CMP_GT_OS);
                    best_us = _mm256_blendv_ps(best_us, us, best_us_choices);
                    pred_best_us = _mm256_blendv_ps(pred_best_us, pred_simd, best_us_choices);

                    let best_ds_choices = _mm256_cmp_ps(ds, best_ds, _CMP_GT_OS);
                    best_ds = _mm256_blendv_ps(best_ds, ds, best_ds_choices);
                    pred_best_ds = _mm256_blendv_ps(pred_best_ds, pred_simd, best_ds_choices);
                }
                best_us = _mm256_add_ps(best_us, us_update);

                best_ds = _mm256_add_ps(best_ds, ds_update);

                let best_choice = _mm256_cmp_ps(best_ds, best_us, _CMP_GT_OS);
                let result = _mm256_blendv_ps(best_us, best_ds, best_choice);

                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);

                pred_best_ds = _mm256_add_ps(pred_best_ds, d_move_simd);
                pred_best_us = _mm256_add_ps(pred_best_us, u_move_simd);

                let dir_result = _mm256_blendv_ps(pred_best_us, pred_best_ds, best_choice);
                _mm256_storeu_ps(path[i].get_unchecked_mut(j), dir_result);
            }
            //update with left value if it is better than upper and diagonal
            for idx in j..j + 8 {
                let l = m[i][idx - 1] + score_matrix.get(&(read[j], '-')).unwrap();
                if l > m[i][idx] {
                    m[i][idx] = l;
                    path[i][idx] = i as f32 + 0.3;
                }
                if m[i][idx] >= m[i][best_col] {
                    best_col = idx
                }
            }
        }
        // set last positions without simd (positions with simd must be multiple of 8)
        if end < right {
            for j in end..right {
                if !graph.nwp[i] {
                    let l = m[i][j - 1] + score_matrix.get(&(read[j], '-')).unwrap();
                    let u = m[i - 1][j] + score_matrix.get(&(graph.lnz[i], '-')).unwrap();
                    let d = m[i - 1][j - 1] + score_matrix.get(&(graph.lnz[i], read[j])).unwrap();
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
                    u += score_matrix.get(&(graph.lnz[i], '-')).unwrap();
                    d += score_matrix.get(&(read[j], graph.lnz[i])).unwrap();

                    let l = m[i][j - 1] + score_matrix.get(&(read[j], '-')).unwrap();

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
    let mut last_row = 0;
    for p in graph.pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        if first {
            best_result = m[*p][read.len() - 1];
            last_row = *p;
            first = false;
        }
        if m[*p][read.len() - 1] > best_result {
            best_result = m[*p][read.len() - 1];
            last_row = *p
        }
    }
    if seq_name.1 != 0 {
        let gaf_struct = gaf_output::gaf_of_global_abpoa_simd(
            &path,
            &read,
            seq_name,
            last_row,
            read.len() - 1,
            amb_mode,
            hofp,
            &graph.lnz,
            best_result,
        );
        (best_result as i32, Some(gaf_struct))
    } else {
        (best_result as i32, None)
    }
}

/// adaptive banded POA without SIMD instructions
pub fn exec(
    // comments on normal version in gap_global_abpoa::exec()
    sequence: &[char],
    seq_name: (&str, usize),
    graph_struct: &LnzGraph,
    score_matrix: &HashMap<(char, char), i32>,
    bta: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) -> (i32, Option<GAFStruct>) {
    let lnz = &graph_struct.lnz;
    let nodes_w_pred = &graph_struct.nwp;
    let pred_hash = &graph_struct.pred_hash;

    let r_values = utils::set_r_values(nodes_w_pred, pred_hash, lnz.len());
    let mut best_scoring_pos = vec![0; lnz.len()];

    let mut m = vec![vec![]; lnz.len()];
    let mut path = vec![vec![]; lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); lnz.len()];

    let simd_version = false;
    for i in 0..lnz.len() - 1 {
        let mut p_arr = &vec![];
        if nodes_w_pred[i] {
            p_arr = pred_hash.get(&i).unwrap()
        }
        let (left, right) = utils::set_ampl_for_row(
            i,
            p_arr,
            r_values[i],
            &best_scoring_pos,
            sequence.len(),
            bta,
            simd_version,
        );
        ampl_for_row[i] = (left, right);
        let mut best_val_pos: usize = 0;
        m[i] = vec![0; right - left];
        let bv = bitvec![u16, Msb0; 0; 32];
        path[i] = vec![bv; right - left];
        for j in 0..right - left {
            if i == 0 && j == 0 {
                m[i][j] = 0;
                path[i][j] = bf::set_path_cell(0, 'O');
            } else if i == 0 {
                //only left
                m[i][j] = m[i][j - 1] + score_matrix.get(&('-', sequence[j + left])).unwrap();
                path[i][j] = bf::set_path_cell(i, 'L');
            } else if j == 0 && left == 0 {
                // only upper
                let best_p = if !nodes_w_pred[i] {
                    i - 1
                } else {
                    *pred_hash.get(&i).unwrap().iter().min().unwrap()
                };
                m[i][j] = m[best_p][j] + score_matrix.get(&('-', lnz[i])).unwrap();
                path[i][j] = bf::set_path_cell(best_p, 'U');
            } else {
                let mut p_arr = &vec![i - 1];
                if nodes_w_pred[i] {
                    p_arr = pred_hash.get(&i).unwrap();
                }
                // get best l
                let l;
                let l_pred;
                match get_best_l(&m, i, j) {
                    Some((l_score, idx)) => {
                        l = l_score + score_matrix.get(&(sequence[j + left], '-')).unwrap();
                        l_pred = idx;
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };
                        l = score_matrix.get(&(sequence[j + left], '-')).unwrap()
                            * (i + left + j) as i32;
                        l_pred = best_p;
                    }
                }
                //get best u
                let u;
                let u_pred;

                match get_best_u(p_arr, &m, &ampl_for_row, i, j) {
                    Some((u_score, idx)) => {
                        u = u_score + score_matrix.get(&(lnz[i], '-')).unwrap();
                        u_pred = idx;
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };
                        u = score_matrix.get(&(lnz[i], '-')).unwrap() * (i + left + j) as i32;
                        u_pred = best_p;
                    }
                }
                //get best d
                let d;
                let d_pred;
                match get_best_d(p_arr, &m, &ampl_for_row, i, j) {
                    Some((d_score, idx)) => {
                        d = d_score + score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        d_pred = idx;
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };
                        d = score_matrix.get(&(lnz[i], '-')).unwrap() * (i + left) as i32;
                        d_pred = best_p;
                    }
                }
                let (best_val, mut dir) = utils::get_max_d_u_l(d, u, l);
                if dir == 'D' && sequence[j + left] != lnz[i] {
                    dir = 'd'
                }
                m[i][j] = best_val;
                path[i][j] = match dir {
                    'D' => bf::set_path_cell(d_pred, 'D'),
                    'd' => bf::set_path_cell(d_pred, 'd'),
                    'U' => bf::set_path_cell(u_pred, 'U'),
                    _ => bf::set_path_cell(l_pred, 'L'),
                };
            }
            if m[i][j] >= m[i][best_val_pos] {
                best_val_pos = j
            }
        }
        best_scoring_pos[i] = best_val_pos + left;
    }
    let mut last_row = m.len() - 2;
    let mut last_col = m[last_row].len() - 1;
    for p in pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        let tmp_last_col = (ampl_for_row[*p].1 - ampl_for_row[*p].0) - 1;
        if m[*p][tmp_last_col] > m[last_row][last_col] {
            last_row = *p;
            last_col = tmp_last_col;
        }
    }
    let check = band_ampl_enough(&path, &ampl_for_row, sequence.len(), last_row, last_col);
    if !check {
        println!("Band length probably too short, maybe try with larger b and f");
    }

    if seq_name.1 != 0 {
        let gaf_struct = gaf_output::gaf_of_global_abpoa(
            &path,
            sequence,
            seq_name,
            //lnz,
            &ampl_for_row,
            last_row,
            last_col,
            amb_mode,
            hofp,
        );
        (m[last_row][last_col], Some(gaf_struct))
    } else {
        (m[last_row][last_col], None)
    }
}
fn band_ampl_enough(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    ampl_for_row: &[(usize, usize)],
    sequence_len: usize,
    start_row: usize,
    start_col: usize,
) -> bool {
    let mut i = start_row;
    let mut j = start_col;
    while bf::dir_from_bitvec(&path[i][j]) != 'O' {
        let (left, right) = ampl_for_row[i];
        if i == 0 || j == 0 && left == 0 {
            return true;
        }
        if (j == 0 && left != 0) || (j == right - left - 1 && right != sequence_len) {
            return false;
        }
        let curr_bv = &path[i][j];
        let pred = bf::pred_from_bitvec(curr_bv);
        let left_p = ampl_for_row[pred].0;
        let j_pos = if left_p < left {
            j + (left - left_p)
        } else {
            j - (left_p - left)
        };
        let dir = bf::dir_from_bitvec(curr_bv);
        match dir {
            'D' => {
                j = j_pos - 1;
                i = pred;
            }
            'd' => {
                j = j_pos - 1;
                i = pred;
            }
            'L' => {
                j -= 1;
            }
            'U' => {
                i = pred;
                j = j_pos;
            }
            _ => {
                panic!();
            }
        }
    }
    true
}
#[inline]
fn get_best_l(m: &[Vec<i32>], i: usize, j: usize) -> Option<(i32, usize)> {
    if j > 0 {
        Some((m[i][j - 1], i))
    } else {
        None
    }
}

#[inline]
fn get_best_d(
    p_arr: &[usize],
    m: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut d = 0;
    let mut d_idx = 0;
    let mut first = true;
    let left_i = ampl_for_row[i].0;

    for p in p_arr.iter() {
        let left_p = ampl_for_row[*p].0;
        if j + left_i > ampl_for_row[*p].0 && j + left_i <= ampl_for_row[*p].1 {
            let j_pos = if left_p < left_i {
                j + (left_i - left_p)
            } else {
                j - (left_p - left_i)
            };
            let curr_d = m[*p][j_pos - 1];
            if first {
                d = curr_d;
                d_idx = *p;
                first = false;
            }
            if curr_d > d {
                d = curr_d;
                d_idx = *p;
            }
        }
    }

    if first {
        // j is too far to be aligned with current predecessors
        None
    } else {
        Some((d, d_idx))
    }
}

#[inline]
fn get_best_u(
    p_arr: &[usize],
    m: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut u = 0;
    let mut u_idx = 0;
    let left_i = ampl_for_row[i].0;
    let mut first = true;
    for p in p_arr.iter() {
        let left_p = ampl_for_row[*p].0;
        if j + left_i >= ampl_for_row[*p].0 && j + left_i < ampl_for_row[*p].1 {
            let j_pos = if left_p < left_i {
                j + (left_i - left_p)
            } else {
                j - (left_p - left_i)
            };
            let current_u = m[*p][j_pos as usize];
            if first {
                first = false;
                u = current_u;
                u_idx = *p;
            }
            if current_u > u {
                u = current_u;
                u_idx = *p;
            }
        }
    }

    if first {
        None
    } else {
        Some((u, u_idx))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use bit_vec::BitVec;

    use crate::graph::LnzGraph;

    #[test]
    fn test1() {
        let s = vec!['$', 'A', 'A', 'A', 'A'];

        let lnz = vec!['$', 'A', 'A', 'A', 'A', 'F'];

        let mut nwp = BitVec::from_elem(6, false);
        nwp.set(1, true);
        nwp.set(5, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(5, vec![4]);

        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('A', '-'), -1);
        score_matrix.insert(('-', 'A'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            100,
            false,
            &HashMap::new(),
        );

        assert_eq!(align.0, 4);
    }
    #[test]
    fn test2() {
        let s = vec!['$', 'A', 'A', 'C', 'A', 'A'];

        let lnz = vec!['$', 'A', 'A', 'C', 'A', 'A', 'A', 'F'];

        let mut nwp = BitVec::from_elem(8, false);
        nwp.set(1, true);
        nwp.set(3, true);
        nwp.set(4, true);
        nwp.set(5, true);
        nwp.set(7, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(3, vec![2]);
        pred_hash.insert(4, vec![2]);
        pred_hash.insert(5, vec![3, 4]);
        pred_hash.insert(7, vec![6]);
        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('A', '-'), -1);
        score_matrix.insert(('-', 'A'), -1);
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('-', 'C'), -1);
        score_matrix.insert(('C', '-'), -1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            4,
            false,
            &HashMap::new(),
        );

        assert_eq!(align.0, 5);
    }
    #[test]
    fn multiple_starts() {
        let s = vec!['$', 'C', 'A', 'C', 'A', 'A'];

        let lnz = vec!['$', 'A', 'C', 'A', 'C', 'C', 'A', 'A', 'F'];

        let mut nwp = BitVec::from_elem(9, false);
        nwp.set(1, true);
        nwp.set(2, true);
        nwp.set(3, true);
        nwp.set(4, true);
        nwp.set(5, true);
        nwp.set(8, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(2, vec![0]);
        pred_hash.insert(3, vec![1, 2]);
        pred_hash.insert(4, vec![1, 2]);
        pred_hash.insert(5, vec![3, 4]);
        pred_hash.insert(8, vec![7]);
        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('A', '-'), -1);
        score_matrix.insert(('-', 'A'), -1);
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('-', 'C'), -1);
        score_matrix.insert(('C', '-'), -1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            4,
            false,
            &HashMap::new(),
        );

        assert_eq!(align.0, 5);
    }

    #[test]
    fn multiple_ends() {
        let s = vec!['$', 'C', 'A', 'C', 'A', 'A'];

        let lnz = vec!['$', 'A', 'C', 'A', 'C', 'C', 'A', 'A', 'C', 'F'];

        let mut nwp = BitVec::from_elem(10, false);
        nwp.set(1, true);
        nwp.set(2, true);
        nwp.set(3, true);
        nwp.set(4, true);
        nwp.set(5, true);
        nwp.set(7, true);
        nwp.set(8, true);
        nwp.set(9, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(2, vec![0]);
        pred_hash.insert(3, vec![1, 2]);
        pred_hash.insert(4, vec![1, 2]);
        pred_hash.insert(5, vec![3, 4]);
        pred_hash.insert(7, vec![6]);
        pred_hash.insert(8, vec![6]);
        pred_hash.insert(9, vec![7, 8]);
        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('A', '-'), -1);
        score_matrix.insert(('-', 'A'), -1);
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('-', 'C'), -1);
        score_matrix.insert(('C', '-'), -1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            4,
            false,
            &HashMap::new(),
        );

        assert_eq!(align.0, 5);
    }
}
