use crate::{bitfield_path as bf, utils};
use crate::{gaf_output, graph::LnzGraph};
use bitvec::prelude::*;
use std::arch::x86_64::*;
use std::cmp;
use std::collections::HashMap;

#[target_feature(enable = "avx2")]
pub unsafe fn exec_simd(
    // comments on simd instructions are in global_abpoa::exec_simd()
    read: &[char],
    seq_name: (&str, usize),
    graph: &LnzGraph,
    scores_matrix: &HashMap<(char, char), f32>,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) -> f32 {
    let mut m: Vec<Vec<f32>> = vec![vec![0f32; read.len()]; graph.lnz.len()];
    let mut path: Vec<Vec<f32>> = vec![vec![0f32; read.len()]; graph.lnz.len()];

    let max_multiple = if read.len() % 8 != 0 {
        (read.len() / 8) * 8
    } else {
        read.len() - 8
    };

    let d_move_simd = _mm256_set1_ps(0.1);
    let u_move_simd = _mm256_set1_ps(0.2);
    let mut best_row = 0;
    let mut best_col = 0;
    for i in 1..graph.lnz.len() - 1 {
        let us_update = _mm256_set1_ps(*scores_matrix.get(&(graph.lnz[i], '-')).unwrap());
        for j in (1..max_multiple + 1).step_by(8) {
            let ds_update = _mm256_set_ps(
                *scores_matrix.get(&(graph.lnz[i], read[j + 7])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 6])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 5])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 4])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 3])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 2])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j + 1])).unwrap(),
                *scores_matrix.get(&(graph.lnz[i], read[j])).unwrap(),
            );
            if !graph.nwp[i] {
                let us = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j)), us_update);

                let ds = _mm256_add_ps(_mm256_loadu_ps(m[i - 1].get_unchecked(j - 1)), ds_update);

                let best_choice = _mm256_cmp_ps(ds, us, _CMP_GT_OS);
                let result = _mm256_blendv_ps(us, ds, best_choice);

                _mm256_storeu_ps(m[i].get_unchecked_mut(j), result);

                let dir_result = _mm256_blendv_ps(u_move_simd, d_move_simd, best_choice);
                let path_update = _mm256_add_ps(_mm256_set1_ps((i - 1) as f32), dir_result);
                _mm256_storeu_ps(path[i].get_unchecked_mut(j), path_update);
            } else {
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

            // update with l for each one
            for idx in j..cmp::min(j + 8, read.len()) {
                let l = m[i][idx - 1] + scores_matrix.get(&(read[j], '-')).unwrap();
                if l > m[i][idx] {
                    m[i][idx] = l;
                    path[i][idx] = i as f32 + 0.3;
                }
                if m[i][idx] <= 0.0 {
                    m[i][idx] = 0.0;
                    path[i][idx] = 0.0;
                }
                if m[i][idx] >= m[best_row][best_col] {
                    best_row = i;
                    best_col = idx;
                }
            }
        }
        for j in max_multiple + 1..read.len() {
            if !graph.nwp[i] {
                let l = m[i][j - 1] + scores_matrix.get(&(read[j], '-')).unwrap();
                let u = m[i - 1][j] + scores_matrix.get(&(graph.lnz[i], '-')).unwrap();
                let d = m[i - 1][j - 1] + scores_matrix.get(&(graph.lnz[i], read[j])).unwrap();
                m[i][j] = [l, u, d].into_iter().reduce(f32::max).unwrap();
                if m[i][j] < 0.0 {
                    m[i][j] = 0.0;
                    path[i][j] = 0.0;
                } else if m[i][j] == d {
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
                u += scores_matrix.get(&(graph.lnz[i], '-')).unwrap();
                d += scores_matrix.get(&(read[j], graph.lnz[i])).unwrap();

                let l = m[i][j - 1] + scores_matrix.get(&(read[j], '-')).unwrap();

                m[i][j] = [l, u, d].into_iter().reduce(f32::max).unwrap();

                if m[i][j] == d {
                    path[i][j] = d_pred as f32 + 0.1;
                } else if m[i][j] == u {
                    path[i][j] = u_pred as f32 + 0.2;
                } else {
                    path[i][j] = i as f32 + 0.3;
                }
            }
            if m[i][j] >= m[best_row][best_col] {
                best_row = i;
                best_col = j;
            }
        }
    }
    if seq_name.1 != 0 {
        gaf_output::gaf_of_local_poa_simd(
            &path, read, seq_name, best_row, best_col, amb_mode, hofp,
        );
    }

    m[best_row][best_col]
}

pub fn exec(
    sequence: &[char],
    seq_name: (&str, usize),
    graph: &LnzGraph,
    scores_matrix: &HashMap<(char, char), i32>,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) -> i32 {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;

    let mut m = vec![vec![0; sequence.len()]; lnz.len()];
    let mut path = vec![vec![bitvec![u16, Msb0; 0; 32]; sequence.len()]; lnz.len()];
    let (mut best_row, mut best_col) = (0, 0);

    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, _) | (_, 0) => path[i][j] = bf::set_path_cell(0, 'O'),
                _ => {
                    let l = m[i][j - 1] + scores_matrix.get(&(sequence[j], '-')).unwrap();
                    let l_idx = i;

                    let mut d;
                    let d_idx;

                    let mut u;
                    let u_idx;
                    if !nodes_with_pred[i] {
                        d = m[i - 1][j - 1] + scores_matrix.get(&(sequence[j], lnz[i])).unwrap();
                        d_idx = i - 1;

                        u = m[i - 1][j] + scores_matrix.get(&('-', lnz[i])).unwrap();
                        u_idx = i - 1;
                    } else {
                        (d, d_idx) = get_best_d(&m, pred_hash.get(&i).unwrap(), j);
                        (u, u_idx) = get_best_u(&m, pred_hash.get(&i).unwrap(), j);
                        d += scores_matrix.get(&(sequence[j], lnz[i])).unwrap();
                        u += scores_matrix.get(&('-', lnz[i])).unwrap();
                    }
                    if d < 0 && l < 0 && u < 0 {
                        m[i][j] = 0;
                        path[i][j] = bf::set_path_cell(0, 'O');
                    } else {
                        let (best_val, mut dir) = utils::get_max_d_u_l(d, u, l);
                        if dir == 'D' && lnz[i] != sequence[j] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' | 'd' => bf::set_path_cell(d_idx, dir),
                            'U' => bf::set_path_cell(u_idx, dir),
                            _ => bf::set_path_cell(l_idx, dir),
                        }
                    }
                }
            }

            if m[i][j] > m[best_row][best_col] {
                best_row = i;
                best_col = j;
            }
        }
    }

    if seq_name.1 != 0 {
        gaf_output::gaf_of_local_poa(
            &path, sequence, seq_name, best_row, best_col, amb_mode, hofp,
        );
    }
    m[best_row][best_col]
}

fn get_best_d(m: &[Vec<i32>], p_arr: &[usize], j: usize) -> (i32, usize) {
    let mut d = 0;
    let mut d_idx = 0;
    let mut first = false;
    for p in p_arr {
        let current_d = m[*p][j - 1];
        if first {
            first = false;
            d = current_d;
            d_idx = *p;
        }
        if current_d > d {
            d = current_d;
            d_idx = *p;
        }
    }
    (d, d_idx)
}

fn get_best_u(m: &[Vec<i32>], p_arr: &[usize], j: usize) -> (i32, usize) {
    let mut u = 0;
    let mut u_idx = 0;
    let mut first = false;
    for p in p_arr {
        let current_u = m[*p][j];
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
    (u, u_idx)
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use bit_vec::BitVec;

    use crate::graph::LnzGraph;

    #[test]
    fn test_local_poa_consider_substrings() {
        let s = vec!['$', 'A', 'A', 'C', 'C', 'C', 'A', 'A'];

        let lnz = vec!['$', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'F'];
        let mut nwp = BitVec::from_elem(lnz.len(), false);
        nwp.set(1, true);
        nwp.set(8, true);
        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(8, vec![7]);
        let graph_struct = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        for c1 in ['A', 'C', 'G', '-'] {
            for c2 in ['A', 'C', 'G', '-'] {
                if c1 == c2 {
                    score_matrix.insert((c1, c2), 1);
                } else {
                    score_matrix.insert((c1, c2), -1);
                }
            }
        }
        let align_score = super::exec(
            &s,
            ("seq", 0),
            &graph_struct,
            &score_matrix,
            false,
            &HashMap::new(),
        );
        assert_eq!(align_score, 3);
    }

    #[test]
    fn local_poa_consider_best_predecessor() {
        let s = vec!['$', 'A', 'A', 'C', 'C', 'C', 'A', 'A'];

        let lnz = vec!['$', 'G', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'F'];
        let mut nwp = BitVec::from_elem(lnz.len(), false);
        nwp.set(1, true);
        nwp.set(6, true);
        nwp.set(9, true);
        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(6, vec![3]);
        pred_hash.insert(9, vec![8, 5]);
        let graph_struct = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        for c1 in ['A', 'C', 'G', '-'] {
            for c2 in ['A', 'C', 'G', '-'] {
                if c1 == c2 {
                    score_matrix.insert((c1, c2), 1);
                } else {
                    score_matrix.insert((c1, c2), -1);
                }
            }
        }
        let align_score = super::exec(
            &s,
            ("seq", 0),
            &graph_struct,
            &score_matrix,
            false,
            &HashMap::new(),
        );
        assert_eq!(align_score, 2);
    }
}
