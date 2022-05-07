use crate::graph::LnzGraph;
use crate::{bitfield_path as bf, gaf_output, utils};
use bitvec::prelude::*;
use std::{cmp::Ordering, collections::HashMap, vec};

pub fn exec(
    sequence: &[char],
    seq_name: (&str, usize),
    graph_struct: &LnzGraph,
    score_matrix: &HashMap<(char, char), i32>,
    bta: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>
) -> i32 {
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
                let (best_val, mut dir) = get_max_d_u_l(d, u, l);
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
        gaf_output::gaf_of_global_abpoa(
            &path,
            sequence,
            seq_name,
            //lnz,
            &ampl_for_row,
            last_row,
            last_col,
            amb_mode,
            hofp
        );
    }

    m[last_row][last_col]
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
#[inline]
fn get_max_d_u_l(d: i32, u: i32, l: i32) -> (i32, char) {
    let (best_val, dir) = match d.cmp(&u) {
        Ordering::Less => match u.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (u, 'U'),
        },
        _ => match d.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (d, 'D'),
        },
    };
    (best_val, dir)
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
            &HashMap::new()
        );

        assert_eq!(align, 4);
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
            &HashMap::new()
        );

        assert_eq!(align, 5);
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
            &HashMap::new()
        );

        assert_eq!(align, 5);
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
            &HashMap::new()
        );

        assert_eq!(align, 5);
    }
}
