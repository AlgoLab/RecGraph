use crate::{bitfield_path as bf, gaf_output, graph::LnzGraph};
use bitvec::prelude::*;
use std::{
    cmp::{self, Ordering},
    collections::HashMap,
    vec,
};

pub fn exec(
    sequence: &[char],
    seq_name: (&str, usize),
    graph_struct: &LnzGraph,
    score_matrix: &HashMap<(char, char), i32>,
    o: i32,
    e: i32,
    bta: usize,
    file_path: &str,
    amb_mode: bool,
) -> i32 {
    let lnz = &graph_struct.lnz;
    let nodes_w_pred = &graph_struct.nwp;
    let pred_hash = &graph_struct.pred_hash;

    let mut m = vec![vec![]; lnz.len()]; // best alignment
    let mut x = vec![vec![]; lnz.len()]; //best alignment final gap in graph
    let mut y = vec![vec![]; lnz.len()]; // best alignment final gap in sequence

    let mut path = vec![vec![]; lnz.len()];
    let mut path_x = vec![vec![]; lnz.len()];
    let mut path_y = vec![vec![]; lnz.len()];

    let r_values = set_r_values(nodes_w_pred, pred_hash, lnz.len());
    let mut best_scoring_pos = vec![0; lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); lnz.len()];

    // TODO: set true if using simd
    let simd_version = false;
    for i in 0..lnz.len() - 1 {
        let mut p_arr = &vec![];
        if nodes_w_pred[i] {
            p_arr = pred_hash.get(&i).unwrap()
        }
        let (left, right) = set_ampl_for_row(
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
        x[i] = vec![0; right - left];
        y[i] = vec![0; right - left];

        path[i] = vec![bitvec![u16, Msb0; 0; 32]; right - left];
        path_x[i] = vec![bitvec![u16, Msb0; 0; 32]; right - left];
        path_y[i] = vec![bitvec![u16, Msb0; 0; 32]; right - left];

        for j in 0..right - left {
            if i == 0 && j == 0 {
                m[i][j] = 0;
                path[i][j] = bf::set_path_cell(0, 'O');
            } else if i == 0 {
                // set y
                y[i][j] = o + e * (j + left) as i32;
                // set m
                m[i][j] = y[i][j];
                path[i][j] = bf::set_path_cell(i, 'L');
            } else if j == 0 && left == 0 {
                // set x

                let best_p = if !nodes_w_pred[i] {
                    i - 1
                } else {
                    *pred_hash.get(&i).unwrap().iter().min().unwrap()
                };

                x[i][j] = o + e * (best_p + 1) as i32;

                // set m
                m[i][j] = x[i][j];
                path[i][j] = bf::set_path_cell(best_p, 'U');
            } else {
                let mut p_arr = &vec![i - 1];
                if nodes_w_pred[i] {
                    p_arr = pred_hash.get(&i).unwrap();
                }

                // try set and get l from x (j > left, if j == left same as j = 0)
                let l_score_idx = get_best_l(&m, &x, i, j, o);
                let l_pred;
                match l_score_idx {
                    Some((l, idx, from_m)) => {
                        x[i][j] = l + e;
                        l_pred = idx;
                        if !from_m {
                            path_x[i][j] = bf::set_path_cell(i, 'X');
                        }
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };

                        x[i][j] = 2 * o + e * (best_p + 1) as i32 + e * (j + left) as i32;
                        l_pred = best_p;
                    }
                }
                // try set and get u from y (pred_left < j < pred_right else None)
                let u_score_idx = get_best_u(p_arr, &m, &y, &ampl_for_row, i, j, o);
                let u_pred;
                match u_score_idx {
                    Some((u, idx, from_m)) => {
                        y[i][j] = u + e;
                        u_pred = idx;
                        if !from_m {
                            path_y[i][j] = bf::set_path_cell(idx, 'Y');
                        }
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };

                        y[i][j] = 2 * o + e * (best_p + 1) as i32 + e * (j + left) as i32;
                        u_pred = best_p;
                    }
                }
                // try get d from m (pred_left < j < pred_right else None)
                let d_score_idx = get_best_d(p_arr, &m, &ampl_for_row, i, j);
                match d_score_idx {
                    Some((mut d, d_idx)) => {
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        let l = x[i][j];
                        let u = y[i][j];
                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => match l.cmp(&u) {
                                Ordering::Less => {
                                    if u_pred == 0 {
                                        path[i][j] = bf::set_path_cell(u_pred, 'u');
                                    } else {
                                        path[i][j] = bf::set_path_cell(u_pred, 'U');
                                    }
                                    u
                                }
                                _ => {
                                    path[i][j] = bf::set_path_cell(l_pred, 'L');
                                    l
                                }
                            },
                            _ => match d.cmp(&u) {
                                Ordering::Less => {
                                    path[i][j] = bf::set_path_cell(u_pred, 'U');
                                    u
                                }
                                _ => {
                                    if lnz[i] == sequence[j + left] {
                                        path[i][j] = bf::set_path_cell(d_idx, 'D');
                                    } else {
                                        path[i][j] = bf::set_path_cell(d_idx, 'd');
                                    }
                                    d
                                }
                            },
                        }
                    }
                    _ => {
                        let l = x[i][j];
                        let u = y[i][j];
                        m[i][j] = match l.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = bf::set_path_cell(u_pred, 'U');
                                u
                            }
                            _ => {
                                path[i][j] = bf::set_path_cell(l_pred, 'L');
                                l
                            }
                        }
                    }
                }
                // set m with best d, u, l != from None
            }
            //update best scoring position
            if m[i][j] >= m[i][best_val_pos] {
                best_val_pos = j
            }
        }
        // set best scoring position for current row
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
    let best_value = m[last_row][last_col];
    let check = band_ampl_enough(
        &path,
        &path_x,
        &path_y,
        last_row,
        last_col,
        &ampl_for_row,
        sequence.len(),
    );
    if !check {
        println!("Band length probably too short, maybe try with larger b and f");
    }
    drop(m);
    drop(x);
    drop(y);

    println!("{}", best_value);
    if seq_name.1 != 0 {
        gaf_output::gaf_of_gap_abpoa(
            &path,
            &path_x,
            &path_y,
            sequence,
            seq_name,
            &ampl_for_row,
            last_row,
            last_col,
            nodes_w_pred,
            file_path,
            amb_mode,
        );
    }

    best_value
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
    y: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
    o: i32,
) -> Option<(i32, usize, bool)> {
    let mut u_m = 0;
    let mut u_y = 0;
    let mut u_m_idx = 0;
    let mut u_y_idx = 0;
    let mut first = true;
    let left_i = ampl_for_row[i].0;
    for p in p_arr.iter() {
        let left_p = ampl_for_row[*p].0;
        if j + left_i >= ampl_for_row[*p].0 && j + left_i < ampl_for_row[*p].1 {
            let j_pos = if left_p < left_i {
                j + (left_i - left_p)
            } else {
                j - (left_p - left_i)
            };
            let current_u_m = m[*p][j_pos as usize] + o;
            let current_u_y = y[*p][j_pos as usize];
            if first {
                first = false;
                u_m = current_u_m;
                u_y = current_u_y;
                u_y_idx = *p;
                u_m_idx = *p
            }
            if current_u_m > u_m {
                u_m = current_u_m;
                u_m_idx = *p;
            }
            if current_u_y > u_y {
                u_y = current_u_y;
                u_y_idx = *p;
            }
        }
    }

    if first {
        None
    } else if u_y > u_m {
        Some((u_y, u_y_idx, false))
    } else {
        Some((u_m, u_m_idx, true))
    }
}

#[inline]
fn get_best_l(
    m: &[Vec<i32>],
    x: &[Vec<i32>],
    i: usize,
    j: usize,
    o: i32,
) -> Option<(i32, usize, bool)> {
    if j > 0 {
        let l_x = x[i][j - 1];
        let l_m = m[i][j - 1] + o;
        if l_x > l_m {
            Some((l_x, i, false))
        } else {
            Some((l_m, i, true))
        }
    } else {
        None
    }
}

#[inline]
fn set_ampl_for_row(
    i: usize,
    p_arr: &[usize],
    r_val: usize,
    best_scoring_pos: &[usize],
    seq_len: usize,
    bta: usize,
    simd_version: bool,
) -> (usize, usize) {
    let ms;
    let me;
    if i == 0 {
        ms = 0;
        me = 0;
    } else if p_arr.is_empty() {
        let pl = best_scoring_pos[i - 1];
        ms = pl + 1;
        me = pl + 1;
    } else {
        let mut pl = 0;
        let mut pr = 0;
        let mut first = true;
        for p in p_arr.iter() {
            let current_best = best_scoring_pos[*p];
            if first {
                pl = current_best;
                pr = current_best;
                first = false;
            }
            if current_best < pl {
                pl = current_best;
            }
            if current_best > pr {
                pr = current_best;
            }
        }
        ms = pl + 1;
        me = pr + 1;
    }
    let tmp_bs = cmp::min(ms as i32, (seq_len as i32 - r_val as i32) - bta as i32);
    let band_start = if tmp_bs < 0 {
        0
    } else {
        cmp::max(0, tmp_bs as usize)
    };
    let band_end = if seq_len > r_val {
        cmp::min(seq_len, cmp::max(me, seq_len - r_val) + bta)
    } else {
        cmp::min(seq_len, me + bta)
    };
    if simd_version {
        set_left_right_x64(band_start, band_end, seq_len)
    } else {
        (band_start, band_end)
    }
}

fn set_left_right_x64(left: usize, right: usize, seq_len: usize) -> (usize, usize) {
    let mut new_right = right;
    let mut new_left = left;
    while (new_right - new_left) % 64 != 0 {
        if (new_right - new_left) % 2 == 0 && new_right < seq_len {
            new_right += 1;
        } else if new_left > 0 {
            new_left -= 1;
        } else {
            break;
        }
    }
    (new_left, new_right)
}
fn set_r_values(
    nwp: &bit_vec::BitVec,
    pred_hash: &HashMap<usize, Vec<usize>>,
    lnz_len: usize,
) -> Vec<usize> {
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[lnz_len - 1] = 0;
    for p in pred_hash.get(&(lnz_len - 1)).unwrap() {
        r_values[*p] = 0;
    }
    for i in (1..lnz_len - 1).rev() {
        if r_values[i] == -1 || r_values[i] > r_values[i + 1] + 1 {
            r_values[i] = r_values[i + 1] + 1;
        }
        if nwp[i] {
            for p in pred_hash.get(&i).unwrap() {
                if r_values[*p] == -1 || r_values[*p] > r_values[i] + 1 {
                    r_values[*p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values.iter().map(|x| *x as usize).collect()
}
fn band_ampl_enough(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_x: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_y: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    start_row: usize,
    start_col: usize,
    ampl_for_row: &[(usize, usize)],
    sequence_len: usize,
) -> bool {
    let mut i = start_row;
    let mut j = start_col;
    while bf::dir_from_bitvec(&path[i][j]) != 'O' {
        let (left, right) = ampl_for_row[i];
        if i == 0 || j == 0 && left == 0 {
            return true;
        }
        if (j == left && left != 0) || (j == right - left - 1 && right != sequence_len) {
            return false;
        }
        let curr_bv = &path[i][j];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(curr_bv);
        match dir {
            'D' => {
                let left_p = ampl_for_row[pred].0;
                let j_pos = if left_p < left {
                    j + (left - left_p)
                } else {
                    j - (left_p - left)
                };
                j = j_pos - 1;
                i = pred;
            }
            'd' => {
                let left_p = ampl_for_row[pred].0;
                let j_pos = if left_p < left {
                    j + (left - left_p)
                } else {
                    j - (left_p - left)
                };
                j = j_pos - 1;
                i = pred;
            }
            'L' => {
                if bf::dir_from_bitvec(&path_x[i][j]) == 'X' {
                    while bf::dir_from_bitvec(&path_x[i][j]) == 'X' && j > 0 {
                        j -= 1;
                    }
                } else {
                    j -= 1
                }
            }
            'U' => {
                if bf::dir_from_bitvec(&path_y[i][j]) == 'Y' {
                    while bf::dir_from_bitvec(&path_y[i][j]) == 'Y' {
                        let left_row = ampl_for_row[i].0;
                        let p = bf::pred_from_bitvec(&path_y[i][j]);
                        let left_p = ampl_for_row[p].0;
                        let j_pos = if left_p < left_row {
                            j + (left_row - left_p)
                        } else {
                            j - (left_p - left_row)
                        };
                        j = j_pos;
                        i = p;
                    }
                } else {
                    let p = bf::pred_from_bitvec(&path[i][j]);
                    let left_p = ampl_for_row[p].0;
                    let j_pos = if left_p < left {
                        j + (left - left_p)
                    } else {
                        j - (left_p - left)
                    };
                    j = j_pos;
                    i = p;
                }
            }
            _ => {
                return false;
            }
        }
    }
    true
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

        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -4,
            -1,
            3,
            "prova.gfa",
            false,
        );

        assert_eq!(align, 4);
    }

    #[test]
    fn gap_correctly_considered() {
        let s = vec!['$', 'A', 'A', 'C', 'A', 'A', 'C'];

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
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);

        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -4,
            -1,
            3,
            "prova.gfa",
            false,
        );

        assert_eq!(align, 0);
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
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -4,
            -1,
            3,
            "prova.gfa",
            false,
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
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -4,
            -1,
            3,
            "prova.gfa",
            false,
        );

        assert_eq!(align, 5);
    }

    #[test]
    fn gap_poa_same_result_as_normal_if_o_0() {
        let s = vec!['$', 'A', 'A', 'C', 'A', 'A', 'C'];

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
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            0,
            -1,
            5,
            "prova.gfa",
            false,
        );

        assert_eq!(align, 4);
    }
    #[test]
    fn gap_open_only_once_if_penalty_high() {
        let s = vec!['$', 'A', 'A', 'A'];

        let lnz = vec!['$', 'A', 'C', 'A', 'C', 'A', 'F'];

        let mut nwp = BitVec::from_elem(7, false);
        nwp.set(1, true);
        nwp.set(6, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(6, vec![5]);
        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -100,
            -1,
            10,
            "prova.gfa",
            false,
        );

        assert_eq!(align, -101);
    }
    #[test]
    fn sequence_longer_than_graph() {
        let s = vec!['$', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'];

        let lnz = vec!['$', 'A', 'A', 'A', 'A', 'A', 'F'];
        let mut nwp = BitVec::from_elem(7, false);
        nwp.set(1, true);
        nwp.set(6, true);

        let mut pred_hash = HashMap::new();
        pred_hash.insert(1, vec![0]);
        pred_hash.insert(6, vec![5]);
        let graph = LnzGraph {
            lnz,
            nwp,
            pred_hash,
        };
        let mut score_matrix = HashMap::new();
        score_matrix.insert(('A', 'A'), 1);
        score_matrix.insert(('C', 'C'), 1);
        score_matrix.insert(('C', 'A'), -1);
        score_matrix.insert(('A', 'C'), -1);
        let align = super::exec(
            &s,
            ("test", 0),
            &graph,
            &score_matrix,
            -4,
            -1,
            7,
            "prova.gfa",
            false,
        );

        assert_eq!(align, -3);
    }
    #[test]
    fn ampl_for_row_multiple_of_64_if_simd() {
        let seq_len = 1000;
        let l1 = 100;
        let r1 = 200;
        let ampl1 = super::set_left_right_x64(l1, r1, seq_len);
        let delta1 = ampl1.1 - ampl1.0;

        let l2 = 5;
        let r2 = 15;
        let ampl2 = super::set_left_right_x64(l2, r2, seq_len);
        let delta2 = ampl1.1 - ampl1.0;

        let l3 = 989;
        let r3 = 990;
        let ampl3 = super::set_left_right_x64(l3, r3, seq_len);
        let delta3 = ampl1.1 - ampl1.0;

        let ampl4 = super::set_left_right_x64(5, 10, 15);

        assert_eq!(delta1 % 64, 0);

        assert_eq!(delta2 % 64, 0);
        assert_eq!(ampl2.0, 0);

        assert_eq!(delta3 % 64, 0);
        assert_eq!(ampl3.1, 1000);

        assert_eq!(ampl4.0, 0);
        assert_eq!(ampl4.1, 15);
    }
}
