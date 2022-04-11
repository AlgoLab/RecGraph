use crate::bitfield_path::{self as bf, pred_from_bitvec};
use crate::{basic_output, graph::LnzGraph};
use bitvec::prelude::*;
use std::{
    cmp::{self, Ordering},
    collections::HashMap,
    vec,
};
//TODO: remove last row, only best value needed
pub fn exec(
    sequence: &[char],
    graph_struct: &LnzGraph,
    score_matrix: &HashMap<(char, char), i32>,
    bta: usize,
) -> i32 {
    let lnz = &graph_struct.lnz;
    let nodes_w_pred = &graph_struct.nwp;
    let pred_hash = &graph_struct.pred_hash;

    let r_values = set_r_values(nodes_w_pred, pred_hash, lnz.len());
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
        let bv = bitvec![u16, Msb0; 0; 32];
        path[i] = vec![bv; right - left];
        for j in 0..right - left {
            match (i, j) {
                (0, 0) => {
                    m[i][j] = 0;
                    path[i][j] = bf::set_path_cell(0, 'O');
                }
                (0, _) => {
                    //only left
                    m[i][j] = m[i][j - 1] + score_matrix.get(&('-', sequence[j + left])).unwrap();
                    path[i][j] = bf::set_path_cell(i, 'L');
                }
                (_, 0)
                    if left == 0 || left_equal_for_every_p(pred_hash.get(&i), &ampl_for_row, i) =>
                {
                    // only upper
                    if !nodes_w_pred[i] {
                        m[i][j] = m[i - 1][j] + score_matrix.get(&(lnz[i], '-')).unwrap();
                        path[i][j] = bf::set_path_cell(i - 1, 'U');
                    } else {
                        let p = pred_hash.get(&i).unwrap().iter().min().unwrap();
                        m[i][j] = m[*p][j] + score_matrix.get(&(lnz[i], '-')).unwrap();
                        path[i][j] = bf::set_path_cell(*p, 'U');
                    }
                }
                (_, 0) if left > 0 => {
                    //only u or d
                    if !nodes_w_pred[i] {
                        let delta;
                        let u;
                        let d;
                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            u = m[i - 1][j - delta] + score_matrix.get(&('-', lnz[i])).unwrap();
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            u = m[i - 1][j + delta] + score_matrix.get(&('-', lnz[i])).unwrap();
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        }
                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = bf::set_path_cell(i - 1, 'U');
                                u
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = bf::set_path_cell(i - 1, 'D');
                                } else {
                                    path[i][j] = bf::set_path_cell(i - 1, 'd');
                                }
                                d
                            }
                        }
                    } else if let (Some((mut d, d_idx)), Some((mut u, u_idx))) = (
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                        get_best_u(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                    ) {
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        u += score_matrix.get(&(lnz[i], '-')).unwrap();
                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = bf::set_path_cell(u_idx, 'U');
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
                        }
                    }
                }
                _ if j == right - left - 1 => {
                    //only d or l
                    let l = m[i][j - 1]
                        + score_matrix
                            .get(&(sequence[j + ampl_for_row[i].0], '-'))
                            .unwrap();
                    if !nodes_w_pred[i] {
                        let d;
                        let delta;

                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        }

                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => {
                                path[i][j] = bf::set_path_cell(i, 'L');
                                l
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = bf::set_path_cell(i - 1, 'D');
                                } else {
                                    path[i][j] = bf::set_path_cell(i - 1, 'd');
                                }
                                d
                            }
                        }
                    } else if let Some((mut d, d_idx)) =
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j)
                    {
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();

                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => {
                                path[i][j] = bf::set_path_cell(i, 'L');
                                l
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = bf::set_path_cell(d_idx, 'D');
                                } else {
                                    path[i][j] = bf::set_path_cell(d_idx, 'd');
                                }
                                d
                            }
                        }
                    }
                }
                _ => {
                    //every value ok
                    let l = m[i][j - 1]
                        + score_matrix
                            .get(&(sequence[j + ampl_for_row[i].0], '-'))
                            .unwrap();
                    if !nodes_w_pred[i] {
                        let delta;
                        let u;
                        let d;
                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            u = m[i - 1][j - delta] + score_matrix.get(&('-', lnz[i])).unwrap();
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            u = m[i - 1][j + delta] + score_matrix.get(&('-', lnz[i])).unwrap();
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();
                        }
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != lnz[i] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => bf::set_path_cell(i - 1, 'D'),
                            'd' => bf::set_path_cell(i - 1, 'd'),
                            'U' => bf::set_path_cell(i - 1, 'U'),
                            _ => bf::set_path_cell(i, 'L'),
                        };
                    } else if let (Some((mut d, d_idx)), Some((mut u, u_idx))) = (
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                        get_best_u(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                    ) {
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        u += score_matrix.get(&(lnz[i], '-')).unwrap();
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != lnz[i] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => bf::set_path_cell(d_idx, 'D'),
                            'd' => bf::set_path_cell(d_idx, 'd'),
                            'U' => bf::set_path_cell(u_idx, 'U'),
                            _ => bf::set_path_cell(i, 'L'),
                        };
                    }
                }
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
    let best_value = m[last_row][last_col];
    let check = ampl_is_enough(&path, &ampl_for_row, sequence.len(), last_row, last_col);
    if !check {
        println!("Band length probably too short, maybe try with larger b and f");
    }
    println!("Alignment mk {:?}", best_value);
    basic_output::write_align_banded_poa(&path, sequence, lnz, &ampl_for_row, last_row, last_col);

    m[last_row][last_col]
}

fn ampl_is_enough(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    ampl_for_row: &[(usize, usize)],
    seq_len: usize,
    last_row: usize,
    last_col: usize,
) -> bool {
    let mut row = last_row;
    let mut col = last_col;

    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        //reached end of path, no need to continue
        if ampl_for_row[row].0 == 0 {
            return true;
        }

        let p_left = ampl_for_row[bf::pred_from_bitvec(&path[row][col])].0;
        let j_pos = if ampl_for_row[row].0 < p_left {
            let delta = p_left - ampl_for_row[row].0;
            col - delta
        } else {
            let delta = ampl_for_row[row].0 - p_left;
            col + delta
        };
        if col == 0
            || (col == ampl_for_row[row].1 - ampl_for_row[row].0 - 1
                && ampl_for_row[row].1 != seq_len - 1)
        {
            // != from path_len because couldn't go larger
            if bf::dir_from_bitvec(&path[row][col]) == 'D' {
                row = pred_from_bitvec(&path[row][col]);
                col = j_pos - 1;
                // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            let curr_bv = &path[row][col];
            let pred = bf::pred_from_bitvec(curr_bv);
            let dir = bf::dir_from_bitvec(&curr_bv);
            match dir {
                'D' | 'd' => {
                    row = pred;
                    col = j_pos - 1;
                }
                'U' => {
                    row = pred;
                    col = j_pos;
                }
                'L' => {
                    col -= 1;
                }
                _ => {
                    return false;
                }
            }
        }
    }
    true
}

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
    let left = ampl_for_row[i].0;
    for p in p_arr.iter() {
        if j + left >= ampl_for_row[*p].0 + 1 && j + left < ampl_for_row[*p].1 + 1 {
            let delta;
            let current_d;
            if ampl_for_row[i].0 < ampl_for_row[*p].0 {
                delta = ampl_for_row[*p].0 - ampl_for_row[i].0;
                current_d = m[*p][j - delta - 1];
            } else {
                delta = ampl_for_row[i].0 - ampl_for_row[*p].0;
                current_d = m[*p][j + delta - 1];
            }
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
    }

    if first {
        // j is too far to be aligned with current predecessors
        println!("ERR best_d: {} {}", i, j); // if no pred found should panic
        None
    } else {
        Some((d, d_idx))
    }
}

fn get_best_u(
    p_arr: &[usize],
    m: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut u = 0;
    let mut u_idx = 0;
    let left = ampl_for_row[i].0;
    let mut first = true;
    for p in p_arr.iter() {
        if j + left >= ampl_for_row[*p].0 && j + left < ampl_for_row[*p].1 {
            let delta;
            let current_u;
            if ampl_for_row[i].0 < ampl_for_row[*p].0 {
                delta = ampl_for_row[*p].0 - ampl_for_row[i].0;
                current_u = m[*p][j - delta];
            } else {
                delta = ampl_for_row[i].0 - ampl_for_row[*p].0;
                current_u = m[*p][j + delta];
            }

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
        println!("ERR best_u");
        None
    } else {
        Some((u, u_idx))
    }
}

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
fn left_equal_for_every_p(
    p_arr: Option<&Vec<usize>>,
    ampl_for_row: &[(usize, usize)],
    i: usize,
) -> bool {
    if let Some(arr) = p_arr {
        let mut check = true;
        for p in arr.iter() {
            if ampl_for_row[*p].0 != ampl_for_row[i].0 {
                check = false
            }
        }
        check
    } else {
        ampl_for_row[i - 1].0 == ampl_for_row[i].0
    }
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
        let align = super::exec(&s, &graph, &score_matrix, 4);

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
        let align = super::exec(&s, &graph, &score_matrix, 4);

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
        let align = super::exec(&s, &graph, &score_matrix, 4);

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
        let align = super::exec(&s, &graph, &score_matrix, 4);

        assert_eq!(align, 5);
    }
}
