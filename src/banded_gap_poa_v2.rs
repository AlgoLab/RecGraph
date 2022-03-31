use crate::{basic_output, graph::LnzGraph};
use std::{
    cmp::{self, Ordering},
    collections::HashMap,
    vec,
};
pub fn exec(
    sequence: &[char],
    graph_struct: &LnzGraph,
    score_matrix: &HashMap<(char, char), i32>,
    o: i32,
    e: i32,
    bta: usize
) {
    let lnz = &graph_struct.lnz;
    let nodes_w_pred = &graph_struct.nwp;
    let pred_hash = &graph_struct.pred_hash;

    let mut m = vec![vec![0; sequence.len()]; lnz.len()]; // best alignment
    let mut x = vec![vec![0; sequence.len()]; lnz.len()]; //best alignment final gap in graph
    let mut y = vec![vec![0; sequence.len()]; lnz.len()]; // best alignment final gap in sequence

    let r_values = set_r_values(lnz.len(), pred_hash);
    let mut best_scoring_pos = vec![0; lnz.len()];

    let mut path = vec![vec![('X', 0); sequence.len()]; lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); lnz.len()];

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
            &ampl_for_row,
            sequence.len(),
            bta
        );
        ampl_for_row[i] = (left, right);
        let mut best_val_pos: usize = left;
        for j in left..right {
            if i == 0 && j == 0 {
                path[i][j] = ('O', 0);
            } else if i == 0 {
                // set y
                y[i][j] = o + e * j as i32;
                // set m
                m[i][j] = y[i][j];
                path[i][j] = ('L', i);
            } else if j == 0 {
                // set x

                let best_p = if !nodes_w_pred[i] {
                    i - 1
                } else {
                    *pred_hash.get(&i).unwrap().iter().min().unwrap()
                };

                x[i][j] = if best_p == 0 { o + e } else { x[best_p][j] + e };

                // set m
                m[i][j] = x[i][j];
                path[i][j] = ('U', best_p);
            } else {
                let mut p_arr = &vec![i - 1];
                if nodes_w_pred[i] {
                    p_arr = pred_hash.get(&i).unwrap();
                }

                // try set and get l from x (j > left, if j == left same as j = 0)
                let l_score_idx = get_best_l(&m, &x, &ampl_for_row, i, j, o);
                match l_score_idx {
                    Some((l, l_idx)) => {
                        x[i][j] = l + e;
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };

                        x[i][j] = if best_p == 0 { o + e } else { x[best_p][j] + e };
                    }
                }
                // try set and get u from y (pred_left < j < pred_right else None)
                let u_score_idx = get_best_u(p_arr, &m, &y, &ampl_for_row, i, j, o);
                match u_score_idx {
                    Some((u, u_idx)) => {
                        y[i][j] = u + e;
                    }
                    _ => {
                        y[i][j] = o + e * j as i32;
                    }
                }
                // try get d from m (pred_left < j < pred_right else None)
                let d_score_idx = get_best_d(p_arr, &m, &ampl_for_row, i, j);
                match d_score_idx {
                    Some((mut d, d_idx)) => {
                        d += score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                        let l = x[i][j];
                        let u = y[i][j];
                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => match l.cmp(&u) {
                                Ordering::Less => u,
                                _ => l,
                            },
                            _ => match d.cmp(&u) {
                                Ordering::Less => u,
                                _ => d,
                            },
                        }
                    }
                    _ => {
                        let l = x[i][j];
                        let u = x[i][j];
                        m[i][j] = match l.cmp(&u) {
                            Ordering::Less => u,
                            _ => l,
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
        best_scoring_pos[i] = best_val_pos;
    }
    //println!("{:?}", ampl_for_row);
    let last_row = m.len() - 2;
    let last_col = m[last_row].len() - 1;
    println!("{}", m[last_row][last_col]);
    println!("{:?}", m[0]);
    println!("{:?}", ampl_for_row[0]);

    ampl_for_row.iter().for_each(|l| {println!("{:?}", l)});
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
    for p in p_arr.iter() {
        if j > ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 + 1 {
            let curr_d = m[*p][j - 1];
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
fn get_best_u(
    p_arr: &[usize],
    m: &[Vec<i32>],
    y: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
    o: i32,
) -> Option<(i32, usize)> {
    let mut u_m = 0;
    let mut u_y = 0;
    let mut u_m_idx = 0;
    let mut u_y_idx = 0;
    let mut first = true;
    for p in p_arr.iter() {
        if j >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 {
            let current_u_m = m[*p][j] + o;
            let current_u_y = y[*p][j];

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
        Some((u_y, u_y_idx))
    } else {
        Some((u_m, u_m_idx))
    }
}

fn get_best_l(
    m: &[Vec<i32>],
    x: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
    o: i32,
) -> Option<(i32, usize)> {
    if j > ampl_for_row[i].0 {
        let l_x = x[i][j - 1];
        let l_m = m[i][j - 1];
        if l_x > l_m {
            Some((l_x, i))
        } else {
            Some((l_m, i))
        }
    } else {
        None
    }
}

fn set_ampl_for_row(
    i: usize,
    p_arr: &[usize],
    r_val: usize,
    best_scoring_pos: &[usize],
    ampl_for_row: &Vec<(usize, usize)>,
    seq_len: usize,
    bta: usize
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
        let mut left = 0;
        let mut first = true;
        for p in p_arr.iter() {
            let current_best = best_scoring_pos[*p];
            if first {
                pl = current_best;
                pr = current_best;
                left = ampl_for_row[*p].0;
                first = false;
            }
            if current_best + ampl_for_row[*p].0 < pl + left {
                pl = current_best;
                left = ampl_for_row[*p].0;
            }
            if current_best + ampl_for_row[*p].0 > pr + left {
                pr = current_best;
                left = ampl_for_row[*p].0;
            }
        }
        ms = pl + 1;
        me = pr + 1;
    }
    let band_start;
    let band_end;
    let tmp_bs = cmp::min(ms as i32, (seq_len as i32 - r_val as i32) - bta as i32);
    if tmp_bs < 0 {
        band_start = 0;
    } else {
        band_start = cmp::max(0, tmp_bs as usize);
    }
    if seq_len > r_val {
        band_end = cmp::min(seq_len, cmp::max(me, seq_len - r_val) + bta);
    } else {
        band_end = cmp::min(seq_len, me + bta);
    }

    (band_start, band_end)
}

fn set_r_values(lnz_len: usize, pred_hash: &HashMap<usize, Vec<usize>>) -> Vec<usize> {
    let mut r_values = vec![0; lnz_len];
    let mut i = lnz_len - 1;
    let mut count = 0;
    while i > 0 {
        match pred_hash.get(&i) {
            Some(arr) => {
                i = *arr.iter().max().unwrap();
            }
            _ => {
                i -= 1;
            }
        }
        count += 1;
    }
    r_values[0] = count;

    for i in 1..r_values.len() {
        match pred_hash.get(&i) {
            Some(arr) => {
                let best_p = arr.iter().max().unwrap();
                r_values[i] = r_values[*best_p] - 1;
            }
            _ => {
                r_values[i] = r_values[i - 1] - 1;
            }
        }
    }
    r_values
}
