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
)  {
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
        let (left, right) = set_ampl_for_row(i, p_arr, r_values[i], &best_scoring_pos, &ampl_for_row, sequence.len());
        ampl_for_row[i] = (left, right);
        let mut best_val_pos:usize = left;
        for j in left..right {
            if i == 0 && j == 0 {
                path[i][j] = ('O', 0);
            } else if i == 0 {
                // set y
                y[i][j] = o + e * i as i32;
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

                x[i][j] = if best_p == 0 {
                    o + e
                } else {
                   x[best_p][j] + e
                };
                
                // set m
                m[i][j] = x[i][j];
                path[i][j] = ('U', best_p);
            } else {
                // try set and get l from x (j > left, if j == left same as j = 0)
                // try set and get u from y (pred_left < j < pred_right else None)
                // try get d from m (pred_left < j < pred_right else None)
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
    println!("{:?}", ampl_for_row);
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
    let left = ampl_for_row[i].0;
    let mut first = true;
    for p in p_arr.iter() {
        if j + left >= ampl_for_row[*p].0 && j + left < ampl_for_row[*p].1 {
            let delta;
            let current_u_m;
            let current_u_y;

            if ampl_for_row[i].0 < ampl_for_row[*p].0 {
                delta = ampl_for_row[*p].0 - ampl_for_row[i].0;
                current_u_m = m[*p][j - delta] + o;
                current_u_y = y[*p][j - delta];
            } else {
                delta = ampl_for_row[i].0 - ampl_for_row[*p].0;
                current_u_m = m[*p][j + delta] + o;
                current_u_y = y[*p][j + delta];
            }

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
        println!("ERR best_u");
        None
    } else if u_y > u_m {
        Some((u_y, u_y_idx))
    } else {
        Some((u_m, u_m_idx))
    }
}

fn set_ampl_for_row(i: usize, p_arr: &[usize], r_val: usize, best_scoring_pos: &[usize],ampl_for_row: &Vec<(usize, usize)>, seq_len: usize) -> (usize, usize){
    let ms;
    let me;
    if i == 0 {
        ms = 0;
        me = 0;
    }
    else if p_arr.is_empty() {
        let pl = best_scoring_pos[i-1];
        ms = pl+1;
        me = pl+1;
        
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
            if  current_best + ampl_for_row[*p].0 < pl + left {
                pl = current_best;
                left = ampl_for_row[*p].0;
            }
            if current_best + ampl_for_row[*p].0 > pr + left {
                pr = current_best;
                left = ampl_for_row[*p].0;
            }
        }
        ms = pl+1;
        me = pr + 1;
    }
    let band_start;
    let band_end;
    if seq_len > r_val +10{
        band_start = cmp::max(0, cmp::min(ms, seq_len -r_val)-10);
        band_end = cmp::min(seq_len, cmp::max(me, seq_len-r_val)+10);
    } else {
        band_start = 0;
        band_end = cmp::min(seq_len, me+10);
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

    for i in 1..r_values.len(){
        match pred_hash.get(&i) {
            Some(arr) => {
                let best_p = arr.iter().max().unwrap();
                r_values[i] = r_values[*best_p] - 1;
            }
            _ => {
                r_values[i] = r_values[i-1] - 1;
            }
        }
    };
    r_values
}

