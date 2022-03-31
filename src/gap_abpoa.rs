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
    bta: usize,
) -> i32 {
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
            bta,
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
                let l_pred;
                match l_score_idx {
                    Some((l, idx)) => {
                        x[i][j] = l + e;
                        l_pred = idx;
                    }
                    _ => {
                        let best_p = if !nodes_w_pred[i] {
                            i - 1
                        } else {
                            *pred_hash.get(&i).unwrap().iter().min().unwrap()
                        };

                        x[i][j] = if best_p == 0 { o + e } else { x[best_p][j] + e };
                        l_pred = best_p;
                    }
                }
                // try set and get u from y (pred_left < j < pred_right else None)
                let u_score_idx = get_best_u(p_arr, &m, &y, &ampl_for_row, j, o);
                let u_pred;
                match u_score_idx {
                    Some((u, idx)) => {
                        y[i][j] = u + e;
                        u_pred = idx;
                    }
                    _ => {
                        y[i][j] = o + e * j as i32;
                        u_pred = 0;
                    }
                }
                // try get d from m (pred_left < j < pred_right else None)
                let d_score_idx = get_best_d(p_arr, &m, &ampl_for_row, j);
                match d_score_idx {
                    Some((mut d, d_idx)) => {
                        d += score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                        let l = x[i][j];
                        let u = y[i][j];
                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => match l.cmp(&u) {
                                Ordering::Less => {
                                    path[i][j] = ('U', u_pred);
                                    u
                                }
                                _ => {
                                    path[i][j] = ('L', l_pred);
                                    l
                                }
                            },
                            _ => match d.cmp(&u) {
                                Ordering::Less => {
                                    path[i][j] = ('U', u_pred);
                                    u
                                }
                                _ => {
                                    if lnz[i] == sequence[j] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx);
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
                                path[i][j] = ('U', u_pred);
                                u
                            }
                            _ => {
                                path[i][j] = ('L', l_pred);
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
        best_scoring_pos[i] = best_val_pos;
    }
    let mut last_row = m.len() - 2;
    let mut last_col = m[last_row].len() - 1;

    for p in pred_hash.get(&(m.len() - 1)).unwrap().iter() {
        let tmp_last_col = ampl_for_row[*p].1 - 1;
        if m[*p][tmp_last_col] > m[last_row][last_col] {
            last_row = *p;
            last_col = tmp_last_col;
        }
    }
    path.iter().for_each(|line| println!("{:?}", line));
    println!("{}", m[last_row][last_col]);
    m[last_row][last_col]
}
fn get_best_d(
    p_arr: &[usize],
    m: &[Vec<i32>],
    ampl_for_row: &[(usize, usize)],
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
        let l_m = m[i][j - 1] + o;
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
    bta: usize,
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

        let align = super::exec(&s, &graph, &score_matrix, -4, -1, 1);

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
        let align = super::exec(&s, &graph, &score_matrix, -4, -1, 10);

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
        let align = super::exec(&s, &graph, &score_matrix, -4, -1, 1);

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
        let align = super::exec(&s, &graph, &score_matrix, -4, -1, 1);

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
        let align = super::exec(&s, &graph, &score_matrix, 0, -1, 1);

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
        let align = super::exec(&s, &graph, &score_matrix, -100, -1, 10);

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
        let align = super::exec(&s, &graph, &score_matrix, -4, -1, 7);
        assert_eq!(align, -3);
    }
}
