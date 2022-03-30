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
    ampl: usize,
    o: i32,
    e: i32,
) -> i32 {
    let lnz = &graph_struct.lnz;
    let nodes_w_pred = &graph_struct.nwp;
    let pred_hash = &graph_struct.pred_hash;
    let mut m = vec![vec![]; lnz.len()]; // best alignment
    let mut x = vec![vec![]; lnz.len()]; //best alignment final gap in graph
    let mut y = vec![vec![]; lnz.len()]; // best alignment final gap in sequence
    let r_values = set_r_values(lnz.len(), pred_hash);
    let best_scoring_pos = vec![0; lnz.len()];

    let mut path = vec![vec![]; lnz.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); lnz.len()];
    for i in 0..lnz.len() {
        let mut p_arr = &vec![];
        if nodes_w_pred[i] {
            p_arr = pred_hash.get(&i).unwrap()
        }
        let (left, right) = set_left_right(ampl, i, p_arr, sequence.len(), &mut ampl_for_row);
        m[i] = vec![0; right - left];
        x[i] = vec![0; right - left];
        y[i] = vec![0; right - left];
        path[i] = vec![('X', 0); right - left];
    }

    for i in 0..lnz.len() - 1 {
        let mut p_arr = &vec![];
        if nodes_w_pred[i] {
            p_arr = pred_hash.get(&i).unwrap()
        }
        set_ampl_for_row_v2(i, p_arr, r_values[i], &best_scoring_pos);
        let (left, right) = ampl_for_row[i];
        for j in 0..right - left {
            match (i, j) {
                (0, 0) => {
                    m[i][j] = 0;
                    x[i][j] = 0;
                    y[i][j] = 0;

                    path[i][j] = ('O', 0);
                }
                (0, _) => {
                    //only left
                    //set y
                    y[i][j] = o + e * (j + ampl_for_row[0].0) as i32;
                    // set m
                    m[i][j] = y[i][j];
                    path[i][j] = ('L', i);
                }
                (_, 0) if left == 0 => {
                    // only upper
                    if !nodes_w_pred[i] {
                        //set x
                        x[i][j] = o + e * i as i32;
                        // set m
                        m[i][j] = x[i][j];
                        path[i][j] = ('U', i - 1);
                    } else {
                        let p = pred_hash.get(&i).unwrap().iter().min().unwrap();
                        // set x
                        x[i][j] = o + e * (*p + 1) as i32;
                        // set m
                        m[i][j] = x[i][j];
                        path[i][j] = ('U', *p);
                    }
                }
                (_, 0) if left_equal_for_every_p(pred_hash.get(&i), &ampl_for_row, i) => {
                    //only upper value
                    if !nodes_w_pred[i] {
                        // set y
                        let j_pos = if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            let delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            j - delta
                        } else {
                            let delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            j + delta
                        };
                        y[i][j] = cmp::max(m[i - 1][j_pos] + o + e, y[i - 1][j_pos] + e);
                        // set m
                        m[i][j] = y[i][j];
                        path[i][j] = ('U', i - 1);
                    } else {
                        let (u, u_idx) =
                            get_best_u(pred_hash.get(&i).unwrap(), &m, &y, &ampl_for_row, i, j, o)
                                .unwrap();
                        // set y
                        y[i][j] = u + e;
                        //set m
                        x[i][j] = y[i][j];
                        path[i][j] = ('U', u_idx);
                    }
                }
                (_, 0) if left > 0 => {
                    //only u or d
                    if !nodes_w_pred[i] {
                        let j_pos = if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            let delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            j - delta
                        } else {
                            let delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            j + delta
                        };

                        // set y
                        y[i][j] = cmp::max(m[i - 1][j_pos] + o + e, y[i - 1][j_pos] + e);

                        // set m
                        let u = y[i][j];
                        let d = m[i - 1][j_pos - 1]
                            + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();

                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = ('U', i - 1);
                                u
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1);
                                }
                                d
                            }
                        }
                    } else if let (Some((mut d, d_idx)), Some((mut u, u_idx))) = (
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                        get_best_u(pred_hash.get(&i).unwrap(), &m, &y, &ampl_for_row, i, j, o),
                    ) {
                        // set y
                        y[i][j] = u + e;

                        //set m
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        u += e;
                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = ('U', u_idx);
                                u
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = ('D', d_idx);
                                } else {
                                    path[i][j] = ('d', d_idx);
                                }
                                d
                            }
                        }
                    }
                }
                _ if j == right - left - 1 => {
                    //only d or l
                    // set x
                    x[i][j] = cmp::max(m[i][j - 1] + o + e, x[i][j - 1] + e);
                    let l = x[i][j];

                    if !nodes_w_pred[i] {
                        let d;
                        let delta;
                        // set m
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
                                path[i][j] = ('L', i);
                                l
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1);
                                }
                                d
                            }
                        }
                    } else if let Some((mut d, d_idx)) =
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j)
                    {
                        // set m
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();

                        m[i][j] = match d.cmp(&l) {
                            Ordering::Less => {
                                path[i][j] = ('L', i);
                                l
                            }
                            _ => {
                                if lnz[i] == sequence[j + left] {
                                    path[i][j] = ('D', d_idx);
                                } else {
                                    path[i][j] = ('d', d_idx);
                                }
                                d
                            }
                        }
                    }
                }
                _ => {
                    //every value ok
                    // set x
                    x[i][j] = cmp::max(m[i][j - 1] + o + e, x[i][j - 1] + e);
                    let l = x[i][j];

                    if !nodes_w_pred[i] {
                        let j_pos = if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            let delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            j - delta
                        } else {
                            let delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            j + delta
                        };

                        // set y
                        y[i][j] = cmp::max(m[i - 1][j_pos] + o + e, y[i - 1][j_pos] + e);

                        // set m
                        let u = y[i][j];
                        let d = m[i - 1][j_pos - 1]
                            + score_matrix.get(&(sequence[j + left], lnz[i])).unwrap();

                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != lnz[i] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => ('D', i - 1),
                            'd' => ('d', i - 1),
                            'U' => ('U', i - 1),
                            _ => ('L', i),
                        };
                    } else if let (Some((mut d, d_idx)), Some((mut u, u_idx))) = (
                        get_best_d(pred_hash.get(&i).unwrap(), &m, &ampl_for_row, i, j),
                        get_best_u(pred_hash.get(&i).unwrap(), &m, &y, &ampl_for_row, i, j, o),
                    ) {
                        // set y
                        y[i][j] = u + e;

                        //set m
                        d += score_matrix.get(&(lnz[i], sequence[j + left])).unwrap();
                        u += e;
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != lnz[i] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => ('D', d_idx),
                            'd' => ('d', d_idx),
                            'U' => ('U', u_idx),
                            _ => ('L', i),
                        };
                    }
                }
            }
        }
    }
    let last_col_f_node = ampl_for_row[m.len() - 1].1 - ampl_for_row[m.len() - 1].0 - 1;
    best_last_node(
        pred_hash.get(&(lnz.len() - 1)).unwrap(),
        &mut m,
        &mut path,
        last_col_f_node,
        &ampl_for_row,
    );
    let last_row = path[m.len() - 1][last_col_f_node].1;
    let last_col = ampl_for_row[last_row].1 - ampl_for_row[last_row].0 - 1;
    match ampl_is_enough(&path, &ampl_for_row, sequence.len()) {
        true => {
            println!("Alignment mk {:?}", m[m.len() - 1][last_col_f_node]);

            basic_output::write_align_banded_poa(
                &path,
                sequence,
                lnz,
                &ampl_for_row,
                last_row,
                last_col,
            );

            m[last_row][last_col]
        }
        false => exec(sequence, graph_struct, score_matrix, ampl * 2, o, e),
    }
}

fn best_last_node(
    p_arr: &[usize],
    m: &mut [Vec<i32>],
    path: &mut [Vec<(char, usize)>],
    j: usize,
    ampl_for_row: &[(usize, usize)],
) {
    let mut best_align = 0;
    let mut best_idx = 0;
    let mut first = true;
    let m_len = m.len();
    let last_row = ampl_for_row.len() - 1;

    for p in p_arr.iter() {
        let last_col = ampl_for_row[*p].1 - ampl_for_row[*p].0 - 1;
        if first {
            best_align = m[*p][last_col];
            best_idx = *p;
            first = false;
        }
        if m[*p][last_col] > best_align {
            best_align = m[*p][last_col];
            best_idx = *p;
        }
    }
    m[m_len - 1][j] = best_align;
    path[m_len - 1][j] = ('F', best_idx);
}

fn set_left_right(
    ampl: usize,
    i: usize,
    pred_array: &[usize],
    sequence_len: usize,
    ampl_for_row: &mut [(usize, usize)],
) -> (usize, usize) {
    let mut left = 0;
    let mut right = sequence_len;
    if ampl < sequence_len {
        if i == 0 {
            right = cmp::min(ampl / 2, sequence_len);
        } else if pred_array.is_empty() {
            if i >= ampl / 2 {
                left = ampl_for_row[i - 1].0 + 1
            }
            right = cmp::min(ampl_for_row[i - 1].1 + 1, sequence_len);
            if right <= left + ampl / 2 {
                (left, right) = (right - ampl / 2, right);
            }
        } else {
            let mut first = true;
            for p in pred_array.iter() {
                let (current_l, current_r);
                if p + 1 < ampl / 2 {
                    current_l = 0;
                } else {
                    current_l = ampl_for_row[*p].0 + 1
                }
                current_r = cmp::min(ampl_for_row[*p].1 + 1, sequence_len);
                if first {
                    first = false;
                    (left, right) = (current_l, current_r)
                }
                if current_l < left {
                    left = current_l;
                }
                if current_r > right {
                    right = current_r;
                }
            }
            if right <= left + ampl / 2 {
                (left, right) = (right - ampl / 2, right);
            }
        }
    }
    ampl_for_row[i] = (left, right);
    (left, right)
}

fn ampl_is_enough(
    path: &[Vec<(char, usize)>],
    ampl_for_row: &[(usize, usize)],
    seq_len: usize,
) -> bool {
    let last_row = path.len() - 1;
    let last_col = ampl_for_row[last_row].1 - ampl_for_row[last_row].0 - 1;
    let mut row = path[last_row][last_col].1;
    let mut col = path[row].len() - 1;

    while path[row][col].0 != 'O' {
        //reached end of path, no need to continue

        if ampl_for_row[row].0 == 0 {
            return true;
        }
        let p_left = ampl_for_row[path[row][col].1].0;
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
            if path[row][col].0 == 'D' {
                row = path[row][col].1 as usize;
                col = j_pos - 1;
                // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            match path[row][col].0 {
                'D' | 'd' => {
                    row = path[row][col].1 as usize;
                    col = j_pos - 1;
                }
                'U' => {
                    row = path[row][col].1 as usize;
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
    ampl_for_row: &Vec<(usize, usize)>,
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
        if ampl_for_row[i - 1].0 != ampl_for_row[i].0 {
            false
        } else {
            true
        }
    }
}

fn set_ampl_for_row_v2(i: usize, p_arr: &[usize], r_val: usize, best_scoring_pos: &[usize],) {
    let lnz_len = best_scoring_pos.len();
    let ms;
    let me;

    if p_arr.is_empty() {
        let pl = best_scoring_pos[i];
        ms = pl+1;
        me = pl+1;
        
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
            if  current_best < pl {
                pl = current_best;
            }
            if current_best > pr {
                pr = current_best;
            }
        }
        ms = pl+1;
        me = pr + 1;
    }
    let band_start = cmp::max(0, cmp::min(ms, lnz_len -r_val));
    let band_end = cmp::min(lnz_len, cmp::max(me, lnz_len-r_val));
    println!("{} {}",band_start, band_end);
}

fn set_r_values(lnz_len: usize, pred_hash: &HashMap<usize, Vec<usize>>) -> Vec<usize> {
    let mut r_values = vec![0; lnz_len];
    let mut i = lnz_len - 1;
    let mut count = 0;
    while i > 0 {
        match pred_hash.get(&i) {
            Some(arr) => {
                i = arr[0];
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
                r_values[i] = r_values[arr[0]] - 1;
            }
            _ => {
                r_values[i] = r_values[i-1] - 1;
            }
        }
    };
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

        let align = super::exec(&s, &graph, &score_matrix, 4, -4, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 4, -4, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 4, -4, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 4, -4, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 4, 0, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 4, -100, -1);

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
        let align = super::exec(&s, &graph, &score_matrix, 10, 0, -1);
        assert_eq!(align, 1);
    }
}
