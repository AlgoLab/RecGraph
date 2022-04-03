use std::{cmp::Ordering, collections::HashMap};

use crate::{basic_output, graph::LnzGraph};

pub fn exec(
    sequence: &[char],
    graph: &LnzGraph,
    scores_matrix: &HashMap<(char, char), i32>,
) -> i32 {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;

    let mut m = vec![vec![0; sequence.len()]; lnz.len()];
    let mut path = vec![vec![('x', 0); sequence.len()]; lnz.len()];
    let (mut best_row, mut best_col) = (0, 0);
    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() - 1 {
            match (i, j) {
                (0, _) | (_, 0) => path[i][j] = ('O', 0),
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
                        path[i][j] = ('O', 0);
                    } else {
                        let (best_val, mut dir) = get_best_d_u_l(d, u, l);
                        if dir == 'D' && lnz[i] != sequence[j] {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' | 'd' => (dir, d_idx),
                            'U' => (dir, u_idx),
                            _ => (dir, l_idx),
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

    println!("Best alignment: {}", m[best_row][best_col]);

    basic_output::write_align_local_poa(&path, sequence, lnz, best_row, best_col);
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

fn get_best_d_u_l(d: i32, u: i32, l: i32) -> (i32, char) {
    match d.cmp(&u) {
        Ordering::Less => match u.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (u, 'U'),
        },
        _ => match d.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (d, 'D'),
        },
    }
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
        let align_score = super::exec(&s, &graph_struct, &score_matrix);
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
        let align_score = super::exec(&s, &graph_struct, &score_matrix);
        assert_eq!(align_score, 2);
    }
}
