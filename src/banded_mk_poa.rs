use crate::basic_output;
use std::{
    cmp::{self, Ordering},
    collections::HashMap,
};
pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
) {
    let mut m = vec![vec![]; graph.len()];
    let mut path = vec![vec![]; graph.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); graph.len()];
    for i in 0..graph.len() {
        let (left, right) = set_left_right(ampl, i, graph, sequence, &mut ampl_for_row);
        m[i] = vec![0; right - left];
        path[i] = vec![('X', 0); right - left];
    }

    for i in 0..graph.len() - 1 {
        let (left, right) = ampl_for_row[i];
        for j in 0..right - left {
            match (i, j) {
                (0, 0) => {
                    m[i][j] = 0;
                    path[i][j] = ('O', 0);
                }
                (0, _) => {
                    //only left
                    m[i][j] = m[i][j - 1] + score_matrix.get(&('-', sequence[j + left])).unwrap();
                    path[i][j] = ('L', j - 1);
                }
                (_, 0) if left == 0 || ampl_for_row[i].0 == ampl_for_row[i - 1].0 => {
                    // only upper
                    if graph[i].1.is_empty() {
                        m[i][j] = m[i - 1][j] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        path[i][j] = ('U', i - 1);
                    } else {
                        let p = graph[i].1.iter().min().unwrap();
                        m[i][j] = m[*p][j] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        path[i][j] = ('U', *p);
                    }
                }
                (_, 0) if left > 0 => {
                    //only u or d
                    if graph[i].1.is_empty() {
                        let delta;
                        let u;
                        let d;
                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            u = m[i - 1][j - delta] + score_matrix.get(&('-', graph[i].0)).unwrap();
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            u = m[i - 1][j + delta] + score_matrix.get(&('-', graph[i].0)).unwrap();
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        }
                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = ('U', i - 1);
                                u
                            }
                            _ => {
                                if graph[i].0 == sequence[j + left] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1);
                                }
                                d
                            }
                        }
                    } else if let (Some((d, d_idx)), Some((u, u_idx))) = (
                        get_best_d(graph, sequence, &m, score_matrix, &ampl_for_row, i, j),
                        get_best_u(graph, &m, score_matrix, &ampl_for_row, i, j),
                    ) {
                        m[i][j] = match d.cmp(&u) {
                            Ordering::Less => {
                                path[i][j] = ('U', u_idx);
                                u
                            }
                            _ => {
                                if graph[i].0 == sequence[j + left] {
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
                    //only d or l, u only if band ends predecessor and current node equals
                    let l = m[i][j - 1]
                        + score_matrix
                            .get(&(sequence[j + ampl_for_row[i].0], '-'))
                            .unwrap();
                    if graph[i].1.is_empty() {
                        let d;
                        let delta;

                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        }
                        let u = if cannot_look_up(i - 1, i, &ampl_for_row) {
                            d - 1
                        } else {
                            // can have u value even if last char of band
                            m[i - 1][j + delta] + score_matrix.get(&('-', graph[i].0)).unwrap()
                        };
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != graph[i].0 {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => ('D', i - 1),
                            'd' => ('d', i - 1),
                            'U' => ('U', i - 1),
                            _ => ('L', j - 1),
                        };
                    } else if let Some((d, d_idx)) =
                        get_best_d(graph, sequence, &m, score_matrix, &ampl_for_row, i, j)
                    {
                        if let Some((u, u_idx)) =
                            get_best_u_special(graph, &m, score_matrix, &ampl_for_row, i, j)
                        // same as case before
                        {
                            let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                            if dir == 'D' && sequence[j + left] != graph[i].0 {
                                dir = 'd'
                            }
                            m[i][j] = best_val;
                            path[i][j] = match dir {
                                'D' => ('D', d_idx),
                                'd' => ('d', d_idx),
                                'U' => ('U', u_idx),
                                _ => ('L', j - 1),
                            };
                        } else {
                            m[i][j] = match d.cmp(&l) {
                                Ordering::Less => {
                                    path[i][j] = ('L', j - 1);
                                    l
                                }
                                _ => {
                                    if graph[i].0 == sequence[j + left] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx);
                                    }
                                    d
                                }
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
                    if graph[i].1.is_empty() {
                        let delta;
                        let u;
                        let d;
                        if ampl_for_row[i].0 < ampl_for_row[i - 1].0 {
                            delta = ampl_for_row[i - 1].0 - ampl_for_row[i].0;
                            u = m[i - 1][j - delta] + score_matrix.get(&('-', graph[i].0)).unwrap();
                            d = m[i - 1][j - delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        } else {
                            delta = ampl_for_row[i].0 - ampl_for_row[i - 1].0;
                            u = m[i - 1][j + delta] + score_matrix.get(&('-', graph[i].0)).unwrap();
                            d = m[i - 1][j + delta - 1]
                                + score_matrix.get(&(sequence[j + left], graph[i].0)).unwrap();
                        }
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != graph[i].0 {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => ('D', i - 1),
                            'd' => ('d', i - 1),
                            'U' => ('U', i - 1),
                            _ => ('L', j - 1),
                        };
                    } else if let (Some((d, d_idx)), Some((u, u_idx))) = (
                        get_best_d(graph, sequence, &m, score_matrix, &ampl_for_row, i, j),
                        get_best_u(graph, &m, score_matrix, &ampl_for_row, i, j),
                    ) {
                        let (best_val, mut dir) = get_max_d_u_l(d, u, l);
                        if dir == 'D' && sequence[j + left] != graph[i].0 {
                            dir = 'd'
                        }
                        m[i][j] = best_val;
                        path[i][j] = match dir {
                            'D' => ('D', d_idx),
                            'd' => ('d', d_idx),
                            'U' => ('U', u_idx),
                            _ => ('L', j - 1),
                        };
                    }
                }
            }
        }
    }

    let last_col = ampl_for_row[m.len() - 1].1 - ampl_for_row[m.len() - 1].0 - 1;
    best_last_node(graph, &mut m, &mut path, last_col, &ampl_for_row);
    let last_row = path[m.len() - 1][last_col].1;
    match ampl_is_enough(&path, &ampl_for_row, sequence.len()) {
        true => {
            println!("Alignment mk {:?}", m[m.len() - 1][last_col]);
            basic_output::write_align_banded_poa(
                &path,
                sequence,
                graph,
                &ampl_for_row,
                last_row,
                last_col,
            );

            //m[graph.len() - 1][sequence.len() - 1]
        }
        false => exec(sequence, graph, score_matrix, ampl * 2),
    }
}

fn best_last_node(
    graph: &[(char, Vec<usize>)],
    m: &mut [Vec<i32>],
    path: &mut [Vec<(char, usize)>],
    j: usize,
    ampl_for_row: &[(usize, usize)],
) {
    let mut best_align = 0;
    let mut best_idx = 0;
    let mut first = true;
    let m_len = m.len();

    for p in graph[graph.len() - 1].1.iter() {
        let delta;
        let j_pos;
        if ampl_for_row[graph.len() - 1].0 >= ampl_for_row[*p].0 {
            delta = ampl_for_row[graph.len() - 1].0 - ampl_for_row[*p].0;
            j_pos = j + delta;
        } else {
            delta = ampl_for_row[*p].0 - ampl_for_row[graph.len() - 1].0;
            j_pos = j - delta;
        }
        if first {
            best_align = m[*p][j_pos];
            best_idx = *p;
            first = false;
        }
        if m[*p][j_pos] > best_align {
            best_align = m[*p][j_pos];
            best_idx = *p;
        }
    }
    m[m_len - 1][j] = best_align;
    path[m_len - 1][j] = ('F', best_idx);
}

fn set_left_right(
    ampl: usize,
    i: usize,
    graph: &[(char, Vec<usize>)],
    sequence: &[char],
    ampl_for_row: &mut [(usize, usize)],
) -> (usize, usize) {
    let mut left;
    let mut right;
    if ampl >= sequence.len() {
        left = 0;
        right = sequence.len();
    } else {
        left = 0;
        right = sequence.len();
        if i == 0 {
            right = cmp::min(ampl / 2, sequence.len());
        } else if graph[i].1.is_empty() {
            if i < ampl / 2 {
                left = 0;
            } else {
                left = ampl_for_row[i - 1].0 + 1
            }
            right = cmp::min(ampl_for_row[i - 1].1 + 1, sequence.len());
            if right <= left + ampl / 2 {
                (left, right) = (right - ampl / 2, right);
            }
        } else {
            let mut first = true;
            for p in graph[i].1.iter() {
                let (current_l, current_r);
                if p + 1 < ampl / 2 {
                    current_l = 0;
                } else {
                    current_l = ampl_for_row[*p].0 + 1
                }
                current_r = cmp::min(ampl_for_row[*p].1 + 1, sequence.len());
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
    graph: &[(char, Vec<usize>)],
    sequence: &[char],
    m: &[Vec<i32>],
    scores_matrix: &HashMap<(char, char), i32>,
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut d = 0;
    let mut d_idx = 0;
    let mut first = true;
    let left = ampl_for_row[i].0;
    for p in graph[i].1.iter() {
        if j + left >= ampl_for_row[*p].0 + 1 && j + left < ampl_for_row[*p].1 + 1 {
            let delta;
            let current_d;
            if ampl_for_row[i].0 < ampl_for_row[*p].0 {
                delta = ampl_for_row[*p].0 - ampl_for_row[i].0;
                current_d = m[*p][j - delta - 1]
                    + scores_matrix
                        .get(&(graph[i].0, sequence[j + left]))
                        .unwrap();
            } else {
                delta = ampl_for_row[i].0 - ampl_for_row[*p].0;
                current_d = m[*p][j + delta - 1]
                    + scores_matrix
                        .get(&(graph[i].0, sequence[j + left]))
                        .unwrap();
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
        None
    } else {
        Some((d, d_idx))
    }
}

fn get_best_u(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    scores_matrix: &HashMap<(char, char), i32>,
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut u = 0;
    let mut u_idx = 0;
    let left = ampl_for_row[i].0;
    let mut first = true;
    for p in graph[i].1.iter() {
        if j + left >= ampl_for_row[*p].0 && j + left < ampl_for_row[*p].1 {
            let delta;
            let current_u;
            if ampl_for_row[i].0 < ampl_for_row[*p].0 {
                delta = ampl_for_row[*p].0 - ampl_for_row[i].0;
                current_u = m[*p][j - delta] + scores_matrix.get(&(graph[i].0, '-')).unwrap();
            } else {
                delta = ampl_for_row[i].0 - ampl_for_row[*p].0;
                current_u = m[*p][j + delta] + scores_matrix.get(&(graph[i].0, '-')).unwrap();
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

fn cannot_look_up(p: usize, i: usize, ampl_for_row: &[(usize, usize)]) -> bool {
    ampl_for_row[i].1 > ampl_for_row[p].1
}

fn get_best_u_special(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    scores_matrix: &HashMap<(char, char), i32>,
    ampl_for_row: &[(usize, usize)],
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut u = 0;
    let mut u_idx = 0;
    let mut first = true;
    for p in graph[i].1.iter() {
        if ampl_for_row[i].1 == ampl_for_row[*p].1 {
            let current_u = m[*p][j] + scores_matrix.get(&(graph[i].0, '-')).unwrap();
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
