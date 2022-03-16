use std::{cmp::Ordering, collections::HashMap};

use crate::basic_output;

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
) {
    //colonne caratteri di sequence, righe caratteri di graph
    let mut m = vec![vec![-10000; ampl]; graph.len()];
    let mut path = vec![vec![('x', 0); ampl]; graph.len()];

    m[0][(ampl / 2) as usize] = 0;
    path[0][(ampl / 2) as usize] = ('O', 0);

    for j in ampl / 2 + 1..ampl {
        // prima riga
        if j - (ampl / 2) < sequence.len() {
            m[0][j] = m[0][j - 1] + score_matrix.get(&('-', sequence[j - ampl / 2])).unwrap();
            path[0][j] = ('L', 0);
        }
    }
    for i in 1..graph.len() {
        for j in 0..ampl {
            if i + j < sequence.len() + ampl / 2 && i + j >= ampl / 2 {
                if i + j == ampl / 2 {
                    // elementi con j = 0 matrice mn
                    if graph[i].1.is_empty() {
                        m[i][j] = m[i - 1][j + 1] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        path[i][j] = ('U', i as i32 - 1);
                    } else {
                        let (u_best, u_idx) = get_best_u_pred(graph, &m, i, j);
                        m[i][j] = u_best + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        path[i][j] = ('U', u_idx);
                    }
                } else if j == 0 {
                    // bordo inf banda, prendo da diag o da sopra
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.is_empty() {
                        // come se fosse sequenza
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
                        let u = m[i - 1][j + 1] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        match d.cmp(&u) {
                            Ordering::Less => {
                                m[i][j] = u;
                                path[i][j] = ('U', i as i32 - 1)
                            }
                            _ => {
                                m[i][j] = d;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', i as i32 - 1)
                                } else {
                                    path[i][j] = ('d', i as i32 - 1)
                                }
                            }
                        }
                    } else {
                        // diversi predecessori
                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);
                        let (mut u_best, u_idx) = get_best_u_pred(graph, &m, i, j);

                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();
                        u_best += score_matrix.get(&('-', graph[i].0)).unwrap();

                        match d_best.cmp(&u_best) {
                            Ordering::Less => {
                                m[i][j] = u_best;
                                path[i][j] = ('U', u_idx)
                            }
                            _ => {
                                m[i][j] = d_best;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', d_idx)
                                } else {
                                    path[i][j] = ('d', d_idx)
                                }
                            }
                        }
                    }
                } else if j == ampl - 1 {
                    // bordo sup banda
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.is_empty() {
                        // come se fosse sequenza
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
                        let l =
                            m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();
                        match d.cmp(&l) {
                            Ordering::Less => {
                                m[i][j] = l;
                                path[i][j] = ('L', 0)
                            }
                            _ => {
                                m[i][j] = d;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', i as i32 - 1)
                                } else {
                                    path[i][j] = ('d', i as i32 - 1)
                                }
                            }
                        }
                    } else {
                        // diversi predecessori
                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);
                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();

                        let l_best =
                            m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();

                        match d_best.cmp(&l_best) {
                            Ordering::Less => {
                                m[i][j] = l_best;
                                path[i][j] = ('L', 0)
                            }
                            _ => {
                                m[i][j] = d_best;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', d_idx)
                                } else {
                                    path[i][j] = ('d', d_idx)
                                }
                            }
                        }
                    }
                } else {
                    // celle interne banda
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.is_empty() {
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
                        let u = m[i - 1][j + 1] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        let l =
                            m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();

                        match d.cmp(&l) {
                            Ordering::Less => match l.cmp(&u) {
                                Ordering::Less => {
                                    m[i][j] = u;
                                    path[i][j] = ('U', i as i32 - 1)
                                }
                                _ => {
                                    m[i][j] = l;
                                    path[i][j] = ('L', 0)
                                }
                            },
                            _ => match d.cmp(&u) {
                                Ordering::Less => {
                                    m[i][j] = u;
                                    path[i][j] = ('U', i as i32 - 1);
                                }
                                _ => {
                                    m[i][j] = d;
                                    if graph[i].0 == sequence[index_of_seq] {
                                        path[i][j] = ('D', i as i32 - 1);
                                    } else {
                                        path[i][j] = ('d', i as i32 - 1)
                                    }
                                }
                            },
                        }
                    } else {
                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);
                        let (mut u_best, u_idx) = get_best_u_pred(graph, &m, i, j);

                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();
                        u_best += score_matrix.get(&('-', graph[i].0)).unwrap();

                        let l_best =
                            m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();

                        match d_best.cmp(&l_best) {
                            Ordering::Less => match l_best.cmp(&u_best) {
                                Ordering::Less => {
                                    m[i][j] = u_best;
                                    path[i][j] = ('U', u_idx)
                                }
                                _ => {
                                    m[i][j] = l_best;
                                    path[i][j] = ('L', 0)
                                }
                            },
                            _ => match d_best.cmp(&u_best) {
                                Ordering::Less => {
                                    m[i][j] = u_best;
                                    path[i][j] = ('U', u_idx);
                                }
                                _ => {
                                    m[i][j] = d_best;
                                    if graph[i].0 == sequence[index_of_seq] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx)
                                    }
                                }
                            },
                        }
                    }
                }
            }
        }
    }
    match ampl_is_enough(&path, (sequence.len() - 1) + (ampl / 2) - (graph.len() - 1)) {
        true => {
            println!(
                "AB Best alignment: {}",
                m[graph.len() - 1][(sequence.len() - 1) + (ampl / 2) - (graph.len() - 1)]
            );
        }
        false => {
            exec(sequence, graph, score_matrix, ampl * 2 + 1);
        }
    }

    basic_output::write_align_ab_poa(&path, sequence, graph);
}

fn get_best_d_pred(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    i: usize,
    j: usize,
) -> (i32, i32) {
    let mut d_best = 0;
    let mut d_idx = 0;

    let mut first = true;
    for p in graph[i].1.iter() {
        let d_align = m[*p][j];
        if first {
            d_best = d_align;
            d_idx = *p;
            first = false;
        }
        if d_align > d_best {
            d_best = d_align;
            d_idx = *p;
        }
    }
    (d_best, d_idx as i32)
}

fn get_best_u_pred(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    i: usize,
    j: usize,
) -> (i32, i32) {
    let mut u_best = 0;
    let mut u_idx = 0;

    let mut first = true;
    for p in graph[i].1.iter() {
        let p_align = m[*p][j + (i - p)]; // not j + 1 because +1 for every line up
        if first {
            u_best = p_align;
            u_idx = *p;
            first = false;
        }
        if p_align > u_best {
            u_best = p_align;
            u_idx = *p;
        }
    }
    (u_best, u_idx as i32)
}

fn ampl_is_enough(path: &[Vec<(char, i32)>], start_col: usize) -> bool {
    let mut row = path.len() - 1;
    let mut col = start_col;
    let col_number = path[0].len();

    while path[row][col].0 != 'O' {
        if col == 0 || col == col_number - 1 {
            if path[row][col].0 == 'D' {
                row -= 1; // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            match path[row][col].0 {
                'D' | 'd' => {
                    row = path[row][col].1 as usize;
                }
                'U' => {
                    let delta = row - path[row][col].1 as usize;
                    row = path[row][col].1 as usize;
                    col += delta;
                }
                'L' => {
                    col -= 1;
                }
                _ => panic!("ampl_is_enough panic"),
            }
        }
    }
    true
}
