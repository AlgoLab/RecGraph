use std::{cmp::Ordering, collections::HashMap};

use crate::basic_output;

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
) {
    let mut m = vec![vec![0; graph.len()]; sequence.len()];
    let mut path = vec![vec![('x', 0 as i32); graph.len()]; sequence.len()];

    for i in 0..sequence.len() {
        for j in 0..graph.len() {
            match (i, j) {
                (0, 0) => path[i][j] = ('O', 0),
                (0, _) => {
                    // internal node in gfa
                    if graph[j].1.len() == 0 {
                        m[i][j] = m[i][j - 1] + scores_matrix.get(&('-', graph[j].0)).unwrap();
                        path[i][j] = ('L', j as i32 - 1);
                    } else {
                        let mut best_al = 0;
                        let mut first = true;
                        // search in other predecessor to find better alignment
                        for p in graph[j].1.iter() {
                            let p_align = m[i][*p];
                            if first {
                                best_al = p_align;
                                path[i][j] = ('L', *p as i32);
                                first = false;
                            }
                            if p_align > best_al {
                                best_al = p_align;
                                path[i][j] = ('L', *p as i32);
                            }
                        }
                        m[i][j] = best_al + scores_matrix.get(&('-', graph[j].0)).unwrap();
                    }
                }
                (_, 0) => {
                    m[i][j] = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
                    path[i][j] = ('U', 0)
                }
                _ => {
                    // same as (0, _)
                    if graph[j].1.len() == 0 {
                        let best_d = m[i - 1][j - 1]
                            + scores_matrix.get(&(graph[j].0, sequence[i])).unwrap();
                        let best_u = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
                        let best_l = m[i][j - 1] + scores_matrix.get(&('-', graph[j].0)).unwrap();

                        match best_d.cmp(&best_u) {
                            Ordering::Less => match best_u.cmp(&best_l) {
                                Ordering::Less => {
                                    m[i][j] = best_l;
                                    path[i][j] = ('L', j as i32 - 1);
                                }
                                _ => {
                                    m[i][j] = best_u;
                                    path[i][j] = ('U', 0);
                                }
                            },
                            _ => match best_d.cmp(&best_l) {
                                Ordering::Less => {
                                    m[i][j] = best_l;
                                    path[i][j] = ('L', j as i32 - 1);
                                }
                                _ => {
                                    m[i][j] = best_d;
                                    if sequence[i] == graph[j].0 {
                                        path[i][j] = ('D', j as i32 - 1)
                                    } else {
                                        path[i][j] = ('d', j as i32 - 1)
                                    }
                                }
                            },
                        }
                    } else {
                        let mut best_d = 0;
                        let mut best_d_idx = 0;

                        let best_u = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();

                        let mut best_l = 0;
                        let mut best_l_idx = 0;

                        let mut first = true;

                        for p in graph[j].1.iter() {
                            let p_d = m[i - 1][*p]
                                + scores_matrix.get(&(graph[j].0, sequence[i])).unwrap();
                            let p_l = m[i][*p] + scores_matrix.get(&('-', graph[j].0)).unwrap();
                            if first {
                                best_d = p_d;
                                best_d_idx = *p;
                                best_l = p_l;
                                best_l_idx = *p;
                                first = false;
                            }
                            if p_d > best_d {
                                best_d = p_d;
                                best_d_idx = *p
                            }

                            if p_l > best_l {
                                best_l = p_l;
                                best_l_idx = *p;
                            }
                        }
                        match best_d.cmp(&best_u) {
                            Ordering::Less => match best_u.cmp(&best_l) {
                                Ordering::Less => {
                                    m[i][j] = best_l;
                                    path[i][j] = ('L', best_l_idx as i32);
                                }
                                _ => {
                                    m[i][j] = best_u;
                                    path[i][j] = ('U', 0);
                                }
                            },
                            _ => match best_d.cmp(&best_l) {
                                Ordering::Less => {
                                    m[i][j] = best_l;
                                    path[i][j] = ('L', best_l_idx as i32);
                                }
                                _ => {
                                    m[i][j] = best_d;
                                    if sequence[i] == graph[j].0 {
                                        path[i][j] = ('D', best_d_idx as i32)
                                    } else {
                                        path[i][j] = ('d', best_d_idx as i32)
                                    }
                                }
                            },
                        }
                    }
                }
            }
        }
    }
    
    println!("Best alignment: {}", m[sequence.len() - 1][graph.len() - 1]);
    
    basic_output::write_align_poa(&path, sequence, graph);
}
