use std::{cmp::Ordering, collections::HashMap};

use crate::basic_output;

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: usize
) {
    //colonne caratteri di sequence, righe caratteri di graph
    let mut m = vec![vec![0; ampl]; graph.len()];
    let mut path = vec![vec![('x', 0 as i32); ampl]; graph.len()];

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
                    if graph[i].1.len() == 0 {
                        m[i][j] = m[i - 1][j + 1] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        path[i][j] = ('U', i as i32 - 1);
                    } else {
                        let mut best_al = 0;
                        let mut first = true;
                        for p in graph[i].1.iter() {
                            let p_align = m[*p][j+1];
                            if first {
                                best_al = p_align;
                                path[i][j] = ('U', *p as i32);
                                first = false;
                            }
                        }
                        m[i][j] = best_al + score_matrix.get(&(graph[i].0, '-')).unwrap();
                    }
                    
                } else if j == 0 {
                    // bordo inf banda, prendo da diag o da sopra
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.len() == 0 { // come se fosse sequenza
                        let d = m[i - 1][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
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
                    } else { // diversi predecessori
                        let mut d_best = 0;
                        let mut u_best = 0;
                        let mut d_idx = 0;
                        let mut u_idx = 0;
                        let mut first = true;
                        for p in graph[i].1.iter() {
                            let u_align = m[*p][j+1]+ score_matrix.get(&('-', graph[i].0)).unwrap();
                            let d_align = m[*p][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
                            if first {
                                d_best = d_align;
                                u_best = u_align;
                                d_idx = *p;
                                u_idx = *p;
                                first = false;
                            }
                            if d_align > d_best {
                                d_best = d_align;
                                d_idx = *p;
                            }
                            if u_align > u_best {
                                u_best = u_align;
                                u_idx = *p;
                            }
                        }
                        match d_best.cmp(&u_best) {
                            Ordering::Less => {
                                m[i][j] = u_best;
                                path[i][j] = ('U', u_idx as i32)
                            }
                            _ => {
                                m[i][j] = d_best;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', d_idx as i32)
                                } else {
                                    path[i][j] = ('d', d_idx as i32)
                                }
                            }
                        }
                    }
                } else if j == ampl - 1 {
                    // bordo sup banda
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.len() == 0 { // come se fosse sequenza
                        let d = m[i - 1][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
                        let l = m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();
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
                    } else { // diversi predecessori
                        let mut d_best = 0;
                        let mut d_idx = 0;
                        let mut first = true;

                        for p in graph[i].1.iter() {
                            let d_align = m[*p][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
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
                        let l_best = m[i][j-1] + score_matrix.get(&(sequence[index_of_seq],'-')).unwrap();

                        match d_best.cmp(&l_best) {
                            Ordering::Less => {
                                m[i][j] = l_best;
                                path[i][j] = ('L', 0)
                            }
                            _ => {
                                m[i][j] = d_best;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', d_idx as i32)
                                } else {
                                    path[i][j] = ('d', d_idx as i32)
                                }
                            }
                        }
                    }
                } else {
                    // celle interne banda
                    let index_of_seq = i + j - ampl / 2;
                    
                    if graph[i].1.len() == 0 {
                        let d = m[i - 1][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
                        let u = m[i - 1][j + 1] + score_matrix.get(&(graph[i].0, '-')).unwrap();
                        let l = m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();

                        match d.cmp(&l) {
                            Ordering::Less => match l.cmp(&u) {
                                Ordering::Less => {
                                    m[i][j] = u;
                                    path[i][j] = ('U', i as i32-1)
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
                        let mut d_best = 0;
                        let mut u_best = 0;
                        let l_best = m[i][j - 1] + score_matrix.get(&(sequence[index_of_seq], '-')).unwrap();

                        let mut d_idx = 0;
                        let mut u_idx = 0;

                        let mut first = true;
                        for p in graph[i].1.iter() {
                            let u_align = m[*p][j+1]+ score_matrix.get(&('-', graph[i].0)).unwrap();
                            let d_align = m[*p][j] + score_matrix.get(&(sequence[index_of_seq], graph[i].0)).unwrap();
                            if first {
                                d_best = d_align;
                                u_best = u_align;
                                d_idx = *p;
                                u_idx = *p;
                                first = false;
                            }
                            if d_align > d_best {
                                d_best = d_align;
                                d_idx = *p;
                            }
                            if u_align > u_best {
                                u_best = u_align;
                                u_idx = *p;
                            }
                        }
                        match d_best.cmp(&l_best) {
                            Ordering::Less => match l_best.cmp(&u_best) {
                                Ordering::Less => {
                                    m[i][j] = u_best;
                                    path[i][j] = ('U', u_idx as i32)
                                }
                                _ => {
                                    m[i][j] = l_best;
                                    path[i][j] = ('L', 0)
                                }
                            },
                            _ => match d_best.cmp(&u_best) {
                                Ordering::Less => {
                                    m[i][j] = u_best;
                                    path[i][j] = ('U', u_idx as i32);
                                }
                                _ => {
                                    m[i][j] = d_best;
                                    if graph[i].0 == sequence[index_of_seq] {
                                        path[i][j] = ('D', d_idx as i32);
                                    } else {
                                        path[i][j] = ('d', d_idx as i32)
                                    }
                                }
                            },
                        }
                    }
                }
            }
        }
    }
    

    println!("Best alignment: {}", m[graph.len()-1][(sequence.len() - 1) + (ampl/2) - (graph.len()-1)]);
    
    //basic_output::write_align_poa(&path, sequence, graph);
}
