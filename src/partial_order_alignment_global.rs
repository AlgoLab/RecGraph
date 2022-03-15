use std::{cmp::Ordering, collections::HashMap};

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
) {
    let mut m = vec![vec![0; graph.len()]; sequence.len()];

    for i in 0..sequence.len() {
        for j in 0..graph.len() {
            if i == 0 && j > 0 {
                if graph[j].1.len() == 0 {
                    m[i][j] = m[i][j - 1] + scores_matrix.get(&('-', graph[j].0)).unwrap();
                } else {
                    // j-1 is always a predecessor
                    let mut best_al = m[i][j - 1];
                    // search in other predecessor to find better alignment
                    for p in graph[j].1.iter() {
                        let p_align = m[i][*p];
                        if p_align > best_al {
                            best_al = p_align
                        }
                    }
                    m[i][j] = best_al + scores_matrix.get(&('-', graph[j].0)).unwrap();
                }
            } else if j == 0 && i > 0 {
                m[i][j] = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
            } else if i > 0 && j > 0 {
                if graph[j].1.len() == 0 {
                    //only j-1 as predex, normal global alignment
                    let d =
                        m[i - 1][j - 1] + scores_matrix.get(&(graph[j].0, sequence[i])).unwrap();
                    let u = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
                    let l = m[i][j - 1] + scores_matrix.get(&('-', graph[j].0)).unwrap();

                    match d.cmp(&u) {
                        Ordering::Less => match u.cmp(&l) {
                            Ordering::Less => m[i][j] = l,
                            _ => m[i][j] = u,
                        },
                        _ => match d.cmp(&l) {
                            Ordering::Less => m[i][j] = l,
                            _ => m[i][j] = d,
                        },
                    }
                } else {
                    let mut best_d =
                        m[i - 1][j - 1] + scores_matrix.get(&(graph[j].0, sequence[i])).unwrap();
                    let best_u = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
                    let mut best_l = m[i][j - 1] + scores_matrix.get(&('-', graph[j].0)).unwrap();
                    // search in other predecessor to find better alignment
                    for p in graph[j].1.iter() {
                        let p_d =
                            m[i - 1][*p] + scores_matrix.get(&(graph[j].0, sequence[i])).unwrap();
                        if p_d > best_d {
                            best_d = p_d
                        }
                        let p_l = m[i][*p] + scores_matrix.get(&('-', graph[j].0)).unwrap();
                        if p_l > best_l {
                            best_l = p_l
                        }
                    }
                    match best_d.cmp(&best_u) {
                        Ordering::Less => match best_u.cmp(&best_l) {
                            Ordering::Less => m[i][j] = best_l,
                            _ => m[i][j] = best_u,
                        },
                        _ => match best_d.cmp(&best_l) {
                            Ordering::Less => m[i][j] = best_l,
                            _ => m[i][j] = best_d,
                        },
                    }
                }
            }
        }
    }
    println!("Best alignment: {}", m[sequence.len()-1][graph.len()-1])
}
