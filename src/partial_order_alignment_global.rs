use std::{cmp::Ordering, collections::HashMap};

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
) {
    let mut m = vec![vec![0; graph.len()]; sequence.len()];

    for i in 0..sequence.len() {
        for j in 0..graph.len() {
            match (i, j) {
                (0, 0) => {}
                (0, _) => {
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
                (_, 0) => {
                    m[i][j] = m[i - 1][j] + scores_matrix.get(&('-', sequence[i])).unwrap();
                }
                _ => {
                    //j-1 always a predecessor, i only has i-1 as predecessor
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
    println!("Best alignment: {}", m[sequence.len() - 1][graph.len() - 1])
}
