use std::{
    cmp::{self, Ordering},
    collections::HashMap,
};

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
    o: i32,
    e: i32,
) -> i32 {
    let mut m = vec![vec![0; sequence.len()]; graph.len()];
    let mut x = vec![vec![0; sequence.len()]; graph.len()];
    let mut y = vec![vec![0; sequence.len()]; graph.len()];
    let mut path = vec![vec![('x', 0); sequence.len()]; graph.len()];

    for i in 0..graph.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {}
                (0, _) => {
                    y[0][j] = o + e * j as i32;
                    m[0][j] = y[0][j];

                    path[0][j] = ('L', j - 1);
                }
                (_, 0) => {
                    if graph[i].1.is_empty() {
                        x[i][0] = o + e * i as i32;
                        m[i][0] = x[i][0];

                        path[i][j] = ('U', i - 1);
                    } else {
                        let min_p = *graph[i].1.iter().min().unwrap() as i32;
                        x[i][j] = o + e * (min_p + 1);
                        m[i][j] = x[i][j];

                        path[i][j] = ('U', min_p as usize);
                    }
                }
                _ => {
                    let d;
                    if graph[i].1.is_empty() {
                        y[i][j] = cmp::max(y[i - 1][j] + e, m[i - 1][j] + o + e);
                        d = m[i - 1][j - 1]
                            + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                    } else {
                        let mut u_m_best = 0;
                        let mut u_y_best = 0;
                        let mut first = true;
                        let mut best_d = 0;
                        for p in graph[i].1.iter() {
                            let u_m_current = m[*p][j] + o + e;
                            let u_y_current = y[*p][j] + e;
                            let current_d = m[*p][j - 1];
                            if first {
                                u_m_best = u_m_current;
                                u_y_best = u_y_current;
                                best_d = current_d;
                                first = false;
                            }
                            if u_m_current > u_m_best {
                                u_m_best = u_m_current;
                            }
                            if u_y_current > u_y_best {
                                u_y_best = u_y_current;
                            }
                            if current_d > best_d {
                                best_d = current_d;
                            }
                        }
                        y[i][j] = cmp::max(u_m_best, u_y_best);
                        d = best_d + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                    }
                    x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);

                    let l = x[i][j];
                    let u = y[i][j];

                    match d.cmp(&l) {
                        Ordering::Less => match l.cmp(&u) {
                            Ordering::Less => {
                                m[i][j] = u;
                                path[i][j] = ('U', i - 1)
                            }
                            _ => {
                                m[i][j] = l;
                                path[i][j] = ('L', j - 1)
                            }
                        },
                        _ => match d.cmp(&u) {
                            Ordering::Less => {
                                m[i][j] = u;
                                path[i][j] = ('U', i - 1);
                            }
                            _ => {
                                m[i][j] = d;
                                if graph[i].0 == sequence[j] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1)
                                }
                            }
                        },
                    }
                }
            }
        }
    }
    for j in 0..sequence.len() {
        best_last_node(graph, &mut m, &mut path, j);
    }

    println!("Best alignment: {}", m[graph.len() - 1][sequence.len() - 1]);
    m[graph.len() - 1][sequence.len() - 1]

    //basic_output::write_align_poa(&path, sequence, graph);
}

fn best_last_node(
    graph: &[(char, Vec<usize>)],
    m: &mut [Vec<i32>],
    path: &mut [Vec<(char, usize)>],
    j: usize,
) {
    let mut best_align = 0;
    let mut best_idx = 0;
    let mut first = true;
    let m_len = m.len();

    for p in graph[graph.len() - 1].1.iter() {
        if first {
            best_align = m[*p][j];
            best_idx = *p;
            first = false;
        }
        if m[*p][j] > best_align {
            best_align = m[*p][j];
            best_idx = *p;
        }
    }
    m[m_len - 1][j] = best_align;
    path[m_len - 1][j] = ('F', best_idx);
}
#[cfg(test)]
mod tests {
    use crate::{graph, partial_order_alignment_global};
    use std::collections::HashMap;

    #[test]
    fn gap_align_gives_correct_result() {
        let sequence: Vec<char> = "$CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAAC"
            .chars()
            .collect();
        let graph = graph::get_linearization("./prova.gfa");
        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                    score_matrix.insert((*c1, *c2), 1);
                } else {
                    score_matrix.insert((*c1, *c2), -1);
                }
            }
        }
        let align_score = super::exec(&sequence, &graph, &score_matrix, -10, -2);
        assert_eq!(align_score, 20);
    }
    #[test]
    fn gap_align_gives_same_result_as_normal_if_o_zero() {
        let sequence: Vec<char> = "$CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAAC"
            .chars()
            .collect();
        let graph = graph::get_linearization("./prova.gfa");
        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T', '-'].iter() {
            for c2 in ['A', 'C', 'G', 'T', '-'].iter() {
                if c1 == c2 {
                    score_matrix.insert((*c1, *c2), 1);
                } else {
                    score_matrix.insert((*c1, *c2), -1);
                }
            }
        }
        let align_score = super::exec(&sequence, &graph, &score_matrix, 0, -1);
        assert_eq!(
            align_score,
            partial_order_alignment_global::exec(&sequence, &graph, &score_matrix)
        );
    }
}
