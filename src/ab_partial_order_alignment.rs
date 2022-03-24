use crate::basic_output;
use std::{
    cmp::{self},
    collections::HashMap,
    panic,
};

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
) -> i32 {
    let mut m = vec![vec![0; sequence.len()]; graph.len()];
    let mut path = vec![vec![('x', 0); sequence.len()]; graph.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0, 0); graph.len()];

    band_poa_align(
        sequence,
        graph,
        scores_matrix,
        ampl,
        &mut m,
        &mut path,
        &mut ampl_for_row,
    )
}
pub fn band_poa_align(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
    m: &mut Vec<Vec<i32>>,
    path: &mut Vec<Vec<(char, usize)>>,
    ampl_for_row: &mut Vec<(usize, usize)>,
) -> i32 {
    for i in 0..graph.len() - 1 {
        let (left, right) = set_left_right(ampl, i, graph, sequence);
        ampl_for_row[i] = (left, right);
        for j in left..right {
            match (i, j) {
                (0, 0) => {
                    m[i][j] = 0;
                    path[i][j] = ('O', 0);
                }
                (_, 0) => {
                    if graph[i].1.is_empty() {
                        m[i][j] = m[i - 1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        path[i][j] = ('U', i - 1);
                    } else {
                        let p = graph[i].1.iter().min().unwrap();
                        if j < ampl_for_row[*p].0 || j > ampl_for_row[*p].1 {
                            panic!("ampl wrong");
                        }
                        m[i][j] = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        path[i][j] = ('U', *p);
                    }
                }
                (0, _) => {
                    m[i][j] = m[i][j - 1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    path[i][j] = ('L', j - 1);
                }
                _ if j == left => {
                    // only look up or diagonal
                    if graph[i].1.is_empty() {
                        let u = m[i - 1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        let d = m[i - 1][j - 1]
                            + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        let u_idx = i - 1;
                        let d_idx = i - 1;
                        m[i][j] = *[d, u].iter().max().unwrap();
                        if m[i][j] == d {
                            if graph[i].0 == sequence[j] {
                                path[i][j] = ('D', d_idx);
                            } else {
                                path[i][j] = ('d', d_idx);
                            }
                        } else {
                            path[i][j] = ('U', u_idx);
                        }
                        
                    } else {
                        match (
                            get_best_d(graph, sequence, &m, scores_matrix, &ampl_for_row, i, j),
                            get_best_u(graph, &m, scores_matrix, &ampl_for_row, i, j),
                        ) {
                            (Some((d, d_idx)), Some((u, u_idx))) => {
                                m[i][j] = *[d, u].iter().max().unwrap();
                                if m[i][j] == d {
                                    if graph[i].0 == sequence[j] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx);
                                    }
                                } else {
                                    path[i][j] = ('U', u_idx);
                                }
                            }
                            _ => {}
                        }
                    }
                }
                _ if j == right - 1 =>{
                    // only left or d
                    let l = m[i][j - 1] + scores_matrix.get(&('-', sequence[j])).unwrap();

                    if graph[i].1.is_empty() {
                        let d = m[i - 1][j - 1]
                            + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        let d_idx = i - 1;
                        m[i][j] = *[d, l].iter().max().unwrap();
                        if m[i][j] == d {
                            if graph[i].0 == sequence[j] {
                                path[i][j] = ('D', d_idx);
                            } else {
                                path[i][j] = ('d', d_idx);
                            }
                        } else {
                            path[i][j] = ('L', j-1);
                        }
                        
                    } else {
                        match 
                            get_best_d(graph, sequence, &m, scores_matrix, &ampl_for_row, i, j)
                         {
                            Some((d, d_idx)) => {
                                m[i][j] = *[d, l].iter().max().unwrap();
                                if m[i][j] == d {
                                    if graph[i].0 == sequence[j] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx);
                                    }
                                } else {
                                    path[i][j] = ('L', j-1);
                                }
                            }
                            _ => {}
                        }
                    }
                }
                _ => {
                    // d u or l
                    let l = m[i][j - 1] + scores_matrix.get(&('-', sequence[j])).unwrap();

                    if graph[i].1.is_empty() {
                        let u = m[i - 1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        let d = m[i - 1][j - 1]
                            + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        let u_idx = i - 1;
                        let d_idx = i - 1;
                        m[i][j] = *[d, u, l].iter().max().unwrap();
                        if m[i][j] == d {
                            if graph[i].0 == sequence[j] {
                                path[i][j] = ('D', d_idx);
                            } else {
                                path[i][j] = ('d', d_idx);
                            }
                        } else if m[i][j] == l {
                            path[i][j] = ('L', j - 1);
                        } else {
                            path[i][j] = ('U', u_idx);
                        }
                    } else {
                        match (
                            get_best_d(graph, sequence, &m, scores_matrix, &ampl_for_row, i, j),
                            get_best_u(graph, &m, scores_matrix, &ampl_for_row, i, j),
                        ) {
                            (Some((d, d_idx)), Some((u, u_idx))) => {
                                m[i][j] = *[d, u, l].iter().max().unwrap();
                                if m[i][j] == d {
                                    if graph[i].0 == sequence[j] {
                                        path[i][j] = ('D', d_idx);
                                    } else {
                                        path[i][j] = ('d', d_idx);
                                    }
                                } else if m[i][j] == l {
                                    path[i][j] = ('L', j - 1);
                                } else {
                                    path[i][j] = ('U', u_idx);
                                }
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
    }
    for j in 0..m[0].len() {
        best_last_node(graph, m, path, j);
    }
    match ampl_is_enough(&path, &ampl_for_row) {
        true => {
            println!("{}", m[graph.len() - 1][sequence.len() - 1]);
            basic_output::write_align_banded_poa(&path, sequence, graph);
            m[graph.len() - 1][sequence.len() - 1]
        }
        false => band_poa_align(
            sequence,
            graph,
            scores_matrix,
            ampl * 2,
            m,
            path,
            ampl_for_row,
        ),
    }
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

fn set_left_right(
    ampl: usize,
    i: usize,
    graph: &[(char, Vec<usize>)],
    sequence: &[char],
) -> (usize, usize) {
    let mut left = 0;
    let mut right = sequence.len();
    if i == 0 {
        right = cmp::min(ampl / 2, sequence.len());
    } else {
        if graph[i].1.is_empty() {
            if i < ampl / 2 {
                left = 0;
            } else {
                left = i - ampl / 2
            }
            right = cmp::min(i + ampl / 2, sequence.len());
            if right <= left {
                (left, right) = (right - 1, right);
            }
        } else {
            let mut first = true;
            for p in graph[i].1.iter() {
                let (current_l, current_r);
                if p + 1 < ampl / 2 {
                    current_l = 0;
                } else {
                    current_l = p + 1 - ampl / 2;
                }
                current_r = cmp::min(p + 1 + ampl / 2, sequence.len());
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
            if right <= left {
                (left, right) = (right - 1, right);
            }
        }
    }
    (left, right)
}

fn ampl_is_enough(path: &[Vec<(char, usize)>], ampl_for_row: &Vec<(usize, usize)>) -> bool {
    let mut col = path[0].len() - 1;
    let mut row = path[path.len() - 1][path[0].len() - 1].1 as usize;

    while path[row][col].0 != 'O' {
        //reached end of path, no need to continue
        if col == 0 {
            return true;
        }

        if col == ampl_for_row[row].0
            || (col == ampl_for_row[row].1 - 1 && !(col == path[0].len() - 1))
        {
            // != from path_len because couldn't go larger
            if path[row][col].0 == 'D' {
                row = path[row][col].1 as usize;
                col -= 1;
                // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            match path[row][col].0 {
                'D' | 'd' => {
                    row = path[row][col].1 as usize;
                    col -= 1;
                }
                'U' => {
                    row = path[row][col].1 as usize;
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
    ampl_for_row: &Vec<(usize, usize)>,
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut d = 0;
    let mut d_idx = 0;
    let mut first = true;
    for p in graph[i].1.iter() {
        if j - 1 >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 + 1 {
            let current_d = m[*p][j - 1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
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
    ampl_for_row: &Vec<(usize, usize)>,
    i: usize,
    j: usize,
) -> Option<(i32, usize)> {
    let mut u = 0;
    let mut u_idx = 0;
    let mut first = true;
    for p in graph[i].1.iter() {
        if j - 1 >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 + 1 {
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
#[cfg(test)]
mod tests {
    use crate::graph;
    use std::collections::HashMap;

    #[test]
    fn ab_align_gives_correct_result() {
        let sequence: Vec<char> = "$CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG"
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
        let align_score = super::exec(&sequence, &graph, &score_matrix, 5);
        assert_eq!(align_score, 48);
    }
}
