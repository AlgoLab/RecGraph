use std::{
    cmp::{self, Ordering},
    collections::HashMap,
};
//FIXME: band doesn't work

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
    o: i32,
    e: i32,
) -> i32 {
    //colonne caratteri di sequence, righe caratteri di graph
    let mut m = vec![vec![0; ampl]; graph.len()];
    let mut x = vec![vec![0; ampl]; graph.len()];
    let mut y = vec![vec![0; ampl]; graph.len()];
    let mut path = vec![vec![('x', 0); ampl]; graph.len()];

    m[0][(ampl / 2) as usize] = 0;
    path[0][(ampl / 2) as usize] = ('O', 0);

    for j in ampl / 2 + 1..ampl {
        // prima riga
        if j - (ampl / 2) < sequence.len() {
            y[0][j] = o + e * (j - ampl / 2) as i32;
            m[0][j] = y[0][j];

            path[0][j] = ('L', j - 1);
        }
    }
    for i in 1..graph.len() - 1 {
        for j in 0..ampl {
            if i + j < sequence.len() + ampl / 2 && i + j >= ampl / 2 {
                if i + j == ampl / 2 {
                    // elementi con j = 0 matrice mn
                    if graph[i].1.is_empty() {
                        x[i][j] = o + e * i as i32;
                        m[i][j] = x[i][j];

                        path[i][j] = ('U', i - 1);
                    } else {
                        //take smallest predecessor of i
                        let min_p = *graph[i].1.iter().min().unwrap() as i32;
                        x[i][j] = o + e * (min_p + 1);
                        m[i][j] = x[i][j];

                        path[i][j] = ('U', min_p as usize);
                    }
                } else if j == 0 {
                    // bordo inf banda, prendo da diag o da sopra
                    let index_of_seq = i + j - ampl / 2;

                    if graph[i].1.is_empty() {
                        // come se fosse sequenza
                        // set y
                        y[i][j] = cmp::max(y[i - 1][j + 1] + e, m[i - 1][j + 1] + o + e);
                        /*
                        if y[i][j] != m[i - 1][j + 1] + o + e {
                            path_y[i - 1][j + 1] = 'Y'
                        }
                        */

                        // set m
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
                        let u = y[i][j];
                        match d.cmp(&u) {
                            cmp::Ordering::Less => {
                                m[i][j] = u;

                                path[i][j] = ('U', i - 1);
                            }
                            _ => {
                                m[i][j] = d;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1);
                                }
                            }
                        }
                    } else {
                        // diversi predecessori
                        // set y
                        let (u_best, u_idx) = get_best_u_pred(graph, &m, &y, i, j, o, e);
                        y[i][j] = u_best;

                        // set m
                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);
                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();

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
                        // set x
                        x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);

                        // set m
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
                        let l = x[i][j];
                        match d.cmp(&l) {
                            Ordering::Less => {
                                m[i][j] = l;
                                path[i][j] = ('L', j - 1);
                            }
                            _ => {
                                m[i][j] = d;
                                if graph[i].0 == sequence[index_of_seq] {
                                    path[i][j] = ('D', i - 1);
                                } else {
                                    path[i][j] = ('d', i - 1);
                                }
                            }
                        }
                    } else {
                        // diversi predecessori
                        // set x
                        x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);

                        // set m
                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);
                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();

                        let l_best = x[i][j];

                        match d_best.cmp(&l_best) {
                            Ordering::Less => {
                                m[i][j] = l_best;
                                path[i][j] = ('L', j - 1);
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
                        //set x
                        x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);

                        //set y
                        y[i][j] = cmp::max(y[i - 1][j + 1] + e, m[i - 1][j + 1] + o + e);

                        //set m
                        let d = m[i - 1][j]
                            + score_matrix
                                .get(&(sequence[index_of_seq], graph[i].0))
                                .unwrap();
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
                                    if graph[i].0 == sequence[index_of_seq] {
                                        path[i][j] = ('D', i - 1);
                                    } else {
                                        path[i][j] = ('d', i - 1)
                                    }
                                }
                            },
                        }
                    } else {
                        // diversi predecessori
                        // set y
                        let (u_best, u_idx) = get_best_u_pred(graph, &m, &y, i, j, o, e);
                        y[i][j] = u_best;

                        // set x
                        x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);

                        // set m

                        let (mut d_best, d_idx) = get_best_d_pred(graph, &m, i, j);

                        d_best += score_matrix
                            .get(&(sequence[index_of_seq], graph[i].0))
                            .unwrap();

                        let l_best = x[i][j];

                        match d_best.cmp(&l_best) {
                            Ordering::Less => match l_best.cmp(&u_best) {
                                Ordering::Less => {
                                    m[i][j] = u_best;
                                    path[i][j] = ('U', u_idx)
                                }
                                _ => {
                                    m[i][j] = l_best;
                                    path[i][j] = ('L', j - 1)
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
    // virtual last node of the graph, take best value from terminal nodes
    for j in 0..ampl / 2 {
        best_last_node(graph, &mut m, &mut path, j);
    }

    match ampl_is_enough(&path, (sequence.len() - 1) + (ampl / 2) - (graph.len() - 1)) {
        true => {
            println!(
                "AB GAP Best alignment: {}",
                m[graph.len() - 1][(sequence.len() - 1) + (ampl / 2) - (graph.len() - 1)]
            );
            return m[graph.len() - 1][(sequence.len() - 1) + (ampl / 2) - (graph.len() - 1)];
        }
        false => exec(sequence, graph, score_matrix, ampl * 2 + 1, o, e),
    }

    //basic_output::write_align_ab_poa(&path, sequence, graph);
}

fn get_best_d_pred(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    i: usize,
    j: usize,
) -> (i32, usize) {
    let mut d_best = 0;
    let mut d_idx = 0;

    let mut first = true;
    for p in graph[i].1.iter() {
        let d_align = m[*p][j + (i - p) - 1];
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
    (d_best, d_idx)
}

fn get_best_u_pred(
    graph: &[(char, Vec<usize>)],
    m: &[Vec<i32>],
    y: &[Vec<i32>],
    i: usize,
    j: usize,
    o: i32,
    e: i32,
) -> (i32, usize) {
    let mut u_m_best = 0;
    let mut u_y_best = 0;
    let mut u_m_idx = 0;
    let mut u_y_idx = 0;
    let mut first = true;

    for p in graph[i].1.iter() {
        let current_m = m[*p][j + (i - p)] + o + e;
        let current_y = y[*p][j + (i - p)] + e;

        if first {
            u_m_best = current_m;
            u_y_best = current_y;

            u_m_idx = *p;
            u_y_idx = *p;
            first = false;
        }
        if current_m > u_m_best {
            u_m_best = current_m;
            u_m_idx = *p;
        }
        if current_y > u_y_best {
            u_y_best = current_y;
            u_y_idx = *p
        }
    }
    match u_m_best.cmp(&u_y_best) {
        Ordering::Less => return (u_y_best, u_y_idx),
        _ => return (u_m_best, u_m_idx),
    };
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
    let i = m.len() - 1;

    for p in graph[graph.len() - 1].1.iter() {
        let delta = i - p;
        let current_align = m[*p][j + delta];
        if first {
            best_align = current_align;
            best_idx = *p;
            first = false;
        }
        if current_align > best_align {
            best_align = current_align;
            best_idx = *p;
        }
    }
    m[m_len - 1][j] = best_align;
    path[m_len - 1][j] = ('F', best_idx);
}

fn ampl_is_enough(path: &[Vec<(char, usize)>], start_col: usize) -> bool {
    let mut col = start_col;
    let mut row = path[path.len() - 1][start_col].1;

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
                    let delta = row - path[row][col].1 as usize;
                    row = path[row][col].1 as usize;
                    col += delta - 1;
                }
                'U' => {
                    let delta = row - path[row][col].1;
                    row = path[row][col].1;
                    col += delta;
                }
                'L' => {
                    col = path[row][col].1;
                }
                _ => panic!("ampl_is_enough panic"),
            }
        }
    }
    true
}
#[cfg(test)]
mod tests {
    use crate::graph;
    use std::collections::HashMap;
    #[test]
    fn if_g_open_penalty_null_equal_to_normal_poa() {
        let graph = graph::get_linearization("./prova.gfa");
        let seq: Vec<char> = "$CAAATAAG".chars().collect();

        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                    score_matrix.insert((*c1, *c2), 2);
                } else {
                    score_matrix.insert((*c1, *c2), -4);
                }
            }
        }
        let align_score = super::exec(
            &seq,
            &graph,
            &score_matrix,
            (graph.len() - seq.len()) * 2 + 1,
            0,
            -4,
        );
        assert_eq!(align_score, -144);
    }
    #[test]
    fn ab_gap_gives_correct_result() {
        let graph = graph::get_linearization("./prova.gfa");
        let seq: Vec<char> = "$CAAATAAG".chars().collect();

        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                    score_matrix.insert((*c1, *c2), 2);
                } else {
                    score_matrix.insert((*c1, *c2), -4);
                }
            }
        }
        let align_score = super::exec(
            &seq,
            &graph,
            &score_matrix,
            (graph.len() - seq.len()) * 2 + 1,
            200,
            -4,
        );
        assert_eq!(align_score, -344);
    }
}
