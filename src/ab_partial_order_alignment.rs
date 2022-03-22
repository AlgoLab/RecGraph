use std::{cmp::{Ordering, self}, collections::HashMap, panic};

use crate::basic_output;

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
) {
    let mut m = vec![vec![0; sequence.len()]; graph.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0,0); graph.len()];
    let mut path = vec![vec![('x', 0); sequence.len()]; graph.len()];
    let ampl = 51; // ampl will be passed as an argument
    for i in 0..graph.len() - 1 {
        let (left, right) = set_left_right(ampl, i, graph, sequence);
        ampl_for_row[i] = (left, right);
        for j in left..right {
            match (i, j) {
                (0,0) => {
                    path[i][j] = ('O', 0);
                },
                (_,0) => {
                    if graph[i].1.is_empty() {
                        m[i][j] = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        path[i][j] = ('U', i-1);
                    } else {
                        let p = graph[i].1.iter().min().unwrap();
                        if j < ampl_for_row[*p].0 || j > ampl_for_row[*p].1 {
                            panic!("ampl wrong");
                        }
                        m[i][j] = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        path[i][j] = ('U', *p);
                    }
                },
                (0,_) => {
                    m[i][j] = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    path[i][j] = ('L', j-1);
                },
                _ if j == left => {
                    // only look up or diagonal
                    let mut u = 0;
                    let mut d= 0;
                    let mut u_idx = 0;
                    let mut d_idx = 0;

                    if graph[i].1.is_empty() {
                        u = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        u_idx = i-1;
                        d_idx = i-1;
                    } else {
                        let mut first = true;
                        for p in graph[i].1.iter(){
                            if j >= ampl_for_row[*p].0 + 1 && j <= ampl_for_row[*p].1 + 1 {
                                let current_u = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                                let current_d =  m[*p][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                                if first {
                                    first = false;
                                    u = current_u;
                                    d = current_d;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_u > u {
                                    u = current_u;
                                    u_idx = *p;
                                }
                                if current_d > d {
                                    d = current_d;
                                    d_idx = *p;
                                }
                            }
                            if j == ampl_for_row[*p].0 {
                                let current_u = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                                if first {
                                    first = false;
                                    u = current_u;
                                    d = u-1;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_u > u {
                                    u = current_u;
                                    u_idx = *p;
                                }
                            }
                            if j == ampl_for_row[*p].1 + 1 {
                                let current_d = m[*p][j-1] + scores_matrix.get(&(sequence[j], graph[i].0)).unwrap();
                                if first {
                                    first = false;
                                    d = current_d;
                                    u = d-1;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_d > d {
                                    d = current_d;
                                    d_idx = *p;
                                }
                            }
                            
                        }
                        if first {
                            panic!("ampl error: {:?}", j);
                        }
                    }
                    match d.cmp(&u) {
                        Ordering::Less => {
                            m[i][j] = u;
                            path[i][j] = ('U',u_idx);
                        }
                        _ => {
                            m[i][j] = d;
                            path[i][j] = ('D', d_idx);
                        }
                    }
                },
                _ if j == right-1 => {
                    //only look left or diagonal
                    let l = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    let mut d= 0;
                    let mut d_idx = 0;
                    if graph[i].1.is_empty() {
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        d_idx = i-1;
                    } else {
                        let mut first = true;
                        for p in graph[i].1.iter(){
                            if j - 1 >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 + 1{
                                let current_d =  m[*p][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
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
                            if j == ampl_for_row[*p].1 + 1{
                                let current_d = m[*p][j-1] + scores_matrix.get(&(sequence[j], graph[i].0)).unwrap();
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
                            panic!("ampl error");
                        }
                    }
                    match d.cmp(&l) {
                        Ordering::Less => {
                            m[i][j] = l;
                            path[i][j] = ('L',j-1);
                        }
                        _ => {
                            m[i][j] = d;
                            path[i][j] = ('D', d_idx);
                        }
                    }
                },
                _ => {
                    // d u or l
                    let l = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    let mut u = 0;
                    let mut d= 0;
                    let mut u_idx = 0;
                    let mut d_idx = 0;

                    if graph[i].1.is_empty() {
                        u = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                        u_idx = i-1;
                        d_idx = i-1;
                    } else {
                        let mut first = true;
                        for p in graph[i].1.iter(){
                            if j - 1 >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 {
                                let current_u = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                                let current_d =  m[*p][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                                if first {
                                    first = false;
                                    u = current_u;
                                    d = current_d;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_u > u {
                                    u = current_u;
                                    u_idx = *p;
                                }
                                if current_d > d {
                                    d = current_d;
                                    d_idx = *p;
                                }
                            }
                            if j == ampl_for_row[*p].0 {
                                let current_u = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                                if first {
                                    first = false;
                                    u = current_u;
                                    d = u-1;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_u > u {
                                    u = current_u;
                                    u_idx = *p;
                                }
                            }
                            if j == ampl_for_row[*p].1 + 1 {
                                let current_d = m[*p][j-1] + scores_matrix.get(&(sequence[j], graph[i].0)).unwrap();
                                if first {
                                    first = false;
                                    d = current_d;
                                    u = d-1;
                                    u_idx = *p;
                                    d_idx = *p;
                                }
                                if current_d > d {
                                    d = current_d;
                                    d_idx = *p;
                                }
                            }
                        }
                        if first {
                            panic!("ampl error");
                        }
                    }
                    m[i][j] = *[d, u, l].iter().max().unwrap();
                    if m[i][j] == d {
                        path[i][j] = ('D', d_idx);
                    } else if m[i][j] == l {
                        path[i][j] = ('L', j-1);
                    } else{
                        path[i][j] = ('U', u_idx);
                    }
                }
            }
        }
    }
    for j in 0..m[0].len() {
        best_last_node(graph, &mut m, &mut path, j);
    }

    match ampl_is_enough(&path) {
        true => {
            /* 
            m.iter().for_each(|line|{println!("{:?}", line)});
            path.iter().for_each(|line|{println!("{:?}", line)});
            */
            println!("{}", m[graph.len()-1][sequence.len()-1]);
        },
        false => {
            println!("{}", m[graph.len()-1][sequence.len()-1]);
            println!("need double ampl");  
        }
    }
    /* 
    println!("Best alignment: {}", m[sequence.len() - 1][graph.len() - 1]);

    basic_output::write_align_poa(&path, sequence, graph);
    */
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

fn set_left_right(ampl:usize, i: usize, graph: &[(char, Vec<usize>)], sequence: &[char]) -> (usize, usize){
    let mut left = 0;
    let mut right = sequence.len();
    if i == 0 {
        right = ampl/2;
    } else {
        if graph[i].1.is_empty() {
            if i < ampl/2 {
                left = 0;
            } else {
                left = i-ampl/2
            }
            right = cmp::min(i+ampl/2, sequence.len());
            if right-ampl/2 <= left {
                (left,right) = (right-ampl/2, right);
            } 
        } else {
            let mut first = true;
            for p in graph[i].1.iter() {
                let (current_l, current_r);
                if p+1 < ampl/2 {
                    current_l = 0;
                } else {
                    current_l = p+1-ampl/2;
                }
                current_r = cmp::min(p+1+ampl/2, sequence.len());
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
            if right-ampl/2 <= left {
                (left, right) = (right-ampl/2, right);
            } 
        }
        
    }
    (left, right)
}

fn ampl_is_enough(path: &[Vec<(char, usize)>]) -> bool {
    let mut col = path[0].len()-1;
    let mut row = path[path.len() - 1][path[0].len()-1].1 as usize;

    let col_number = path[0].len();

    while path[row][col].0 != 'O' {
        if col == 0 || col == col_number - 1 {
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
