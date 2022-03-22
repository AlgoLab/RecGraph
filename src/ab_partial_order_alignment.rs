use std::{cmp::{Ordering, self}, collections::HashMap, panic};

use crate::basic_output;

pub fn exec(
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
    scores_matrix: &HashMap<(char, char), i32>,
) {
    let mut m = vec![vec![0; sequence.len()]; graph.len()];
    let mut ampl_for_row: Vec<(usize, usize)> = vec![(0,0); graph.len()];
    //let mut path = vec![vec![('x', 0); sequence.len()]; graph.len()];
    let ampl = 41; // ampl will be passed as an argument
    for i in 0..graph.len() - 1 {
        let (left, right) = set_left_right(ampl, i, graph, sequence);
        ampl_for_row[i] = (left, right);
        println!("{} {} ampl: {}", left, right, right-left);
        for j in left..right {
            match (i, j) {
                (0,0) => {},
                (_,0) => {
                    if graph[i].1.is_empty() {
                        m[i][j] = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                    } else {
                        let p = graph[i].1.iter().min().unwrap();
                        if j < ampl_for_row[*p].0 || j > ampl_for_row[*p].1 {
                            panic!("ampl wrong");
                        }
                        m[i][j] = m[*p][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                    }
                },
                (0,_) => {
                    m[i][j] = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                },
                _ if j == left => {
                    // only look up or diagonal
                    let mut u = 0;
                    let mut d= 0;
                    if graph[i].1.is_empty() {
                        u = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
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
                                }
                                if current_u > u {
                                    u = current_u;
                                }
                                if current_d > d {
                                    d = current_d;
                                }
                            }
                            
                        }
                        if first {
                            panic!("ampl error");
                        }
                    }
                    m[i][j] = *[d, u].iter().max().unwrap();
                },
                _ if j == right-1 => {
                    //only look left or diagonal
                    let l = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    let mut d= 0;
                    if graph[i].1.is_empty() {
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                    } else {
                        let mut first = true;
                        for p in graph[i].1.iter(){
                            if j - 1 >= ampl_for_row[*p].0 && j <= ampl_for_row[*p].1 + 1{
                                let current_d =  m[*p][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
                                if first {
                                    first = false;
                                    d = current_d;
                                }
                                if current_d > d {
                                    d = current_d;
                                }
                            }
                            
                        }
                        if first {
                            panic!("ampl error");
                        }
                    }
                    m[i][j] = *[d, l].iter().max().unwrap();
                },
                _ => {
                    // d u or l
                    let l = m[i][j-1] + scores_matrix.get(&('-', sequence[j])).unwrap();
                    let mut u = 0;
                    let mut d= 0;
                    if graph[i].1.is_empty() {
                        u = m[i-1][j] + scores_matrix.get(&('-', graph[i].0)).unwrap();
                        d = m[i-1][j-1] + scores_matrix.get(&(graph[i].0, sequence[j])).unwrap();
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
                                }
                                if current_u > u {
                                    u = current_u;
                                }
                                if current_d > d {
                                    d = current_d;
                                }
                            }
                            
                        }
                        if first {
                            panic!("ampl error");
                        }
                    }
                    m[i][j] = *[d, u, l].iter().max().unwrap();
                }
            }
        }
    }
    
    for j in 0..m[0].len() {
        best_last_node(graph, &mut m, j);
    }
    println!("{:?}", m);

    /* 
    println!("Best alignment: {}", m[sequence.len() - 1][graph.len() - 1]);

    basic_output::write_align_poa(&path, sequence, graph);
    */
}

fn best_last_node(
    graph: &[(char, Vec<usize>)],
    m: &mut [Vec<i32>],
    j: usize,
) {
    let mut best_align = 0;
    let mut first = true;
    let m_len = m.len();

    for p in graph[graph.len() - 1].1.iter() {
        if first {
            best_align = m[*p][j];
            first = false;
        }
        if m[*p][j] > best_align {
            best_align = m[*p][j];
        }
    }
    m[m_len - 1][j] = best_align;
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