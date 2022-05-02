use crate::{bitfield_path as bf, graph::LnzGraph};
use bitvec::prelude::*;
use std::collections::{HashMap, HashSet};
// TODO: edge must be in at least 1 path to be considered
pub fn exec(
    sequence: &[char],
    graph: &LnzGraph,
    path_node: &[Vec<usize>],
    score_matrix: &HashMap<(char, char), i32>,
    path_number: usize,
) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;

    let mut dpm = vec![vec![vec![-10000; path_number]; sequence.len()]; lnz.len()];
    let mut path = vec![vec![bitvec![u16, Msb0; 0; 32]; sequence.len()]; lnz.len()];

    let mut alphas = vec![path_number + 1; lnz.len()];
    
    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {
                    alphas[0] = 0;
                    dpm[i][j] = vec![0; path_number];
                    path[i][j] = bf::set_path_cell(0, 'O');
                }
                (_, 0) => {
                    if !nodes_with_pred[i] {
                        let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                        let pred_paths = path_node[i-1].iter().collect::<HashSet<_>>();
                        let x = curr_node_paths
                            .intersection(&pred_paths)
                            .collect::<Vec<_>>();
                        if x.contains(&&&alphas[i-1]) {
                            alphas[i] = alphas[i-1];
                            for path in x.iter() {
                                if ***path == alphas[i] {
                                    dpm[i][j][***path] = dpm[i-1][j][***path] + score_matrix.get(&(lnz[i], '-')).unwrap();
                                } else {
                                    dpm[i][j][***path] = dpm[i-1][j][***path];
                                }
                            }

                        } else {
                            alphas[i] = ***x.iter().min().unwrap();
                            dpm[i][j][alphas[i]] = dpm[i-1][j][alphas[i]] + dpm[i-1][j][alphas[i-1]] + score_matrix.get(&(lnz[i], '-')).unwrap();
                            for path in x.iter() {
                                if ***path != alphas[i] {
                                    dpm[i][j][***path] = dpm[i-1][j][***path] - dpm[i-1][j][alphas[i]]
                                }
                            }
                        }
                    } else {
                        let mut alphas_deltas = HashMap::new();
                        for p in pred_hash.get(&i).unwrap() {
                            let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                            let pred_paths = path_node[*p].iter().collect::<HashSet<_>>();
                            let x = curr_node_paths
                                .intersection(&pred_paths)
                                .collect::<Vec<_>>();
                            if alphas[i] == path_number + 1 {
                                if x.contains(&&&alphas[*p]) {
                                    alphas[i] = alphas[*p];
                                } else {
                                    alphas[i] = ***x.iter().min().unwrap();
                                }
                            }
                            if x.contains(&&&alphas[*p]) {
                                let paths = x.iter().map(|p| ***p).collect::<Vec<usize>>();
                                alphas_deltas.insert(alphas[*p], paths);
                                
                                dpm[i][j][alphas[*p]] = dpm[*p][j][alphas[*p]] + score_matrix.get(&(lnz[i], '-')).unwrap();
                                for path in x.iter() {
                                    if path != &&&alphas[*p] {
                                        dpm[i][j][***path] = dpm[*p][j][***path];
                                    }
                                }
                            } else {
                                //set new alpha
                                let temp_alpha = x.iter().min().unwrap();
                                let paths = x.iter().map(|p| ***p).collect::<Vec<usize>>();
                                alphas_deltas.insert(***temp_alpha, paths);
                                
                                dpm[i][j][***temp_alpha] = dpm[*p][j][alphas[*p]] + dpm[*p][j][***temp_alpha] + score_matrix.get(&(lnz[i], '-')).unwrap();
                                for path in x.iter() {
                                    if path != temp_alpha {
                                        dpm[i][j][***path] = dpm[*p][j][***path] - dpm[*p][j][***temp_alpha];
                                    }
                                }
                            }
                        }
                        // remove multiple alpha
                        if alphas_deltas.keys().len() > 1 {
                            for (a, delta) in alphas_deltas.iter() {
                                if *a != alphas[i] {
                                    dpm[i][j][*a] -= dpm[i][j][alphas[i]];
                                    for path in delta.iter() {
                                        if path != a {
                                            dpm[i][j][*path] -= dpm[i][j][*a];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                (0, _) => {
                    dpm[i][j][alphas[0]] =
                        dpm[i][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                    for k in alphas[0] + 1..path_number {
                        dpm[i][j][k] = dpm[i][j - 1][k];
                    }
                    path[i][j] = bf::set_path_cell(i, 'L');
                }
                _ => {
                    let l =
                        dpm[i][j - 1][alphas[i]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                    if !nodes_with_pred[i] {
                        let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                        let pred_paths = path_node[i - 1].iter().collect::<HashSet<_>>();
                        let x = curr_node_paths
                            .intersection(&pred_paths)
                            .collect::<Vec<_>>();

                        let u;
                        let d;
                        let mut d_delta_correction = 0;
                        let mut u_delta_correction = 0;
                        if x.contains(&&&alphas[i - 1]) {
                            u = dpm[i - 1][j][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                            d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                        } else {
                            u = dpm[i - 1][j][alphas[i]]
                                + dpm[i - 1][j][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                            d = dpm[i - 1][j - 1][alphas[i]]
                                + dpm[i - 1][j - 1][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            d_delta_correction = dpm[i - 1][j - 1][alphas[i]];
                            u_delta_correction = dpm[i - 1][j][alphas[i]];
                        }
                        dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();
                        if dpm[i][j][alphas[i]] == l {
                            for path in curr_node_paths {
                                if *path != alphas[i] {
                                    dpm[i][j][*path] = dpm[i][j - 1][*path];
                                }
                            }
                        } else {
                            for delta in x.iter() {
                                if dpm[i][j][alphas[i]] == d {
                                    dpm[i][j][***delta] =
                                        dpm[i - 1][j - 1][***delta] - d_delta_correction;
                                } else if dpm[i][j][alphas[i]] == u {
                                    dpm[i][j][***delta] =
                                        dpm[i - 1][j][***delta] - u_delta_correction;
                                }
                            }
                        }
                    } else {
                        //consider multiple alfa
                        let mut alphas_deltas = HashMap::new();
                        for p in pred_hash.get(&i).unwrap() {
                            let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                            let pred_paths = path_node[*p].iter().collect::<HashSet<_>>();
                            let x = curr_node_paths
                                .intersection(&pred_paths)
                                .collect::<Vec<_>>();
                            let u;
                            let d;
                            if x.contains(&&&alphas[*p]) {
                                let paths = x.iter().map(|p| ***p).collect::<Vec<usize>>();
                                alphas_deltas.insert(alphas[*p], paths);
                                u = dpm[*p][j][alphas[*p]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                d = dpm[*p][j - 1][alphas[*p]]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                                dpm[i][j][alphas[*p]] = *[d, u, l].iter().max().unwrap();

                                for path in x.iter() {
                                    if ***path != alphas[*p] {
                                        if dpm[i][j][alphas[*p]] == d {
                                            dpm[i][j][***path] = dpm[*p][j - 1][***path];
                                        } else if dpm[i][j][alphas[*p]] == u {
                                            dpm[i][j][***path] = dpm[*p][j][***path];
                                        } else {
                                            dpm[i][j][***path] = dpm[i][j - 1][***path];
                                        }
                                    }
                                }
                            } else {
                                //set new alpha
                                let temp_alpha = x.iter().min().unwrap();
                                let paths = x.iter().map(|p| ***p).collect::<Vec<usize>>();
                                alphas_deltas.insert(***temp_alpha, paths);
                                u = dpm[*p][j][alphas[*p]]
                                    + dpm[*p][j][***temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                d = dpm[*p][j - 1][alphas[*p]]
                                    + dpm[*p][j - 1][***temp_alpha]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                                dpm[i][j][***temp_alpha] = *[d, u, l].iter().max().unwrap();
                                for path in x.iter() {
                                    if path != temp_alpha {
                                        if dpm[i][j][***temp_alpha] == d {
                                            dpm[i][j][***path] = dpm[*p][j - 1][***path]
                                                - dpm[*p][j - 1][***temp_alpha];
                                        } else if dpm[i][j][***temp_alpha] == u {
                                            dpm[i][j][***path] =
                                                dpm[*p][j][***path] - dpm[*p][j][***temp_alpha];
                                        } else {
                                            dpm[i][j][***path] = dpm[i][j - 1][***path]
                                                - dpm[i][j - 1][***temp_alpha];
                                        }
                                    }
                                }
                            }
                            // if multiple alpha in i, j remove all but alphas[i]
                            if alphas_deltas.keys().len() > 1 {
                                for (a, delta) in alphas_deltas.iter() {
                                    if *a != alphas[i] {
                                        dpm[i][j][*a] -= dpm[i][j][alphas[i]];
                                        for path in delta.iter() {
                                            if path != a {
                                                dpm[i][j][*path] -= dpm[i][j][*a];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    dpm[..dpm.len() - 1]
        .iter()
        .for_each(|l| println!("{:?}", l[l.len()-1]));
    println!("{:?}", alphas);
    println!("{:?}", path_node); 
     
}
