use std::collections::HashMap;

use crate::pathwise_graph::PathGraph;

pub fn exec(sequence: &[char], graph: &PathGraph, score_matrix: &HashMap<(char, char), i32>) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let paths_nodes = &graph.paths_nodes;
    let paths_number = graph.paths_number;
    let alphas = &graph.alphas;
    let mut dpm = vec![vec![vec![0; paths_number]; sequence.len()]; lnz.len()];

    for i in 1..lnz.len() - 1 {
        if !nodes_with_pred[i] {
            //single predecessor
            let mut common_paths = paths_nodes[i].clone();
            common_paths.and(&paths_nodes[i - 1]);

            if common_paths[alphas[i - 1]] {
                // same alpha as predecessor
                common_paths.iter().enumerate().for_each(|(path, is_in)| {
                    if is_in {
                        if path == alphas[i] {
                            dpm[i][0][path] =
                                dpm[i - 1][0][path] + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path];
                        }
                    }
                });
            } else {
                // different alpha for node i, update score
                common_paths.iter().enumerate().for_each(|(path, is_in)| {
                    if is_in {
                        if path == alphas[i] {
                            dpm[i][0][alphas[i]] = dpm[i - 1][0][alphas[i]]
                                + dpm[i - 1][0][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path] - dpm[i - 1][0][alphas[i]];
                        }
                    }
                });
            }
        } else {
            //multiple predecessor
            let mut temp_alphas = HashMap::new();
            for pred in pred_hash.get(&i).unwrap() {
                let mut common_paths = paths_nodes[i].clone();
                common_paths.and(&paths_nodes[*pred]);
                if alphas[*pred] == alphas[i] {
                    // same alpha, alpha in common paths
                    common_paths.iter().enumerate().for_each(|(path, is_in)| {
                        if is_in {
                            if path == alphas[i] {
                                dpm[i][0][path] =
                                    dpm[*pred][0][path] + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path];
                            }
                        }
                    });
                } else if common_paths[alphas[i]] {
                    //update alphas[i] val and relative deltas
                    common_paths.iter().enumerate().for_each(|(path, is_in)| {
                        if is_in {
                            if path == alphas[i] {
                                dpm[i][0][path] = dpm[*pred][0][alphas[*pred]]
                                    + dpm[*pred][0][path]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path] - dpm[*pred][0][alphas[i]];
                            }
                        }
                    })
                } else if common_paths[alphas[*pred]] {
                    // as if normal alpha but update temp alphas for later
                    let deltas = common_paths
                        .iter()
                        .enumerate()
                        .filter_map(|(path_id, is_in)| match is_in && path_id != alphas[*pred] {
                            true => Some(path_id),
                            false => None,
                        })
                        .collect::<Vec<usize>>();
                    temp_alphas.insert(alphas[*pred], deltas);
                    common_paths.iter().enumerate().for_each(|(path, is_in)| {
                        if is_in {
                            if path == alphas[*pred] {
                                dpm[i][0][path] = dpm[*pred][0][alphas[*pred]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path];
                            }
                        }
                    })
                } else {
                    //temp_alpha needed, restore before leaving the node
                    let mut first = true;
                    let mut temp_alpha = 0;
                    for path in 0..paths_number {
                        if common_paths[path] {
                            if first {
                                first = false;
                                temp_alpha = path;
                                let deltas = common_paths
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(path_id, is_in)| {
                                        match is_in && path_id != temp_alpha {
                                            true => Some(path_id),
                                            false => None,
                                        }
                                    })
                                    .collect::<Vec<usize>>();
                                temp_alphas.insert(temp_alpha, deltas);
                                //update temp alpha score
                                dpm[i][0][temp_alpha] = dpm[*pred][0][alphas[*pred]]
                                    + dpm[*pred][0][temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path] - dpm[*pred][0][temp_alpha];
                            }
                        }
                    }
                }
            }
            // remove temp alpha if needed

            if temp_alphas.keys().len() > 0 {
                for (alpha, deltas) in temp_alphas.iter() {
                    dpm[i][0][*alpha] -= dpm[i][0][alphas[i]];
                    for path in deltas.iter() {
                        dpm[i][0][*path] += dpm[i][0][*alpha];
                    }
                }
            }
        }
    }

    for j in 1..sequence.len() {
        dpm[0][j][alphas[0]] =
            dpm[0][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
        for k in 0..paths_number {
            if k != alphas[0] {
                dpm[0][j][k] = dpm[0][j - 1][k];
            }
        }
    }

    for i in 1..lnz.len() - 1 {
        for j in 1..sequence.len() {}
    }
    dpm.iter().for_each(|line| println!("{:?}", line[0]));
    println!("SCORE");
    println!("{}", dpm[lnz.len() - 2][0][alphas[lnz.len() - 2]]);
    for path in 0..paths_number {
        if path != alphas[lnz.len() - 2] {
            println!(
                "{}",
                dpm[lnz.len() - 2][0][alphas[lnz.len() - 2]] + dpm[lnz.len() - 2][0][alphas[path]]
            );
        }
    }
}
