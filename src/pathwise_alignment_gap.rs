use crate::pathwise_alignment_output::build_alignment_gap;
use crate::pathwise_graph::PathGraph;
use std::collections::HashMap;
pub fn exec(
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
    o: i32,
    e: i32,
) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let path_number = graph.paths_number;
    let path_node = &graph.paths_nodes;

    let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let mut x = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let mut y = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let alphas = &graph.alphas;

    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {}
                (0, _) => {
                    //set y
                    y[i][j][alphas[0]] = o + e * j as i32;
                    dpm[i][j][alphas[0]] = y[i][j][alphas[0]];
                    for k in alphas[0] + 1..path_number {
                        y[i][j][k] = y[i][j - 1][k];
                        dpm[i][j][k] = y[i][j][k];
                    }
                }
                (_, 0) => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path == alphas[i] {
                                        if i == 1 {
                                            x[i][j][path] = o + e;
                                        } else {
                                            x[i][j][path] = x[i - 1][j][path] + e;
                                        }
                                    } else {
                                        x[i][j][path] = x[i - 1][j][path]
                                    }
                                    dpm[i][j][path] = x[i][j][path];
                                }
                            }
                        } else {
                            x[i][j][alphas[i]] = if i != 1 {
                                x[i - 1][j][alphas[i]] + x[i - 1][j][alphas[i - 1]] + e
                            } else {
                                o + e
                            };
                            dpm[i][j][alphas[i]] = x[i][j][alphas[i]];

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in && path != alphas[i] {
                                    x[i][j][path] = x[i - 1][j][path] - x[i - 1][j][alphas[i]];
                                    dpm[i][j][path] = x[i][j][path];
                                }
                            }
                        }
                    } else {
                        let mut alphas_deltas = HashMap::new();
                        for (p, p_paths) in pred_hash.get_preds_and_paths(i) {
                            let mut common_paths = path_node[i].clone();
                            common_paths.and(&p_paths);

                            if common_paths[alphas[p]] {
                                let paths = common_paths
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(path_id, is_in)| match is_in {
                                        true => Some(path_id),
                                        false => None,
                                    })
                                    .collect::<Vec<usize>>();
                                alphas_deltas.insert(alphas[p], paths);

                                x[i][j][alphas[p]] = if p == 0 {
                                    o + e
                                } else {
                                    x[p][j][alphas[p]] + e
                                };
                                dpm[i][j][alphas[p]] = x[i][j][alphas[p]];

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[p] {
                                        x[i][j][path] = x[p][j][path];
                                        dpm[i][j][path] = x[i][j][path];
                                    }
                                }
                            } else {
                                let temp_alpha = if common_paths[alphas[i]] {
                                    alphas[i]
                                } else {
                                    common_paths.iter().position(|is_in| is_in).unwrap()
                                };
                                let paths = common_paths
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(path_id, is_in)| match is_in {
                                        true => Some(path_id),
                                        false => None,
                                    })
                                    .collect::<Vec<usize>>();
                                alphas_deltas.insert(temp_alpha, paths);

                                x[i][j][temp_alpha] = if p == 0 {
                                    o + e
                                } else {
                                    x[p][j][temp_alpha] + x[p][j][alphas[p]] + e
                                };
                                dpm[i][j][temp_alpha] = x[i][j][temp_alpha];

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != temp_alpha {
                                            x[i][j][path] = x[p][j][path] - x[p][j][temp_alpha];
                                            dpm[i][j][path] = x[i][j][path]
                                        }
                                    }
                                }
                            }
                        }
                        // remove multiple alpha
                        if alphas_deltas.keys().len() > 0 {
                            for (a, delta) in alphas_deltas.iter() {
                                if *a != alphas[i] {
                                    x[i][j][*a] -= x[i][j][alphas[i]];
                                    dpm[i][j][*a] = x[i][j][*a];
                                    for path in delta.iter() {
                                        if path != a {
                                            x[i][j][*path] += x[i][j][*a];
                                            dpm[i][j][*path] = x[i][j][*path]
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }
}
