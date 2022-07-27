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
                (0, 0) => {
                    dpm[i][j] = vec![0; path_number];
                    x[i][j] = vec![0; path_number];
                    y[i][j] = vec![0; path_number];
                }
                (_, 0) => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path == alphas[i] {
                                        x[i][j][path] =
                                            if i == 1 { o + e } else { x[i - 1][j][path] + e };
                                        dpm[i][j][path] = x[i][j][path];
                                    } else {
                                        x[i][j][path] = x[i - 1][j][path];
                                        dpm[i][j][path] = x[i][j][path];
                                    }
                                }
                            }
                        } else {
                            x[i][j][alphas[i]] = if i == 1 {
                                o + e
                            } else {
                                x[i - 1][j][alphas[i]] + x[i - 1][j][alphas[i - 1]] + e
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
                        //TODO: restart here
                    }
                }
                (0, _) => {
                    //first row, only left move
                    y[i][j][alphas[0]] = o + e * j as i32;
                    for k in alphas[0] + 1..path_number {
                        y[i][j][k] = y[i][j - 1][k];
                    }

                    // set dpm
                    dpm[i][j][alphas[0]] = y[i][j][alphas[0]];
                    for k in alphas[0] + 1..path_number {
                        dpm[i][j][k] = y[i][j][k];
                    }
                }
                _ => {}
            }
        }
    }
}
