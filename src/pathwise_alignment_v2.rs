use std::collections::HashMap;

use crate::pathwise_graph::PathGraph;

pub fn exec(sequence: &[char], graph: &PathGraph, score_matrix: &HashMap<(char, char), i32>) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let paths_nodes = &graph.paths_nodes;
    let paths_number = graph.paths_number;

    let mut dpm = vec![vec![vec![0; paths_number]; sequence.len()]; lnz.len()];
    let mut alphas = vec![paths_number + 1; lnz.len()];
    alphas[0] = 0;

    for i in 1..lnz.len() - 1 {
        if !nodes_with_pred[i] {
            //single predecessor
            let mut common_paths = paths_nodes[i].clone();
            common_paths.and(&paths_nodes[i - 1]);

            if common_paths[alphas[i - 1]] {
                // same alpha as predecessor
                alphas[i] = alphas[i - 1];
                for path in 0..paths_number {
                    if common_paths[path] {
                        if path == alphas[i] {
                            dpm[i][0][path] =
                                dpm[i - 1][0][path] + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path];
                        }
                    }
                }
            } else {
                // define new alpha for node i, update score
                let mut first = true;
                for path in 0..paths_number {
                    if common_paths[path] {
                        if first {
                            first = false;
                            alphas[i] = path;
                            dpm[i][0][alphas[i]] = dpm[i - 1][0][alphas[i]]
                                + dpm[i - 1][0][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path] - dpm[i - 1][0][alphas[i]];
                        }
                    }
                }
            }
        } else {
            //multiple predecessor
            //TODO: restart here
        }
    }

    for j in 1..sequence.len() {
        dpm[0][j][alphas[0]] =
            dpm[0][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
        for k in alphas[0] + 1..paths_number {
            dpm[0][j][k] = dpm[0][j - 1][k];
        }
    }

    for i in 1..lnz.len() - 1 {
        for j in 1..sequence.len() {}
    }
}
