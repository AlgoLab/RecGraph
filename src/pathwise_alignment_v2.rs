use crate::pathwise_graph::PathGraph;
use std::collections::HashMap;
pub fn exec(sequence: &[char], graph: &PathGraph, score_matrix: &HashMap<(char, char), i32>) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let path_number = graph.paths_number;
    let path_node = &graph.paths_nodes;

    let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let alphas = &graph.alphas;
    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {
                    dpm[i][j] = vec![0; path_number];
                }
                (_, 0) => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path == alphas[i] {
                                        dpm[i][j][path] = dpm[i - 1][j][path]
                                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                                    } else {
                                        dpm[i][j][path] = dpm[i - 1][j][path];
                                    }
                                }
                            }
                        } else {
                            dpm[i][j][alphas[i]] = dpm[i - 1][j][alphas[i]]
                                + dpm[i - 1][j][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in && path != alphas[i] {
                                    dpm[i][j][path] = dpm[i - 1][j][path] - dpm[i - 1][j][alphas[i]]
                                }
                            }
                        }
                    } else {
                        let mut alphas_deltas = HashMap::new();
                        for p in pred_hash.get(&i).unwrap() {
                            let mut common_paths = path_node[i].clone();
                            common_paths.and(&path_node[*p]);

                            if common_paths[alphas[*p]] {
                                let paths = common_paths
                                    .iter()
                                    .enumerate()
                                    .filter_map(|(path_id, is_in)| match is_in {
                                        true => Some(path_id),
                                        false => None,
                                    })
                                    .collect::<Vec<usize>>();
                                alphas_deltas.insert(alphas[*p], paths);

                                dpm[i][j][alphas[*p]] = dpm[*p][j][alphas[*p]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[*p] {
                                        dpm[i][j][path] = dpm[*p][j][path];
                                    }
                                }
                            } else {
                                //set new alpha
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

                                dpm[i][j][temp_alpha] = dpm[*p][j][alphas[*p]]
                                    + dpm[*p][j][temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != temp_alpha {
                                            dpm[i][j][path] =
                                                dpm[*p][j][path] - dpm[*p][j][temp_alpha];
                                        }
                                    }
                                }
                            }
                        }
                        // remove multiple alpha
                        if alphas_deltas.keys().len() > 0 {
                            for (a, delta) in alphas_deltas.iter() {
                                if *a != alphas[i] {
                                    dpm[i][j][*a] -= dpm[i][j][alphas[i]];
                                    for path in delta.iter() {
                                        if path != a {
                                            dpm[i][j][*path] += dpm[i][j][*a];
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
                }
                _ => {}
            }
        }
    }

    dpm[..dpm.len() - 1]
        .iter()
        .enumerate()
        .for_each(|(_i, line)| println!("{:?}", line[0]));

    println!("SCORE:");
    println!("{:?}", dpm[dpm.len() - 2][0]);
    println!("{}", alphas[alphas.len() - 2]);
}
