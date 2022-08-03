use crate::pathwise_graph::PathGraph;
use std::{cmp::Ordering, collections::HashMap};
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

                                x[i][j][alphas[p]] = if i == 1 {
                                    o + e
                                } else {
                                    x[p][j][alphas[p]] + e
                                };
                                dpm[i][j][alphas[p]] = x[i][j][alphas[p]];
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[p] {
                                        x[i][j][path] = x[p][j][path];
                                        dpm[i][j][path] = x[p][j][path];
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

                                x[i][j][temp_alpha] = if i == 1 {
                                    o + e
                                } else {
                                    x[p][j][alphas[p]] + x[p][j][temp_alpha] + e
                                };
                                dpm[i][j][temp_alpha] = x[i][j][temp_alpha];

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != temp_alpha {
                                            x[i][j][path] = x[p][j][path] - x[p][j][temp_alpha];
                                            dpm[i][j][path] = x[i][j][path];
                                        }
                                    }
                                }
                            }
                        }
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
                _ => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            //set x
                            x[i][j][alphas[i]] = match (x[i][j - 1][alphas[i]] + e)
                                .cmp(&(dpm[i][j - 1][alphas[i]] + o + e))
                            {
                                Ordering::Greater => x[i][j - 1][alphas[i]] + e,
                                _ => dpm[i][j - 1][alphas[i]] + o + e,
                            };
                            let l = x[i][j][alphas[i]];

                            //set y
                            y[i][j][alphas[i]] = match (y[i - 1][j][alphas[i - 1]] + e)
                                .cmp(&(dpm[i - 1][j][alphas[i - 1]] + o + e))
                            {
                                Ordering::Greater => y[i - 1][j][alphas[i - 1]] + e,
                                _ => dpm[i - 1][j][alphas[i - 1]] + o + e,
                            };
                            let u = y[i][j][alphas[i]];

                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path != alphas[i] {
                                        if dpm[i][j][alphas[i]] == d {
                                            dpm[i][j][path] = dpm[i - 1][j - 1][path];
                                        } else if dpm[i][j][alphas[i]] == u {
                                            dpm[i][j][path] = dpm[i - 1][j][path];
                                        } else {
                                            dpm[i][j][path] = dpm[i][j - 1][path];
                                        }
                                    }
                                }
                            }
                        } else {
                            //TODO: restart here, alphas i-1 != alphas
                            //set x
                            x[i][j][alphas[i]] = match (x[i][j - 1][alphas[i]] + e)
                                .cmp(&(dpm[i][j - 1][alphas[i]] + o + e))
                            {
                                Ordering::Greater => x[i][j - 1][alphas[i]] + e,
                                _ => dpm[i][j - 1][alphas[i]] + o + e,
                            };
                            let l = x[i][j][alphas[i]];

                            //set y
                            y[i][j][alphas[i]] =
                                match (y[i - 1][j][alphas[i - 1]] + y[i - 1][j][alphas[i]] + e)
                                    .cmp(&(y[i - 1][j][alphas[i - 1]] + y[i - 1][j][alphas[i]] + e))
                                {
                                    Ordering::Greater => {
                                        y[i - 1][j][alphas[i - 1]] + y[i - 1][j][alphas[i]] + e
                                    }
                                    _ => y[i - 1][j][alphas[i - 1]] + y[i - 1][j][alphas[i]] + e,
                                };
                            let u = y[i][j][alphas[i]];

                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + dpm[i - 1][j - 1][alphas[i]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();
                        }
                    } else {
                        //multiple alphas
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

                                let u = dpm[p][j][alphas[p]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                let d = dpm[p][j - 1][alphas[p]]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                                let l = if alphas[i] == alphas[p] {
                                    dpm[i][j - 1][alphas[p]]
                                        + score_matrix.get(&(sequence[j], '-')).unwrap()
                                } else {
                                    dpm[i][j - 1][alphas[p]]
                                        + dpm[i][j - 1][alphas[i]]
                                        + score_matrix.get(&(sequence[j], '-')).unwrap()
                                };
                                dpm[i][j][alphas[p]] = *[d, u, l].iter().max().unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[p] {
                                            if dpm[i][j][alphas[p]] == d {
                                                dpm[i][j][path] = dpm[p][j - 1][path];
                                            } else if dpm[i][j][alphas[p]] == u {
                                                dpm[i][j][path] = dpm[p][j][path];
                                            } else {
                                                dpm[i][j][path] = dpm[i][j - 1][path];
                                            }
                                        }
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

                                let u = dpm[p][j][alphas[p]]
                                    + dpm[p][j][temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                let d = dpm[p][j - 1][alphas[p]]
                                    + dpm[p][j - 1][temp_alpha]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                                let l = if alphas[i] == temp_alpha {
                                    dpm[i][j - 1][temp_alpha]
                                        + score_matrix.get(&(sequence[j], '-')).unwrap()
                                } else {
                                    dpm[i][j - 1][temp_alpha]
                                        + dpm[i][j - 1][alphas[i]]
                                        + score_matrix.get(&(sequence[j], '-')).unwrap()
                                };
                                dpm[i][j][temp_alpha] = *[d, u, l].iter().max().unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if path != temp_alpha {
                                        if is_in {
                                            if dpm[i][j][temp_alpha] == d {
                                                dpm[i][j][path] =
                                                    dpm[p][j - 1][path] - dpm[p][j - 1][temp_alpha];
                                            } else if dpm[i][j][temp_alpha] == u {
                                                dpm[i][j][path] =
                                                    dpm[p][j][path] - dpm[p][j][temp_alpha];
                                            } else {
                                                dpm[i][j][path] = dpm[i][j - 1][path];
                                            }
                                        }
                                    }
                                }
                            }
                        }
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
            }
        }
    }
}