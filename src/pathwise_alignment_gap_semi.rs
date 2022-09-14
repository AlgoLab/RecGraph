use crate::pathwise_alignment_output::build_alignment_semiglobal_gap;
use crate::pathwise_graph::PathGraph;
use std::collections::HashMap;

pub fn exec(
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
    o: i32,
    e: i32,
) -> usize {
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
                }
                (_, 0) => dpm[0][0] = vec![0; path_number],
                (0, _) => {
                    //set y
                    y[i][j][alphas[0]] = o + e * j as i32;
                    dpm[i][j][alphas[0]] = y[i][j][alphas[0]];
                    for k in alphas[0] + 1..path_number {
                        y[i][j][k] = y[i][j - 1][k];
                        dpm[i][j][k] = y[i][j][k];
                    }
                }
                _ => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            //set y
                            let u_y = y[i - 1][j][alphas[i - 1]] + e;
                            let u_dpm = dpm[i - 1][j][alphas[i - 1]] + o + e;

                            y[i][j][alphas[i]] = if u_dpm >= u_y {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            y[i][j][path] = dpm[i - 1][j][path];
                                        }
                                    }
                                }
                                u_dpm
                            } else {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            y[i][j][path] = y[i - 1][j][path];
                                        }
                                    }
                                }
                                u_y
                            };

                            let u = y[i][j][alphas[i]];

                            //set x
                            let l_x = x[i][j - 1][alphas[i]] + e;
                            let l_dpm = dpm[i][j - 1][alphas[i]] + o + e;

                            x[i][j][alphas[i]] = if l_dpm >= l_x {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            x[i][j][path] = dpm[i][j - 1][path];
                                        }
                                    }
                                }
                                l_dpm
                            } else {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            x[i][j][path] = x[i][j - 1][path];
                                        }
                                    }
                                }
                                l_x
                            };

                            let l = x[i][j][alphas[i]];

                            //set dpm
                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();

                            dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path != alphas[i] {
                                        if dpm[i][j][alphas[i]] == d {
                                            dpm[i][j][path] = dpm[i - 1][j - 1][path];
                                        } else if dpm[i][j][alphas[i]] == u {
                                            dpm[i][j][path] = y[i][j][path];
                                        } else {
                                            dpm[i][j][path] = x[i][j][path];
                                        }
                                    }
                                }
                            }
                        } else {
                            //set y
                            let u_y = y[i - 1][j][alphas[i - 1]] + y[i - 1][j][alphas[i]] + e;
                            let u_dpm =
                                dpm[i - 1][j][alphas[i - 1]] + dpm[i - 1][j][alphas[i]] + o + e;

                            y[i][j][alphas[i]] = if u_dpm >= u_y {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            y[i][j][path] =
                                                dpm[i - 1][j][path] - dpm[i - 1][j][alphas[i]];
                                        }
                                    }
                                }
                                u_dpm
                            } else {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            y[i][j][path] =
                                                y[i - 1][j][path] - y[i - 1][j][alphas[i]];
                                        }
                                    }
                                }
                                u_y
                            };

                            let u = y[i][j][alphas[i]];

                            //set x
                            let l_x = x[i][j - 1][alphas[i]] + e;
                            let l_dpm = dpm[i][j - 1][alphas[i]] + o + e;

                            x[i][j][alphas[i]] = if l_dpm >= l_x {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            x[i][j][path] = dpm[i][j - 1][path];
                                        }
                                    }
                                }
                                l_dpm
                            } else {
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[i] {
                                            x[i][j][path] = x[i][j - 1][path];
                                        }
                                    }
                                }
                                l_x
                            };

                            let l = x[i][j][alphas[i]];

                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + dpm[i - 1][j - 1][alphas[i]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();

                            dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path != alphas[i] {
                                        if dpm[i][j][alphas[i]] == d {
                                            dpm[i][j][path] = dpm[i - 1][j - 1][path]
                                                - dpm[i - 1][j - 1][alphas[i]];
                                        } else if dpm[i][j][alphas[i]] == u {
                                            dpm[i][j][path] = y[i][j][path]
                                        } else {
                                            dpm[i][j][path] = x[i][j][path];
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        // multiple alphas possible
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

                                //set y
                                let u_y = y[p][j][alphas[p]] + e;
                                let u_dpm = dpm[p][j][alphas[p]] + o + e;

                                y[i][j][alphas[p]] = if u_dpm >= u_y {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != alphas[p] {
                                                y[i][j][path] = dpm[p][j][path];
                                            }
                                        }
                                    }
                                    u_dpm
                                } else {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != alphas[i] {
                                                y[i][j][path] = y[p][j][path];
                                            }
                                        }
                                    }
                                    u_y
                                };

                                let u = y[i][j][alphas[p]];

                                //set x
                                let l_x = if alphas[p] == alphas[i] {
                                    x[i][j - 1][alphas[p]] + e
                                } else {
                                    x[i][j - 1][alphas[p]] + x[i][j - 1][alphas[i]] + e
                                };
                                let l_dpm = if alphas[p] == alphas[i] {
                                    dpm[i][j - 1][alphas[p]] + o + e
                                } else {
                                    dpm[i][j - 1][alphas[i]] + dpm[i][j - 1][alphas[p]] + o + e
                                };

                                x[i][j][alphas[p]] = if l_dpm >= l_x {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != alphas[p] {
                                                if alphas[p] == alphas[i] {
                                                    x[i][j][path] = dpm[i][j - 1][path];
                                                } else {
                                                    x[i][j][path] = dpm[i][j - 1][path]
                                                        - dpm[i][j - 1][alphas[p]];
                                                }
                                            }
                                        }
                                    }
                                    l_dpm
                                } else {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != alphas[p] {
                                                if alphas[p] == alphas[i] {
                                                    x[i][j][path] = x[i][j - 1][path];
                                                } else {
                                                    x[i][j][path] =
                                                        x[i][j - 1][path] - x[i][j - 1][alphas[p]];
                                                }
                                            }
                                        }
                                    }
                                    l_x
                                };

                                let l = x[i][j][alphas[p]];

                                //set dpm
                                let d = dpm[p][j - 1][alphas[p]]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();

                                dpm[i][j][alphas[p]] = *[d, u, l].iter().max().unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != alphas[p] {
                                            if dpm[i][j][alphas[p]] == d {
                                                dpm[i][j][path] = dpm[p][j - 1][path];
                                            } else if dpm[i][j][alphas[p]] == u {
                                                dpm[i][j][path] = y[i][j][path];
                                            } else {
                                                dpm[i][j][path] = x[i][j][path];
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

                                //set y
                                let u_y = y[p][j][alphas[p]] + y[p][j][temp_alpha] + e;
                                let u_dpm = dpm[p][j][alphas[p]] + dpm[p][j][temp_alpha] + o + e;

                                y[i][j][temp_alpha] = if u_dpm >= u_y {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != temp_alpha {
                                                y[i][j][path] =
                                                    dpm[p][j][path] - dpm[p][j][temp_alpha];
                                            }
                                        }
                                    }
                                    u_dpm
                                } else {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != temp_alpha {
                                                y[i][j][path] = y[p][j][path] - y[p][j][temp_alpha];
                                            }
                                        }
                                    }
                                    u_y
                                };

                                let u = y[i][j][temp_alpha];

                                //set x
                                let l_x = if alphas[i] == temp_alpha {
                                    x[i][j - 1][alphas[i]] + e
                                } else {
                                    x[i][j - 1][alphas[i]] + x[i][j - 1][temp_alpha] + e
                                };
                                let l_dpm = if alphas[i] == temp_alpha {
                                    dpm[i][j - 1][alphas[i]] + o + e
                                } else {
                                    dpm[i][j - 1][alphas[i]] + dpm[i][j - 1][temp_alpha] + o + e
                                };

                                x[i][j][temp_alpha] = if l_dpm >= l_x {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != temp_alpha {
                                                if temp_alpha == alphas[i] {
                                                    x[i][j][path] = dpm[i][j - 1][path]
                                                } else {
                                                    x[i][j][path] = dpm[i][j - 1][path]
                                                        - dpm[i][j - 1][temp_alpha];
                                                }
                                            }
                                        }
                                    }
                                    l_dpm
                                } else {
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != temp_alpha {
                                                if temp_alpha == alphas[i] {
                                                    x[i][j][path] = x[i][j - 1][path]
                                                } else {
                                                    x[i][j][path] =
                                                        x[i][j - 1][path] - x[i][j - 1][temp_alpha];
                                                }
                                            }
                                        }
                                    }
                                    l_x
                                };

                                let l = x[i][j][temp_alpha];

                                let d = dpm[p][j - 1][alphas[p]]
                                    + dpm[p][j - 1][temp_alpha]
                                    + score_matrix.get(&(lnz[i], sequence[j])).unwrap();

                                dpm[i][j][temp_alpha] = *[d, u, l].iter().max().unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if path != temp_alpha {
                                        if is_in {
                                            if dpm[i][j][temp_alpha] == d {
                                                dpm[i][j][path] =
                                                    dpm[p][j - 1][path] - dpm[p][j - 1][temp_alpha];
                                            } else if dpm[i][j][temp_alpha] == u {
                                                dpm[i][j][path] = y[i][j][path];
                                            } else {
                                                dpm[i][j][path] = x[i][j][path];
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
                                    x[i][j][*a] -= x[i][j][alphas[i]];
                                    y[i][j][*a] -= y[i][j][alphas[i]];

                                    for path in delta.iter() {
                                        if path != a {
                                            dpm[i][j][*path] += dpm[i][j][*a];
                                            x[i][j][*path] += x[i][j][*a];
                                            y[i][j][*path] += y[i][j][*a];
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
    let (final_node, best_path) = best_ending_node(&dpm, graph);
    let cigar_output = build_alignment_semiglobal_gap(
        &dpm,
        &x,
        &y,
        &alphas,
        best_path,
        &pred_hash,
        &nodes_with_pred,
        final_node,
    );
    println!("{}", cigar_output);
    best_path
}

fn best_ending_node(dpm: &Vec<Vec<Vec<i32>>>, graph: &PathGraph) -> (usize, usize) {
    let mut max: Option<i32> = None;
    let mut ending_node: usize = 0;
    let mut chosen_path: usize = 0;
    for i in 0..dpm.len() - 1 {
        let paths = graph.paths_nodes[i].clone();
        let mut absolute_scores = dpm[i][dpm[i].len() - 1].clone();
        for (path, is_in) in paths.iter().enumerate() {
            if is_in {
                if path != graph.alphas[i] {
                    absolute_scores[path] =
                        absolute_scores[path] + absolute_scores[graph.alphas[i]];
                }
            }
        }
        let best_path = absolute_scores
            .iter()
            .enumerate()
            .map(|(path, score)| (score, path))
            .max();
        if max.is_none() || best_path.unwrap().0 > &max.unwrap() {
            max = Some(*best_path.unwrap().0);
            ending_node = i;
            chosen_path = best_path.unwrap().1;
        }
    }
    (ending_node, chosen_path)
}
