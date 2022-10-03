use bit_vec::BitVec;

use crate::{
    gaf_output::GAFStruct,
    pathwise_graph::{self, PathGraph, PredHash},
    recombination_output,
};
use std::collections::HashMap;
pub fn get_node_offset(nodes_handles: &Vec<u64>, curr_node: usize) -> i32 {
    let handle = nodes_handles[curr_node];
    if handle == 0 {
        0
    } else {
        let mut counter = curr_node;
        let mut offset = 0;
        while nodes_handles[counter - 1] == handle {
            counter -= 1;
            offset += 1;
        }
        offset
    }
}
pub fn exec(
    aln_mode: i32,
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
    base_rec_cost: i32,
    multi_rec_cost: f32,
    displacement_matrix: &Vec<Vec<i32>>,
) -> GAFStruct {
    let forward_matrix = align(aln_mode, sequence, graph, score_matrix);

    let rev_graph = pathwise_graph::create_reverse_path_graph(graph);
    let rev_sequence = get_rev_sequence(&sequence);
    let reverse_matrix = rev_align(aln_mode, &rev_sequence, &rev_graph, score_matrix);

    let (forw_ending_node, rev_starting_node, forw_best_path, rev_best_path, recombination_col) =
        best_alignment(
            &forward_matrix,
            &reverse_matrix,
            &displacement_matrix,
            base_rec_cost,
            multi_rec_cost,
            aln_mode,
            &graph.paths_nodes,
            &graph.pred_hash,
        );
    let gaf;
    if aln_mode == 8 {
        if forw_best_path == rev_best_path {
            gaf = recombination_output::gaf_output_global_no_rec(
                &forward_matrix,
                &graph,
                &sequence,
                &score_matrix,
                forw_best_path,
            );
        } else {
            gaf = recombination_output::gaf_output_global_rec(
                &forward_matrix,
                &reverse_matrix,
                &sequence,
                &score_matrix,
                forw_best_path,
                rev_best_path,
                forw_ending_node,
                rev_starting_node,
                recombination_col,
                &graph.lnz,
                &graph.pred_hash,
                &rev_graph.pred_hash,
                &graph.nwp,
                &rev_graph.nwp,
                &graph.nodes_id_pos,
            );
        }
    } else {
        if forw_best_path == rev_best_path {
            let ending_node = ending_node(&forward_matrix, forw_best_path, &graph.paths_nodes);
            gaf = recombination_output::gaf_output_semiglobal_no_rec(
                &forward_matrix,
                &graph.lnz,
                &sequence,
                &score_matrix,
                forw_best_path,
                &graph.pred_hash,
                &graph.nwp,
                &graph.nodes_id_pos,
                ending_node,
            )
        } else {
            gaf = recombination_output::gaf_output_semiglobal_rec(
                &forward_matrix,
                &reverse_matrix,
                &graph.lnz,
                &sequence,
                &score_matrix,
                forw_best_path,
                rev_best_path,
                &graph.pred_hash,
                &rev_graph.pred_hash,
                &graph.nwp,
                &rev_graph.nwp,
                &graph.nodes_id_pos,
                forw_ending_node,
                rev_starting_node,
                recombination_col,
            );
        }
    }
    gaf
}

fn rev_align(
    aln_mode: i32,
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
) -> Vec<Vec<Vec<i32>>> {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let path_number = graph.paths_number;
    let path_node = &graph.paths_nodes;

    let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let alphas = &graph.alphas;

    let last_node_pos = lnz.len() - 1;
    let last_char_pos = sequence.len() - 1;
    for i in (1..last_node_pos + 1).rev() {
        for j in (1..=last_char_pos).rev() {
            if i == last_node_pos && j == last_char_pos {
                dpm[i][j] = vec![0; path_number];
            } else if i == last_node_pos {
                dpm[i][j][alphas[i]] =
                    dpm[i][j + 1][alphas[i]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                for k in alphas[i] + 1..path_number {
                    dpm[i][j][k] = dpm[i][j + 1][k];
                }
            } else if j == last_char_pos {
                if aln_mode == 9 {
                    dpm[i][j] = vec![0; path_number];
                } else {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i + 1]);

                        if common_paths[alphas[i + 1]] {
                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path == alphas[i] {
                                        dpm[i][j][path] = dpm[i + 1][j][path]
                                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                                    } else {
                                        dpm[i][j][path] = dpm[i + 1][j][path];
                                    }
                                }
                            }
                        } else {
                            dpm[i][j][alphas[i]] = dpm[i + 1][j][alphas[i]]
                                + dpm[i + 1][j][alphas[i + 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in && path != alphas[i] {
                                    dpm[i][j][path] = dpm[i + 1][j][path] - dpm[i + 1][j][alphas[i]]
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

                                dpm[i][j][alphas[p]] = dpm[p][j][alphas[p]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[p] {
                                        dpm[i][j][path] = dpm[p][j][path];
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

                                dpm[i][j][temp_alpha] = dpm[p][j][alphas[p]]
                                    + dpm[p][j][temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in {
                                        if path != temp_alpha {
                                            dpm[i][j][path] =
                                                dpm[p][j][path] - dpm[p][j][temp_alpha];
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
            } else {
                if !nodes_with_pred[i] {
                    let mut common_paths = path_node[i].clone();
                    common_paths.and(&path_node[i + 1]);

                    if common_paths[alphas[i + 1]] {
                        let u = dpm[i + 1][j][alphas[i + 1]]
                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                        let d = dpm[i + 1][j + 1][alphas[i + 1]]
                            + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                        let l = dpm[i][j + 1][alphas[i]]
                            + score_matrix.get(&(sequence[j], '-')).unwrap();

                        dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                        for (path, is_in) in common_paths.iter().enumerate() {
                            if is_in {
                                if path != alphas[i] {
                                    if dpm[i][j][alphas[i]] == d {
                                        dpm[i][j][path] = dpm[i + 1][j + 1][path];
                                    } else if dpm[i][j][alphas[i]] == u {
                                        dpm[i][j][path] = dpm[i + 1][j][path];
                                    } else {
                                        dpm[i][j][path] = dpm[i][j + 1][path];
                                    }
                                }
                            }
                        }
                    } else {
                        let u = dpm[i + 1][j][alphas[i + 1]]
                            + dpm[i + 1][j][alphas[i]]
                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                        let d = dpm[i + 1][j + 1][alphas[i + 1]]
                            + dpm[i + 1][j + 1][alphas[i]]
                            + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                        let l = dpm[i][j + 1][alphas[i]]
                            + score_matrix.get(&(sequence[j], '-')).unwrap();
                        dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                        for (path, is_in) in common_paths.iter().enumerate() {
                            if is_in {
                                if path != alphas[i] {
                                    if dpm[i][j][alphas[i]] == d {
                                        dpm[i][j][path] =
                                            dpm[i + 1][j + 1][path] - dpm[i + 1][j + 1][alphas[i]];
                                    } else if dpm[i][j][alphas[i]] == u {
                                        dpm[i][j][path] =
                                            dpm[i + 1][j][path] - dpm[i + 1][j][alphas[i]];
                                    } else {
                                        dpm[i][j][path] = dpm[i][j + 1][path];
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

                            let u =
                                dpm[p][j][alphas[p]] + score_matrix.get(&(lnz[i], '-')).unwrap();
                            let d = dpm[p][j + 1][alphas[p]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            let l = if alphas[i] == alphas[p] {
                                dpm[i][j + 1][alphas[p]]
                                    + score_matrix.get(&(sequence[j], '-')).unwrap()
                            } else {
                                dpm[i][j + 1][alphas[p]]
                                    + dpm[i][j + 1][alphas[i]]
                                    + score_matrix.get(&(sequence[j], '-')).unwrap()
                            };
                            dpm[i][j][alphas[p]] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path != alphas[p] {
                                        if dpm[i][j][alphas[p]] == d {
                                            dpm[i][j][path] = dpm[p][j + 1][path];
                                        } else if dpm[i][j][alphas[p]] == u {
                                            dpm[i][j][path] = dpm[p][j][path];
                                        } else {
                                            if alphas[p] == alphas[i] {
                                                dpm[i][j][path] = dpm[i][j + 1][path];
                                            } else {
                                                dpm[i][j][path] =
                                                    dpm[i][j + 1][path] - dpm[i][j + 1][alphas[p]];
                                            }
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
                            let d = dpm[p][j + 1][alphas[p]]
                                + dpm[p][j + 1][temp_alpha]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            let l = if alphas[i] == temp_alpha {
                                dpm[i][j + 1][temp_alpha]
                                    + score_matrix.get(&(sequence[j], '-')).unwrap()
                            } else {
                                dpm[i][j + 1][temp_alpha]
                                    + dpm[i][j + 1][alphas[i]]
                                    + score_matrix.get(&(sequence[j], '-')).unwrap()
                            };
                            dpm[i][j][temp_alpha] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if path != temp_alpha {
                                    if is_in {
                                        if dpm[i][j][temp_alpha] == d {
                                            dpm[i][j][path] =
                                                dpm[p][j + 1][path] - dpm[p][j + 1][temp_alpha];
                                        } else if dpm[i][j][temp_alpha] == u {
                                            dpm[i][j][path] =
                                                dpm[p][j][path] - dpm[p][j][temp_alpha];
                                        } else {
                                            if temp_alpha == alphas[i] {
                                                dpm[i][j][path] = dpm[i][j + 1][path];
                                            } else {
                                                dpm[i][j][path] =
                                                    dpm[i][j + 1][path] - dpm[i][j + 1][temp_alpha];
                                            }
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

    absolute_scores(&mut dpm, alphas, path_node);
    dpm
}
fn align(
    aln_mode: i32,
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
) -> Vec<Vec<Vec<i32>>> {
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
                    if aln_mode == 9 {
                        dpm[i][j] = vec![0; path_number]
                    } else {
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
                                        dpm[i][j][path] =
                                            dpm[i - 1][j][path] - dpm[i - 1][j][alphas[i]]
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

                                    dpm[i][j][alphas[p]] = dpm[p][j][alphas[p]]
                                        + score_matrix.get(&(lnz[i], '-')).unwrap();
                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in && path != alphas[p] {
                                            dpm[i][j][path] = dpm[p][j][path];
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

                                    dpm[i][j][temp_alpha] = dpm[p][j][alphas[p]]
                                        + dpm[p][j][temp_alpha]
                                        + score_matrix.get(&(lnz[i], '-')).unwrap();

                                    for (path, is_in) in common_paths.iter().enumerate() {
                                        if is_in {
                                            if path != temp_alpha {
                                                dpm[i][j][path] =
                                                    dpm[p][j][path] - dpm[p][j][temp_alpha];
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
                }
                (0, _) => {
                    dpm[i][j][alphas[0]] =
                        dpm[i][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                    for k in alphas[0] + 1..path_number {
                        dpm[i][j][k] = dpm[i][j - 1][k];
                    }
                }
                _ => {
                    if !nodes_with_pred[i] {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&path_node[i - 1]);

                        if common_paths[alphas[i - 1]] {
                            let u = dpm[i - 1][j][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            let l = dpm[i][j - 1][alphas[i]]
                                + score_matrix.get(&(sequence[j], '-')).unwrap();

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
                            let u = dpm[i - 1][j][alphas[i - 1]]
                                + dpm[i - 1][j][alphas[i]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                            let d = dpm[i - 1][j - 1][alphas[i - 1]]
                                + dpm[i - 1][j - 1][alphas[i]]
                                + score_matrix.get(&(lnz[i], sequence[j])).unwrap();
                            let l = dpm[i][j - 1][alphas[i]]
                                + score_matrix.get(&(sequence[j], '-')).unwrap();
                            dpm[i][j][alphas[i]] = *[d, u, l].iter().max().unwrap();

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in {
                                    if path != alphas[i] {
                                        if dpm[i][j][alphas[i]] == d {
                                            dpm[i][j][path] = dpm[i - 1][j - 1][path]
                                                - dpm[i - 1][j - 1][alphas[i]];
                                        } else if dpm[i][j][alphas[i]] == u {
                                            dpm[i][j][path] =
                                                dpm[i - 1][j][path] - dpm[i - 1][j][alphas[i]];
                                        } else {
                                            dpm[i][j][path] = dpm[i][j - 1][path];
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
                                                if alphas[p] == alphas[i] {
                                                    dpm[i][j][path] = dpm[i][j - 1][path];
                                                } else {
                                                    dpm[i][j][path] = dpm[i][j - 1][path]
                                                        - dpm[i][j - 1][alphas[p]];
                                                }
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
                                                if temp_alpha == alphas[i] {
                                                    dpm[i][j][path] = dpm[i][j - 1][path];
                                                } else {
                                                    dpm[i][j][path] = dpm[i][j - 1][path]
                                                        - dpm[i][j - 1][temp_alpha];
                                                }
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

    absolute_scores(&mut dpm, alphas, path_node);
    dpm
}

fn absolute_scores(dpm: &mut Vec<Vec<Vec<i32>>>, alphas: &Vec<usize>, paths_nodes: &Vec<BitVec>) {
    for i in 0..dpm.len() - 1 {
        for j in 0..dpm[i].len() {
            for path in 0..dpm[i][j].len() {
                if path != alphas[i] && paths_nodes[i][path] {
                    dpm[i][j][path] += dpm[i][j][alphas[i]]
                }
            }
        }
    }
}

fn best_alignment(
    m: &Vec<Vec<Vec<i32>>>,
    w: &Vec<Vec<Vec<i32>>>,
    dms: &Vec<Vec<i32>>,
    brc: i32,
    mrc: f32,
    aln_mode: i32,
    nodes_path: &Vec<BitVec>,
    pred_hash: &PredHash,
) -> (usize, usize, usize, usize, usize) {
    let mut forw_ending_node = 0;
    let mut rev_starting_node = 0;
    let mut recombination_col = 0;

    let mut max = None;
    let mut best_path = None;
    if aln_mode == 8 {
        let preds = pred_hash.get_preds_and_paths(m.len() - 1);
        for (pred, paths) in preds.iter() {
            for (path, is_in) in paths.iter().enumerate() {
                if is_in {
                    if max.is_none() || max.unwrap() < m[*pred][m[*pred].len() - 1][path] {
                        max = Some(m[*pred][m[*pred].len() - 1][path]);
                        best_path = Some(path)
                    }
                }
            }
        }
    } else {
        for i in 0..m.len() - 1 {
            for path in 0..m[i][m[i].len() - 1].len() {
                if nodes_path[i][path] {
                    if max.is_none() || max.unwrap() < m[i][m[i].len() - 1][path] {
                        max = Some(m[i][m[i].len() - 1][path]);
                        best_path = Some(path);
                    }
                }
            }
        }
    }
    let mut curr_best_score = max.unwrap();
    let mut forw_best_path = best_path.unwrap();
    let mut rev_best_path = best_path.unwrap();

    for j in 1..m[0].len() - 1 {
        for i in 1..m.len() - 1 {
            for rev_i in 1..m.len() - 1 {
                let forw_path = m[i][j]
                    .iter()
                    .enumerate()
                    .map(|(path, score)| (score, path))
                    .max()
                    .unwrap()
                    .1;
                let rev_path = w[rev_i][j]
                    .iter()
                    .enumerate()
                    .map(|(path, score)| (score, path))
                    .max()
                    .unwrap()
                    .1;
                if nodes_path[i][forw_path] && nodes_path[rev_i][rev_path] {
                    let penalty = if forw_path == rev_path {
                        0
                    } else {
                        brc + (mrc * dms[i][rev_i] as f32) as i32
                    };
                    if m[i][j][forw_path] + w[rev_i][j][rev_path] - penalty >= curr_best_score {
                        curr_best_score = m[i][j][forw_path] + w[rev_i][j][rev_path];
                        forw_ending_node = i;
                        rev_starting_node = rev_i;
                        forw_best_path = forw_path;
                        rev_best_path = rev_path;
                        recombination_col = j;
                    }
                }
            }
        }
    }
    (
        forw_ending_node,
        rev_starting_node,
        forw_best_path,
        rev_best_path,
        recombination_col,
    )
}

pub fn get_rev_sequence(seq: &[char]) -> Vec<char> {
    let mut rev_seq = Vec::new();
    for c in seq.iter() {
        rev_seq.push(*c)
    }
    rev_seq.push('F');
    rev_seq.remove(0);
    rev_seq
}

fn ending_node(dpm: &Vec<Vec<Vec<i32>>>, best_path: usize, paths_nodes: &Vec<BitVec>) -> usize {
    let mut best_score = None;
    let mut best_node = 0;
    for i in 1..dpm.len() - 1 {
        if paths_nodes[i][best_path] {
            if best_score.is_none() || dpm[i][dpm[0].len() - 1][best_path] > best_score.unwrap() {
                best_score = Some(dpm[i][dpm[0].len() - 1][best_path]);
                best_node = i;
            }
        }
    }
    best_node
}
