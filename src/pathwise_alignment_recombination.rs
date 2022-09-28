use bit_vec::BitVec;

use crate::{
    gaf_output::GAFStruct,
    pathwise_alignment_output::{self, build_cigar},
    pathwise_graph::{self, PathGraph, PredHash},
};
use std::collections::HashMap;
fn get_node_offset(nodes_handles: &Vec<u64>, curr_node: usize) -> i32 {
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
) -> GAFStruct {
    let forward_matrix = align(aln_mode, sequence, graph, score_matrix);

    let rev_graph = pathwise_graph::create_reverse_path_graph(graph);

    let reverse_matrix = rev_align(aln_mode, sequence, &rev_graph, score_matrix);

    let displacement_matrix = pathwise_graph::nodes_displacement_matrix(&graph.paths_nodes);

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
    let mut gaf = GAFStruct::new();
    if aln_mode == 8 {
        if forw_best_path == rev_best_path {
            gaf = pathwise_alignment_output::build_alignment(
                &forward_matrix,
                &graph,
                &sequence,
                &score_matrix,
                forw_best_path,
            );
        } else {
            //TODO: set gaf for global recombination
            let fen_handle = graph.nodes_id_pos[forw_ending_node];
            let fen_offset = get_node_offset(&graph.nodes_id_pos, forw_ending_node);

            let rsn_handle = graph.nodes_id_pos[rev_starting_node];
            let rsn_offset = get_node_offset(&graph.nodes_id_pos, rev_starting_node);

            println!("Recombination between paths {forw_best_path} and {rev_best_path}");
            println!("Recombination edge between node {fen_handle} (position: {fen_offset}) and node {rsn_handle} (position: {rsn_offset})");
            println!()
        }
    } else {
        if forw_best_path == rev_best_path {
            let ending_node = ending_node(&forward_matrix, forw_best_path, &graph.paths_nodes);
            gaf = gaf_output_semiglobal_no_rec(
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
            gaf = gaf_output_semiglobal_rec(
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
        for j in (1..last_char_pos + 1).rev() {
            if i == last_node_pos && j == last_char_pos {
                dpm[i][j] = vec![0; path_number];
            } else if i == last_node_pos {
                dpm[i][j][alphas[0]] =
                    dpm[i][j + 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                for k in alphas[0] + 1..path_number {
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
                                            dpm[i][j][path] = dpm[i][j + 1][path];
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
                                            dpm[i][j][path] = dpm[i][j + 1][path];
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
    let mut results = vec![0; path_number];
    for (pred, paths) in pred_hash.get_preds_and_paths(0) {
        for (path, is_in) in paths.iter().enumerate() {
            if is_in {
                if path == alphas[pred] {
                    results[path] = dpm[pred][dpm[pred].len() - 1][path]
                } else {
                    results[path] = dpm[pred][dpm[pred].len() - 1][path]
                        + dpm[pred][dpm[pred].len() - 1][alphas[pred]]
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

    for j in 0..m[0].len() - 1 {
        for i in 0..m.len() - 1 {
            for rev_i in 0..m.len() - 1 {
                for forw_path in 0..m[0][0].len() {
                    if nodes_path[i][forw_path] {
                        for rev_path in 0..m[0][0].len() {
                            if nodes_path[rev_i][rev_path] {
                                let penalty = if forw_path == rev_path {
                                    0
                                } else {
                                    brc + (mrc * dms[i][rev_i] as f32) as i32
                                };
                                if m[i][j][forw_path] + w[rev_i][j + 1][rev_path] - penalty
                                    > curr_best_score
                                {
                                    curr_best_score =
                                        m[i][j][forw_path] + w[rev_i][j + 1][rev_path];
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

fn gaf_output_semiglobal_rec(
    dpm: &Vec<Vec<Vec<i32>>>,
    rev_dpm: &Vec<Vec<Vec<i32>>>,
    lnz: &Vec<char>,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    best_path: usize,
    rev_best_path: usize,
    pred_hash: &PredHash,
    rev_pred_hash: &PredHash,
    nwp: &BitVec,
    rev_nwp: &BitVec,
    handles_nodes_id: &Vec<u64>,
    forward_ending_node: usize,
    reverse_starting_node: usize,

    rec_col: usize,
) -> GAFStruct {
    let mut cigar = Vec::new();
    let mut path_length: usize = 0;
    let mut i = reverse_starting_node;
    let mut j = rec_col + 1;
    let mut handle_id_alignment = Vec::new();
    // reverse alignment
    while i > 0 && i < dpm.len() - 1 && j < dpm[0].len() - 1 {
        let curr_score = rev_dpm[i][j][rev_best_path];

        let mut predecessor = None;
        let (d, u, l) = if !rev_nwp[i] {
            (
                rev_dpm[i + 1][j + 1][rev_best_path] + scores.get(&(lnz[i], seq[j])).unwrap(),
                rev_dpm[i + 1][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap(),
                rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', seq[j])).unwrap(),
            )
        } else {
            let preds = rev_pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[rev_best_path] {
                    predecessor = Some(*pred);
                    d = rev_dpm[*pred][j + 1][rev_best_path]
                        + scores.get(&(lnz[i], seq[j])).unwrap();
                    u = rev_dpm[*pred][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    l = rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', seq[j])).unwrap();
                }
            }
            (d, u, l)
        };

        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            j += 1;
            path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        } else {
            cigar.push('L');
            j += 1;
        }
    }
    while j < dpm[0].len() - 1 {
        cigar.push('L');
        j += 1;
    }

    let rev_ending_node = i;
    let mut temp_cigar = Vec::new();
    let mut temp_handle_id_alignment = Vec::new();
    i = forward_ending_node;
    j = rec_col;
    if lnz[i] == seq[j] {
        cigar.push('D')
    } else {
        cigar.push('d')
    }
    while i > 0 && j > 0 {
        let curr_score = dpm[i][j][best_path];
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                dpm[i - 1][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap(),
                dpm[i - 1][j][best_path] + scores.get(&(lnz[i], '-')).unwrap(),
                dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap(),
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    d = dpm[*pred][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap();
                    u = dpm[*pred][j][best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    l = dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap();
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                temp_cigar.push('d');
            } else {
                temp_cigar.push('D');
            }
            temp_handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
            path_length += 1;
        } else if max == u {
            temp_cigar.push('U');
            temp_handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        } else {
            temp_cigar.push('L');
            j -= 1;
        }
    }
    while j > 0 {
        temp_cigar.push('L');
        j -= 1;
    }

    temp_cigar.reverse();
    temp_cigar.append(&mut cigar);

    temp_handle_id_alignment.reverse();
    temp_handle_id_alignment.append(&mut handle_id_alignment);
    temp_handle_id_alignment.dedup();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    let path: Vec<usize> = temp_handle_id_alignment
        .iter()
        .map(|id| *id as usize)
        .collect();
    let path_start = get_node_offset(handles_nodes_id, if i == 0 { i } else { i + 1 }) as usize; // first letter used in first node of alignment
    let path_end = get_node_offset(handles_nodes_id, rev_ending_node) as usize;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let recombination = if best_path == rev_best_path {
        format!("No recombination, best path: {}", best_path)
    } else {
        format!(
            "recombination path {} {}, nodes {} {}",
            best_path,
            rev_best_path,
            handles_nodes_id[forward_ending_node],
            handles_nodes_id[reverse_starting_node]
        )
    };
    let comments = format!("{}, {}", build_cigar(&temp_cigar), recombination);

    let gaf_output = GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_length,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}

fn ending_node(dpm: &Vec<Vec<Vec<i32>>>, best_path: usize, paths_nodes: &Vec<BitVec>) -> usize {
    let mut best_score = None;
    let mut best_node = 0;
    for i in 0..dpm.len() - 1 {
        if paths_nodes[i][best_path] {
            if best_score.is_none() || dpm[i][dpm[0].len() - 1][best_path] > best_score.unwrap() {
                best_score = Some(dpm[i][dpm[0].len() - 1][best_path]);
                best_node = i;
            }
        }
    }
    best_node
}

fn gaf_output_semiglobal_no_rec(
    dpm: &Vec<Vec<Vec<i32>>>,
    lnz: &Vec<char>,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    best_path: usize,
    pred_hash: &PredHash,
    nwp: &BitVec,
    handles_nodes_id: &Vec<u64>,
    ending_node: usize,
) -> GAFStruct {
    let mut cigar = Vec::new();
    let mut path_length: usize = 0;
    let mut i = ending_node;
    let mut j = dpm[i].len() - 1;
    let mut handle_id_alignment = Vec::new();

    while i > 0 && j > 0 {
        let curr_score = dpm[i][j][best_path];
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                dpm[i - 1][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap(),
                dpm[i - 1][j][best_path] + scores.get(&(lnz[i], '-')).unwrap(),
                dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap(),
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    d = dpm[*pred][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap();
                    u = dpm[*pred][j][best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    l = dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap();
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
            path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        } else {
            cigar.push('L');
            j -= 1;
        }
    }
    while j > 0 {
        cigar.push('L');
        j -= 1;
    }

    cigar.reverse();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    let path: Vec<usize> = handle_id_alignment.iter().map(|id| *id as usize).collect();
    //path length already set
    let path_start = get_node_offset(handles_nodes_id, if i == 0 { i } else { i + 1 }) as usize; // first letter used in first node of alignment
    let path_end = get_node_offset(handles_nodes_id, ending_node) as usize; // last letter used in last node of alignment

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = format!("{}, best path: {}", build_cigar(&cigar), best_path);
    let gaf_output = GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_length,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}
