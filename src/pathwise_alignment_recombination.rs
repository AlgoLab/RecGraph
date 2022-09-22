use bit_vec::BitVec;

use crate::pathwise_graph::{self, PathGraph, PredHash};
use std::collections::HashMap;

pub fn exec(
    aln_mode: i32,
    sequence: &[char],
    graph: &PathGraph,
    score_matrix: &HashMap<(char, char), i32>,
    base_rec_cost: i32,
    multi_rec_cost: f32,
) {
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
        );

    if aln_mode == 8 {
        if forw_best_path == rev_best_path {
            println!("No recombination, best path: {forw_best_path}")
        } else {
            let fen_handle = graph.nodes_id_pos[forw_ending_node];
            let mut counter = forw_ending_node;
            let mut fen_offset = 0;
            while graph.nodes_id_pos[counter] == fen_handle {
                counter -= 1;
                fen_offset += 1;
            }

            let rsn_handle = graph.nodes_id_pos[rev_starting_node];
            let mut counter = rev_starting_node;
            let mut rsn_offset = 0;
            while graph.nodes_id_pos[counter] == rsn_handle {
                counter -= 1;
                rsn_offset += 1;
            }

            println!("Recombination between paths {forw_best_path} and {rev_best_path}");
            println!("Recombination edge between node {fen_handle} (position: {fen_offset}) and node {rsn_handle} (position: {rsn_offset})");
            println!()
        }
    } else {
        let (starting_node, ending_node) = get_start_ending_position(
            recombination_col,
            &graph.nwp,
            &graph.pred_hash,
            &forward_matrix,
            forw_best_path,
            forw_ending_node,
            &rev_graph.nwp,
            &rev_graph.pred_hash,
            &reverse_matrix,
            rev_best_path,
            rev_starting_node,
        );
        let start_handle = graph.nodes_id_pos[starting_node];
        let mut counter = starting_node;
        let mut start_offset = 0;
        while graph.nodes_id_pos[counter] == start_handle {
            counter -= 1;
            start_offset += 1;
        }

        let end_handle = graph.nodes_id_pos[ending_node];
        let mut counter = ending_node;
        let mut end_offset = 0;
        while graph.nodes_id_pos[counter] == end_handle {
            counter -= 1;
            end_offset += 1;
        }

        let fen_handle = graph.nodes_id_pos[forw_ending_node];
        let mut counter = forw_ending_node;
        let mut fen_offset = 0;
        while graph.nodes_id_pos[counter] == fen_handle {
            counter -= 1;
            fen_offset += 1;
        }

        let rsn_handle = graph.nodes_id_pos[rev_starting_node];
        let mut counter = rev_starting_node;
        let mut rsn_offset = 0;
        while graph.nodes_id_pos[counter] == rsn_handle {
            counter -= 1;
            rsn_offset += 1;
        }
        if forw_best_path == rev_best_path {
            println!("No recombination, best path: {forw_best_path} (start: handle {start_handle} [{start_offset}]\tend: handle {end_handle} [{end_offset}])")
        } else {
            println!("Recombination between paths {forw_best_path} and {rev_best_path} (start: handle {start_handle} [{start_offset}]\tend: handle {end_handle} [{end_offset}])");
            println!("Recombination edge between node {fen_handle} (position: {fen_offset}) and node {rsn_handle} (position: {rsn_offset})");
        }
    }
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
    let mut results = vec![0; path_number];
    for (pred, paths) in pred_hash.get_preds_and_paths(lnz.len() - 1) {
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

fn absolute_scores(dpm: &mut Vec<Vec<Vec<i32>>>, alphas: &Vec<usize>, paths_nodes: &Vec<BitVec>) {
    for i in 0..dpm.len() {
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
) -> (usize, usize, usize, usize, usize) {
    let mut curr_best_score = 0;
    let mut forw_ending_node = 0;
    let mut rev_starting_node = 0;
    let mut forw_best_path = 0;
    let mut rev_best_path = 0;
    let mut recombination_col = 0;

    for j in 0..m[0].len() - 1 {
        for i in 0..m.len() - 1 {
            for rev_i in 0..m.len() {
                for forw_path in 0..m[0][0].len() {
                    for rev_path in 0..m[0][0].len() {
                        let penalty = if forw_path == rev_path {
                            0
                        } else {
                            brc + (mrc * dms[i][rev_i] as f32) as i32
                        };
                        if m[i][j][forw_path] + w[rev_i][j + 1][rev_path] - penalty
                            > curr_best_score
                        {
                            curr_best_score = m[i][j][forw_path] + w[rev_i][j + 1][rev_path];
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
    (
        forw_ending_node,
        rev_starting_node,
        forw_best_path,
        rev_best_path,
        recombination_col,
    )
}

fn get_start_ending_position(
    rec_col: usize,
    nwp: &BitVec,
    pred_hash: &PredHash,
    dpm: &Vec<Vec<Vec<i32>>>,
    best_path: usize,
    for_ending_node: usize,
    rev_nwp: &BitVec,
    rev_pred_hash: &PredHash,
    rev_dpm: &Vec<Vec<Vec<i32>>>,
    rev_path: usize,
    rev_ending_node: usize,
) -> (usize, usize) {
    // starting position
    let mut i = for_ending_node;
    let mut j = rec_col;
    while i > 0 && j > 0 {
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                dpm[i - 1][j - 1][best_path],
                dpm[i - 1][j][best_path],
                dpm[i][j - 1][best_path],
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    d = dpm[*pred][j - 1][best_path];
                    u = dpm[*pred][j][best_path];

                    l = dpm[i][j - 1][best_path];
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
        } else if max == u {
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
        } else {
            j -= 1;
        }
    }

    // ending position
    let mut rev_i = rev_ending_node;
    let mut j = rec_col + 1;
    while rev_i < rev_dpm.len() - 1 && j < rev_dpm[rev_i].len() {
        let mut predecessor = None;
        let (d, u, l) = if !rev_nwp[i] {
            (
                dpm[rev_i + 1][j + 1][rev_path],
                dpm[rev_i + 1][j][rev_path],
                dpm[rev_i][j + 1][rev_path],
            )
        } else {
            let preds = rev_pred_hash.get_preds_and_paths(rev_i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[rev_path] {
                    predecessor = Some(*pred);
                    d = dpm[*pred][j + 1][best_path];
                    u = dpm[*pred][j][best_path];

                    l = dpm[rev_i][j + 1][best_path];
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            rev_i = if predecessor.is_none() {
                rev_i + 1
            } else {
                predecessor.unwrap()
            };
            j += 1;
        } else if max == u {
            rev_i = if predecessor.is_none() {
                rev_i + 1
            } else {
                predecessor.unwrap()
            };
        } else {
            j += 1;
        }
    }
    (i, rev_i)
}
