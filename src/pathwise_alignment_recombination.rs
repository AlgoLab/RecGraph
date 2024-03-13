use bit_vec::BitVec;
use bstr::BString;

use crate::{
    dp_matrix::{DpDeltas, DpMatrix},
    gaf_output::GAFStruct,
    pathwise_graph::{PathGraph, PredHash},
    recombination_output,
    utils::{get_abs_val, idx},
};

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
    is_local: bool,
    sequence: &BString,
    graph: &PathGraph,
    rev_graph: &PathGraph,
    score_matrix: &Vec<i32>,
    base_rec_cost: i32,
    multi_rec_cost: f32,
    displacement_matrix: &Vec<Vec<i32>>,
    rec_number: i32,
) -> GAFStruct {
    let (forw_dpm, forw_deltas) = align(is_local, sequence, graph, score_matrix);

    let mut rev_dpm = DpMatrix::empty_new();
    let mut rev_deltas = DpDeltas::empty_new();
    let mut res = BestAlignStruct::new();

    if rec_number > 0 {
        let rev_sequence = get_rev_sequence(sequence);
        (rev_dpm, rev_deltas) = rev_align(is_local, &rev_sequence, rev_graph, score_matrix);
        res = best_alignment(
            &forw_dpm,
            &forw_deltas,
            &rev_dpm,
            &rev_deltas,
            &graph.alphas,
            &graph.paths_nodes,
            displacement_matrix,
            base_rec_cost,
            multi_rec_cost,
            is_local,
            &graph.pred_hash,
            &graph.nodes_id_pos,
        );
    }

    if rec_number > 0 && res.forw_best_path != res.rev_best_path {
        recombination_output::build_alignment_path_rec(
            &res,
            &forw_dpm,
            &forw_deltas,
            &rev_dpm,
            &rev_deltas,
            score_matrix,
            graph,
            rev_graph,
            sequence,
            is_local,
        )
    } else {
        let (max_score, best_path) = best_alignment_no_rec(
            &forw_dpm,
            &forw_deltas,
            &graph.alphas,
            &graph.paths_nodes,
            is_local,
            &graph.pred_hash,
        );
        recombination_output::build_alignment_path_no_rec(
            best_path,
            &forw_dpm,
            &forw_deltas,
            score_matrix,
            graph,
            sequence,
            (max_score as f32, 0),
            is_local,
        )
    }
}

fn rev_align(
    is_local: bool,
    sequence: &BString,
    graph: &PathGraph,
    score_matrix: &Vec<i32>,
) -> (DpMatrix, DpDeltas) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let path_number = graph.paths_number;
    let path_node = &graph.paths_nodes;

    let mut dpm = DpMatrix::new(lnz.len(), sequence.len());
    let mut deltas = DpDeltas::new(lnz.len(), sequence.len());
    let alphas = &graph.alphas;

    let last_node_pos = lnz.len() - 1;
    let last_char_pos = sequence.len() - 1;

    for i in (1..last_node_pos + 1).rev() {
        for j in (0..=last_char_pos).rev() {
            if i == last_node_pos && j == last_char_pos {
                dpm.set(i, j, 0);
                deltas.set_vals(i, j, vec![0; path_number]);
            } else if i == last_node_pos {
                let value = dpm.get(i, j + 1) + score_matrix[idx(sequence[j], b'-')];
                dpm.set(i, j, value);
                deltas.set_coord(i, j, (last_node_pos, last_char_pos));
            } else if j == last_char_pos {
                if is_local {
                    dpm.set(i, j, 0);
                    deltas.set_coord(i, j, (last_node_pos, last_char_pos));
                } else if !nodes_with_pred[i] {
                    dpm.set(i, j, dpm.get(i + 1, j) + score_matrix[idx(lnz[i], b'-')]);
                    deltas.set_coord(i, j, deltas.get_pred_coord(i + 1, j));
                } else {
                    let mut absolute_scores = vec![0; path_number];

                    for (p, p_paths) in pred_hash.get_preds_and_paths(i) {
                        let mut common_paths = path_node[i].clone();
                        common_paths.and(&p_paths);

                        if common_paths[alphas[p]] {
                            absolute_scores[alphas[p]] =
                                dpm.get(p, j) + score_matrix[idx(b'-', lnz[i])];

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in && path != alphas[p] {
                                    absolute_scores[path] =
                                        get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', lnz[i])];
                                }
                            }
                        } else {
                            //set new alpha

                            let temp_alpha = if common_paths[alphas[i]] {
                                alphas[i]
                            } else {
                                common_paths.iter().position(|is_in| is_in).unwrap()
                            };

                            absolute_scores[temp_alpha] =
                                get_abs_val(p, j, alphas, temp_alpha, &dpm, &deltas)
                                    + score_matrix[idx(b'-', lnz[i])];

                            for (path, is_in) in common_paths.iter().enumerate() {
                                if is_in && path != temp_alpha {
                                    absolute_scores[path] =
                                        get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', lnz[i])];
                                }
                            }
                        }
                    }
                    // restore scores
                    dpm.set(i, j, absolute_scores[alphas[i]]);
                    for (path, is_in) in path_node[i].iter().enumerate() {
                        if is_in && path != alphas[i] {
                            absolute_scores[path] -= absolute_scores[alphas[i]];
                        }
                    }
                    absolute_scores[alphas[i]] = 0;

                    deltas.set_vals(i, j, absolute_scores);
                }
            } else if !nodes_with_pred[i] {
                let u = dpm.get(i + 1, j) + score_matrix[idx(b'-', lnz[i])];
                let d = dpm.get(i + 1, j + 1) + score_matrix[idx(sequence[j], lnz[i])];
                let l = dpm.get(i, j + 1) + score_matrix[idx(b'-', sequence[j])];

                let max = *[d, u, l].iter().max().unwrap();
                dpm.set(i, j, max);

                let coord = if max == d {
                    deltas.get_pred_coord(i + 1, j + 1)
                } else if max == u {
                    deltas.get_pred_coord(i + 1, j)
                } else {
                    deltas.get_pred_coord(i, j + 1)
                };
                deltas.set_coord(i, j, coord);
            } else {
                // multiple alphas possible
                let mut absolute_scores = vec![0; path_number];

                for (p, p_paths) in pred_hash.get_preds_and_paths(i) {
                    let mut common_paths = path_node[i].clone();
                    common_paths.and(&p_paths);

                    if common_paths[alphas[p]] {
                        let u = dpm.get(p, j) + score_matrix[idx(b'-', lnz[i])];

                        let d = dpm.get(p, j + 1) + score_matrix[idx(sequence[j], lnz[i])];

                        let l = get_abs_val(i, j + 1, alphas, alphas[p], &dpm, &deltas)
                            + score_matrix[idx(b'-', sequence[j])];

                        let max = *[d, u, l].iter().max().unwrap();

                        absolute_scores[alphas[p]] = max;

                        for (path, is_in) in common_paths.iter().enumerate() {
                            if is_in && path != alphas[p] {
                                if max == d {
                                    absolute_scores[path] =
                                        get_abs_val(p, j + 1, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(sequence[j], lnz[i])];
                                } else if max == u {
                                    absolute_scores[path] =
                                        get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', lnz[i])];
                                } else {
                                    absolute_scores[path] =
                                        get_abs_val(i, j + 1, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', sequence[j])];
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

                        let u = get_abs_val(p, j, alphas, temp_alpha, &dpm, &deltas)
                            + score_matrix[idx(b'-', lnz[i])];

                        let d = get_abs_val(p, j + 1, alphas, temp_alpha, &dpm, &deltas)
                            + score_matrix[idx(sequence[j], lnz[i])];

                        let l = get_abs_val(i, j + 1, alphas, temp_alpha, &dpm, &deltas)
                            + score_matrix[idx(b'-', sequence[j])];

                        let max = *[d, u, l].iter().max().unwrap();
                        absolute_scores[temp_alpha] = max;

                        for (path, is_in) in common_paths.iter().enumerate() {
                            if path != temp_alpha && is_in {
                                if max == d {
                                    absolute_scores[path] =
                                        get_abs_val(p, j + 1, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(sequence[j], lnz[i])];
                                } else if max == u {
                                    absolute_scores[path] =
                                        get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', lnz[i])];
                                } else {
                                    absolute_scores[path] =
                                        get_abs_val(i, j + 1, alphas, path, &dpm, &deltas)
                                            + score_matrix[idx(b'-', sequence[j])];
                                }
                            }
                        }
                    }
                }
                dpm.set(i, j, absolute_scores[alphas[i]]);
                for (path, is_in) in path_node[i].iter().enumerate() {
                    if is_in && path != alphas[i] {
                        absolute_scores[path] -= absolute_scores[alphas[i]];
                    }
                }
                absolute_scores[alphas[i]] = 0;
                deltas.set_vals(i, j, absolute_scores);
            }
        }
    }
    (dpm, deltas)
}

fn align(
    is_local: bool,
    sequence: &BString,
    graph: &PathGraph,
    score_matrix: &Vec<i32>,
) -> (DpMatrix, DpDeltas) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let path_number = graph.paths_number;
    let path_node = &graph.paths_nodes;

    //let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let mut dpm = DpMatrix::new(lnz.len(), sequence.len());
    let mut deltas = DpDeltas::new(lnz.len(), sequence.len());

    let alphas = &graph.alphas;
    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {
                    dpm.set(i, j, 0);
                    deltas.set_vals(i, j, vec![0; path_number]);
                }
                (_, 0) => {
                    if is_local {
                        dpm.set(i, j, 0);
                        deltas.set_coord(i, j, (0, 0));
                    } else if !nodes_with_pred[i] {
                        dpm.set(i, j, dpm.get(i - 1, j) + score_matrix[idx(lnz[i], b'-')]);
                        deltas.set_coord(i, j, deltas.get_pred_coord(i - 1, j));
                    } else {
                        let mut absolute_scores = vec![0; path_number];

                        for (p, p_paths) in pred_hash.get_preds_and_paths(i) {
                            let mut common_paths = path_node[i].clone();
                            common_paths.and(&p_paths);

                            if common_paths[alphas[p]] {
                                absolute_scores[alphas[p]] =
                                    dpm.get(p, j) + score_matrix[idx(b'-', lnz[i])];

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[p] {
                                        absolute_scores[path] =
                                            get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                                + score_matrix[idx(b'-', lnz[i])];
                                    }
                                }
                            } else {
                                //set new alpha

                                let temp_alpha = if common_paths[alphas[i]] {
                                    alphas[i]
                                } else {
                                    common_paths.iter().position(|is_in| is_in).unwrap()
                                };

                                absolute_scores[temp_alpha] =
                                    get_abs_val(p, j, alphas, temp_alpha, &dpm, &deltas)
                                        + score_matrix[idx(b'-', lnz[i])];

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != temp_alpha {
                                        absolute_scores[path] =
                                            get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                                + score_matrix[idx(b'-', lnz[i])];
                                    }
                                }
                            }
                        }
                        // restore scores
                        dpm.set(i, j, absolute_scores[alphas[i]]);
                        for (path, is_in) in path_node[i].iter().enumerate() {
                            if is_in && path != alphas[i] {
                                absolute_scores[path] -= absolute_scores[alphas[i]];
                            }
                        }
                        absolute_scores[alphas[i]] = 0;

                        deltas.set_vals(i, j, absolute_scores);
                    }
                }
                (0, _) => {
                    let value = dpm.get(i, j - 1) + score_matrix[idx(sequence[j], b'-')];
                    dpm.set(i, j, value);
                    deltas.set_coord(i, j, (0, 0));
                }
                _ => {
                    if !nodes_with_pred[i] {
                        let u = dpm.get(i - 1, j) + score_matrix[idx(b'-', lnz[i])];

                        let d = dpm.get(i - 1, j - 1) + score_matrix[idx(sequence[j], lnz[i])];

                        let l = dpm.get(i, j - 1) + score_matrix[idx(b'-', sequence[j])];

                        let max = *[d, u, l].iter().max().unwrap();
                        dpm.set(i, j, max);

                        let coord = if max == d {
                            deltas.get_pred_coord(i - 1, j - 1)
                        } else if max == u {
                            deltas.get_pred_coord(i - 1, j)
                        } else {
                            deltas.get_pred_coord(i, j - 1)
                        };
                        deltas.set_coord(i, j, coord);
                    } else {
                        // multiple alphas possible
                        let mut absolute_scores = vec![0; path_number];

                        for (p, p_paths) in pred_hash.get_preds_and_paths(i) {
                            let mut common_paths = path_node[i].clone();
                            common_paths.and(&p_paths);

                            if common_paths[alphas[p]] {
                                let u = dpm.get(p, j) + score_matrix[idx(b'-', lnz[i])];

                                let d = dpm.get(p, j - 1) + score_matrix[idx(sequence[j], lnz[i])];

                                let l = get_abs_val(i, j - 1, alphas, alphas[p], &dpm, &deltas)
                                    + score_matrix[idx(b'-', sequence[j])];

                                let max = *[d, u, l].iter().max().unwrap();

                                absolute_scores[alphas[p]] = max;
                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if is_in && path != alphas[p] {
                                        if max == d {
                                            absolute_scores[path] =
                                                get_abs_val(p, j - 1, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(sequence[j], lnz[i])];
                                        } else if max == u {
                                            absolute_scores[path] =
                                                get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(b'-', lnz[i])];
                                        } else {
                                            absolute_scores[path] =
                                                get_abs_val(i, j - 1, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(b'-', sequence[j])];
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

                                let u = get_abs_val(p, j, alphas, temp_alpha, &dpm, &deltas)
                                    + score_matrix[idx(b'-', lnz[i])];

                                let d = get_abs_val(p, j - 1, alphas, temp_alpha, &dpm, &deltas)
                                    + score_matrix[idx(sequence[j], lnz[i])];

                                let l = get_abs_val(i, j - 1, alphas, temp_alpha, &dpm, &deltas)
                                    + score_matrix[idx(b'-', sequence[j])];

                                let max = *[d, u, l].iter().max().unwrap();
                                absolute_scores[temp_alpha] = max;

                                for (path, is_in) in common_paths.iter().enumerate() {
                                    if path != temp_alpha && is_in {
                                        if max == d {
                                            absolute_scores[path] =
                                                get_abs_val(p, j - 1, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(sequence[j], lnz[i])];
                                        } else if max == u {
                                            absolute_scores[path] =
                                                get_abs_val(p, j, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(b'-', lnz[i])];
                                        } else {
                                            absolute_scores[path] =
                                                get_abs_val(i, j - 1, alphas, path, &dpm, &deltas)
                                                    + score_matrix[idx(b'-', sequence[j])];
                                        }
                                    }
                                }
                            }
                        }
                        // restore scores
                        dpm.set(i, j, absolute_scores[alphas[i]]);
                        for (path, is_in) in path_node[i].iter().enumerate() {
                            if is_in && path != alphas[i] {
                                absolute_scores[path] -= absolute_scores[alphas[i]];
                            }
                        }
                        absolute_scores[alphas[i]] = 0;
                        deltas.set_vals(i, j, absolute_scores);
                    }
                }
            }
        }
    }
    (dpm, deltas)
}

fn best_scores(
    dpm: &DpMatrix,
    deltas: &DpDeltas,
    alphas: &Vec<usize>,
    paths_nodes: &Vec<BitVec>,
) -> Vec<(usize, i32)> {
    let mut result: Vec<(usize, i32)> = vec![(0, 0); dpm.rows * dpm.cols];
    for i in 1..dpm.rows - 1 {
        for j in 1..dpm.cols {
            let mut max = dpm.get(i, j);
            let mut best_path = alphas[i];

            for (path, is_in) in paths_nodes[i].iter().enumerate() {
                if is_in && deltas.get_val(i, j, path) > 0 {
                    let abs_val = get_abs_val(i, j, alphas, path, dpm, deltas);
                    if abs_val > max {
                        max = abs_val;
                        best_path = path;
                    }
                }
            }
            result[j * dpm.rows + i] = (best_path, max);
        }
    }
    result
}

fn best_alignment(
    dpm: &DpMatrix,
    deltas: &DpDeltas,
    rev_dpm: &DpMatrix,
    rev_deltas: &DpDeltas,
    alphas: &Vec<usize>,
    paths_nodes: &Vec<BitVec>,
    dms: &Vec<Vec<i32>>,
    brc: i32,
    mrc: f32,
    is_local: bool,
    pred_hash: &PredHash,
    nodes_id_pos: &Vec<u64>,
) -> BestAlignStruct {
    let seq_len = dpm.cols;
    let lnz_len = dpm.rows;

    let m = best_scores(dpm, deltas, alphas, paths_nodes);
    let w = best_scores(rev_dpm, rev_deltas, alphas, paths_nodes);

    let mut max = None;
    let mut best_path = None;
    if !is_local {
        let preds = pred_hash.get_preds_and_paths(lnz_len - 1);
        for (pred, _) in preds.iter() {
            if max.is_none() || max.unwrap() < m[(seq_len - 1) * lnz_len + pred].1 {
                max = Some(m[(seq_len - 1) * lnz_len + pred].1);
                best_path = Some(m[(seq_len - 1) * lnz_len + pred].0);
            }
        }
    } else {
        for i in 0..lnz_len - 1 {
            if max.is_none() || max.unwrap() < m[(seq_len - 1) * lnz_len + i].1 {
                max = Some(m[(seq_len - 1) * lnz_len + i].1);
                best_path = Some(m[(seq_len - 1) * lnz_len + i].0);
            }
        }
    }

    let mut result = BestAlignStruct::new();
        result.curr_best_score = max.unwrap() as f32;
        result.forw_best_path = best_path.unwrap();
        result.rev_best_path = best_path.unwrap();

    for j in 0..seq_len {
        // check recomb only if score increment is possible
        let forw_is = &m[j * lnz_len + 1..(j + 1) * lnz_len - 1];
        let rev_is = &w[j * lnz_len + 1..(j + 1) * lnz_len - 1];
        let forw_i_max = forw_is.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap().1;
        let rev_i_max = rev_is.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap().1;

        if (forw_i_max  + rev_i_max  - brc) as f32> result.curr_best_score {
            for (forw_idx, (forw_path, forw_score)) in forw_is.iter().enumerate() {
                let i = forw_idx + 1;
                for (rev_idx, (rev_path, rev_score)) in rev_is.iter().enumerate() {
                    let rev_i = rev_idx + 1;
                    if nodes_id_pos[i] != nodes_id_pos[rev_i] && forw_path != rev_path {
                        let penalty = brc as f32 + (mrc * dms[i][rev_i] as f32);
                        let new_score = (forw_score + rev_score) as f32 - penalty;
                        {
                        if new_score > result.curr_best_score
                            || (new_score == result.curr_best_score
                                && !result.onedge
                                && (i + 1 == m.len() || nodes_id_pos[i] != nodes_id_pos[i + 1])
                                && nodes_id_pos[rev_i] != nodes_id_pos[rev_i - 1])
                        {
                            let onedge = (i + 1 == m.len() || nodes_id_pos[i] != nodes_id_pos[i + 1])
                                && nodes_id_pos[rev_i] != nodes_id_pos[rev_i - 1];
                            result.update(new_score, i, rev_i, forw_path, rev_path, j, penalty, onedge);
                        }
                        }
                        
                    }
                }
            }
        }
    }
    result
}

pub fn get_rev_sequence(seq: &BString) -> BString {
    let mut rev_seq = seq.clone();
    rev_seq.remove(0);
    rev_seq.push(b'F');

    rev_seq
}

fn best_alignment_no_rec(
    dpm: &DpMatrix,
    deltas: &DpDeltas,
    alphas: &Vec<usize>,
    paths_nodes: &Vec<BitVec>,
    is_local: bool,
    pred_hash: &PredHash,
) -> (i32, usize) {
    let seq_len = dpm.cols;
    let lnz_len = dpm.rows;

    let m = best_scores(dpm, deltas, alphas, paths_nodes);

    let mut max = None;
    let mut best_path = None;
    if !is_local {
        let preds = pred_hash.get_preds_and_paths(lnz_len - 1);
        for (pred, _) in preds.iter() {
            if max.is_none() || max.unwrap() < m[(seq_len - 1) * lnz_len + pred].1 {
                max = Some(m[(seq_len - 1) * lnz_len + pred].1);
                best_path = Some(m[(seq_len - 1) * lnz_len + pred].0);
            }
        }
    } else {
        for i in 0..lnz_len - 1 {
            if max.is_none() || max.unwrap() < m[(seq_len - 1) * lnz_len + i].1 {
                max = Some(m[(seq_len - 1) * lnz_len + i].1);
                best_path = Some(m[(seq_len - 1) * lnz_len + i].0);
            }
        }
    }

    (max.unwrap(), best_path.unwrap())
}

#[derive(Debug)]
pub struct BestAlignStruct {
    pub forw_ending_node: usize,
    pub rev_starting_node: usize,
    pub forw_best_path: usize,
    pub rev_best_path: usize,
    pub recombination_col: usize,
    pub curr_best_score: f32,
    pub rec_penalty: f32,
    pub onedge: bool,
}

impl BestAlignStruct {
    pub fn new() -> BestAlignStruct {
        BestAlignStruct {
            forw_ending_node: 0,
            rev_starting_node: 0,
            forw_best_path: 0,
            rev_best_path: 0,
            recombination_col: 0,
            curr_best_score: 0f32,
            rec_penalty: 0f32,
            onedge: false,
        }
    }

    pub fn update(
        &mut self,
        new_score: f32,
        i: usize,
        rev_i: usize,
        forw_path: &usize,
        rev_path: &usize,
        j: usize,
        rec_penalty: f32,
        onedge: bool,
    ) {
        self.curr_best_score = new_score;
        self.forw_ending_node = i;
        self.rev_starting_node = rev_i;
        self.forw_best_path = *forw_path;
        self.rev_best_path = *rev_path;
        self.recombination_col = j;
        self.rec_penalty = rec_penalty;
        self.onedge = onedge;
    }
}
