use bit_vec::BitVec;
use bstr::BString;

use crate::{
    build_cigar::build_cigar,
    dp_matrix::{DpDeltas, DpMatrix},
    gaf_output::GAFStruct,
    pathwise_alignment_recombination::{get_node_offset, get_rev_sequence, BestAlignStruct},
    pathwise_graph::PathGraph,
    utils::{self, get_abs_val, idx},
};

pub fn build_alignment_path_rec(
    res: &BestAlignStruct,
    forw_dpm: &DpMatrix,
    forw_deltas: &DpDeltas,
    rev_dpm: &DpMatrix,
    rev_deltas: &DpDeltas,
    scores: &Vec<i32>,
    graph: &PathGraph,
    rev_graph: &PathGraph,
    seq: &BString,
    is_local: bool,
) -> GAFStruct {
    let seq_len = seq.len();
    let lnz_len = graph.lnz.len();

    // REC
    // reverse alignment
    let r_seq = &get_rev_sequence(seq);

    let mut cigar = Vec::new();
    let mut rev_path_length: usize = 0;
    let mut i = res.rev_starting_node;
    let mut j = res.recombination_col;
    let mut handle_id_alignment = Vec::new();
    let mut rev_ending_node = i;
    let mut path_sequence = Vec::new();

    while i > 0 && i < rev_graph.lnz.len() - 1 && j < seq.len() - 1 {
        let curr_score = get_abs_val(
            i,
            j,
            &rev_graph.alphas,
            res.rev_best_path,
            rev_dpm,
            rev_deltas,
        );

        let (d, u, l) = if rev_graph.nwp[i] {
            let l_val_pos = (
                get_abs_val(
                    i,
                    j + 1,
                    &rev_graph.alphas,
                    res.rev_best_path,
                    rev_dpm,
                    rev_deltas,
                ) + scores[idx(b'-', r_seq[j])],
                (i, j + 1),
            );
            let mut u_val_pos = (0, (0, 0));
            let mut d_val_pos = (0, (0, 0));
            for (pred, paths) in rev_graph.pred_hash.get_preds_and_paths(i) {
                if paths[res.rev_best_path] {
                    u_val_pos = (
                        get_abs_val(
                            pred,
                            j,
                            &rev_graph.alphas,
                            res.rev_best_path,
                            rev_dpm,
                            rev_deltas,
                        ) + scores[idx(b'-', rev_graph.lnz[i])],
                        (pred, j),
                    );
                    d_val_pos = (
                        get_abs_val(
                            pred,
                            j + 1,
                            &rev_graph.alphas,
                            res.rev_best_path,
                            rev_dpm,
                            rev_deltas,
                        ) + scores[idx(r_seq[j], rev_graph.lnz[i])],
                        (pred, j + 1),
                    );
                }
            }
            (d_val_pos, u_val_pos, l_val_pos)
        } else {
            let u_val_pos = (
                get_abs_val(
                    i + 1,
                    j,
                    &rev_graph.alphas,
                    res.rev_best_path,
                    rev_dpm,
                    rev_deltas,
                ) + scores[idx(b'-', rev_graph.lnz[i])],
                (i + 1, j),
            );
            let d_val_pos = (
                get_abs_val(
                    i + 1,
                    j + 1,
                    &rev_graph.alphas,
                    res.rev_best_path,
                    rev_dpm,
                    rev_deltas,
                ) + scores[idx(r_seq[j], rev_graph.lnz[i])],
                (i + 1, j + 1),
            );
            let l_val_pos = (
                get_abs_val(
                    i,
                    j + 1,
                    &rev_graph.alphas,
                    res.rev_best_path,
                    rev_dpm,
                    rev_deltas,
                ) + scores[idx(b'-', r_seq[j])],
                (i, j + 1),
            );
            (d_val_pos, u_val_pos, l_val_pos)
        };

        rev_ending_node = i;

        (i, j, rev_path_length) = update_alignment(
            i,
            j,
            curr_score,
            d,
            u,
            l,
            rev_graph,
            r_seq,
            &mut cigar,
            rev_path_length,
            &mut path_sequence,
            &mut handle_id_alignment,
        );
    }

    while j < seq_len - 1 {
        cigar.push('L');
        j += 1;
    }

    // only if global
    if !is_local {
        while i < lnz_len - 1 {
            cigar.push('U');
            handle_id_alignment.push(rev_graph.nodes_id_pos[i]);
            path_sequence.push(rev_graph.lnz[i]);
            let mut predecessor = None;
            if rev_graph.nwp[i] {
                let preds = rev_graph.pred_hash.get_preds_and_paths(i);
                for (pred, paths) in preds.iter() {
                    if paths[res.rev_best_path] {
                        predecessor = Some(*pred);
                    }
                }
            }
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            rev_path_length += 1;
        }
    }

    // forward alignment

    let mut path_length: usize = 0;
    let mut temp_cigar = Vec::new();
    let mut temp_handle_id_alignment = Vec::new();
    let mut temp_path_sequence = Vec::new();

    i = res.forw_ending_node;
    j = res.recombination_col;

    while i > 0 && j > 0 {
        let curr_score = get_abs_val(
            i,
            j,
            &graph.alphas,
            res.forw_best_path,
            forw_dpm,
            forw_deltas,
        );

        let (d, u, l) = if graph.nwp[i] {
            let l_val_pos = (
                get_abs_val(
                    i,
                    j - 1,
                    &graph.alphas,
                    res.forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', seq[j])],
                (i, j - 1),
            );
            let mut u_val_pos = (0, (0, 0));
            let mut d_val_pos = (0, (0, 0));
            for (pred, paths) in graph.pred_hash.get_preds_and_paths(i) {
                if paths[res.forw_best_path] {
                    u_val_pos = (
                        get_abs_val(
                            pred,
                            j,
                            &graph.alphas,
                            res.forw_best_path,
                            forw_dpm,
                            forw_deltas,
                        ) + scores[idx(b'-', graph.lnz[i])],
                        (pred, j),
                    );
                    d_val_pos = (
                        get_abs_val(
                            pred,
                            j - 1,
                            &graph.alphas,
                            res.forw_best_path,
                            forw_dpm,
                            forw_deltas,
                        ) + scores[idx(seq[j], graph.lnz[i])],
                        (pred, j - 1),
                    );
                }
            }
            (d_val_pos, u_val_pos, l_val_pos)
        } else {
            let u_val_pos = (
                get_abs_val(
                    i - 1,
                    j,
                    &graph.alphas,
                    res.forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', graph.lnz[i])],
                (i - 1, j),
            );
            let d_val_pos = (
                get_abs_val(
                    i - 1,
                    j - 1,
                    &graph.alphas,
                    res.forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(seq[j], graph.lnz[i])],
                (i - 1, j - 1),
            );
            let l_val_pos = (
                get_abs_val(
                    i,
                    j - 1,
                    &graph.alphas,
                    res.forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', seq[j])],
                (i, j - 1),
            );
            (d_val_pos, u_val_pos, l_val_pos)
        };

        (i, j, path_length) = update_alignment(
            i,
            j,
            curr_score,
            d,
            u,
            l,
            graph,
            seq,
            &mut temp_cigar,
            path_length,
            &mut temp_path_sequence,
            &mut temp_handle_id_alignment,
        );
    }
    while j > 0 {
        temp_cigar.push('L');
        j -= 1;
    }

    // only if global
    if !is_local {
        while i > 0 {
            temp_cigar.push('U');
            temp_handle_id_alignment.push(graph.nodes_id_pos[i]);
            temp_path_sequence.push(graph.lnz[i]);
            let mut predecessor = None;
            if graph.nwp[i] {
                let preds = graph.pred_hash.get_preds_and_paths(i);
                for (pred, paths) in preds.iter() {
                    if paths[res.forw_best_path] {
                        predecessor = Some(*pred);
                    }
                }
            }
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        }
    }

    let rec_edge = temp_path_sequence.len() - 1;

    temp_cigar.reverse();
    temp_cigar.append(&mut cigar);

    temp_handle_id_alignment.reverse();
    temp_handle_id_alignment.append(&mut handle_id_alignment);
    temp_handle_id_alignment.dedup();

    temp_path_sequence.reverse();
    temp_path_sequence.append(&mut path_sequence);

    let path_sequence_string: String = String::from_utf8(temp_path_sequence).unwrap();

    let query_name = String::from("Temp");
    let seq_length = seq.len() - 1;
    let query_start = 0;
    let query_end = seq.len() - 2;
    let strand = '+';
    let path: Vec<usize> = temp_handle_id_alignment
        .iter()
        .map(|id| *id as usize)
        .collect();
    let start = if i == 0 { i } else { i + 1 };
    let end = rev_ending_node;
    let (path_len, path_start, path_end) = utils::get_rec_path_len_start_end(
        &graph.nodes_id_pos,
        res.forw_ending_node,
        res.rev_starting_node,
        start,
        end,
        path_length,
        rev_path_length,
    );
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let recombination = if res.forw_best_path == res.rev_best_path {
        format!("No recombination, best path: {}", res.forw_best_path)
    } else {
        let fen_offset = get_node_offset(&graph.nodes_id_pos, res.forw_ending_node);
        let rsn_offset = get_node_offset(&graph.nodes_id_pos, res.rev_starting_node);
        format!(
            "recombination path {} {}, nodes {}[{}] {}[{}], score: {}, displacement: {}\t{}\t{}",
            res.forw_best_path,
            res.rev_best_path,
            graph.nodes_id_pos[res.forw_ending_node],
            fen_offset,
            graph.nodes_id_pos[res.rev_starting_node],
            rsn_offset,
            res.curr_best_score,
            res.rec_penalty,
            path_sequence_string,
            rec_edge
        )
    };
    let comments = format!("{}, {}", build_cigar(&temp_cigar), recombination);

    GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_len,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    )
}

pub fn build_alignment_path_no_rec(
    forw_best_path: usize,
    forw_dpm: &DpMatrix,
    forw_deltas: &DpDeltas,
    scores: &Vec<i32>,
    graph: &PathGraph,
    seq: &BString,
    best_score: (f32, i32),
    is_local: bool,
) -> GAFStruct {
    let seq_len = seq.len();

    let mut path_length: usize = 0;
    let mut temp_cigar = Vec::new();
    let mut temp_handle_id_alignment = Vec::new();
    let mut temp_path_sequence = Vec::new();

    let mut i = 0;
    let ending_nodes = graph.pred_hash.get_preds_and_paths(graph.lnz.len() - 1);
    for (node, paths) in ending_nodes.iter() {
        if paths[forw_best_path] {
            i = *node;
        }
    }
    if is_local {
        i = ending_node(
            forw_dpm,
            forw_deltas,
            forw_best_path,
            &graph.paths_nodes,
            seq_len,
            graph,
        );
    }
    let ending_node = i;
    let mut j = seq.len() - 1;

    while i > 0 && j > 0 {
        let curr_score = get_abs_val(i, j, &graph.alphas, forw_best_path, forw_dpm, forw_deltas);

        let (d, u, l) = if graph.nwp[i] {
            let l_val_pos = (
                get_abs_val(
                    i,
                    j - 1,
                    &graph.alphas,
                    forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', seq[j])],
                (i, j - 1),
            );
            let mut u_val_pos = (0, (0, 0));
            let mut d_val_pos = (0, (0, 0));
            for (pred, paths) in graph.pred_hash.get_preds_and_paths(i) {
                if paths[forw_best_path] {
                    u_val_pos = (
                        get_abs_val(
                            pred,
                            j,
                            &graph.alphas,
                            forw_best_path,
                            forw_dpm,
                            forw_deltas,
                        ) + scores[idx(b'-', graph.lnz[i])],
                        (pred, j),
                    );
                    d_val_pos = (
                        get_abs_val(
                            pred,
                            j - 1,
                            &graph.alphas,
                            forw_best_path,
                            forw_dpm,
                            forw_deltas,
                        ) + scores[idx(seq[j], graph.lnz[i])],
                        (pred, j - 1),
                    );
                }
            }
            (d_val_pos, u_val_pos, l_val_pos)
        } else {
            let u_val_pos = (
                get_abs_val(
                    i - 1,
                    j,
                    &graph.alphas,
                    forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', graph.lnz[i])],
                (i - 1, j),
            );
            let d_val_pos = (
                get_abs_val(
                    i - 1,
                    j - 1,
                    &graph.alphas,
                    forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(seq[j], graph.lnz[i])],
                (i - 1, j - 1),
            );
            let l_val_pos = (
                get_abs_val(
                    i,
                    j - 1,
                    &graph.alphas,
                    forw_best_path,
                    forw_dpm,
                    forw_deltas,
                ) + scores[idx(b'-', seq[j])],
                (i, j - 1),
            );
            (d_val_pos, u_val_pos, l_val_pos)
        };

        (i, j, path_length) = update_alignment(
            i,
            j,
            curr_score,
            d,
            u,
            l,
            graph,
            seq,
            &mut temp_cigar,
            path_length,
            &mut temp_path_sequence,
            &mut temp_handle_id_alignment,
        );
        //ln!("i: {}, j: {}\tPATH:{:?}", i, j, temp_handle_id_alignment);
    }
    while j > 0 {
        temp_cigar.push('L');
        j -= 1;
    }

    // only if global
    if !is_local {
        while i > 0 {
            temp_cigar.push('U');
            temp_handle_id_alignment.push(graph.nodes_id_pos[i]);
            temp_path_sequence.push(graph.lnz[i]);
            let mut predecessor = None;
            if graph.nwp[i] {
                let preds = graph.pred_hash.get_preds_and_paths(i);
                for (pred, paths) in preds.iter() {
                    if paths[forw_best_path] {
                        predecessor = Some(*pred);
                    }
                }
            }
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        }
    }
    temp_cigar.reverse();
    temp_path_sequence.reverse();

    let path_sequence_string: String = String::from_utf8(temp_path_sequence).unwrap();

    let query_name = String::from("Temp");
    let seq_length = seq.len() - 1;
    let query_start = 0;
    let query_end = seq.len() - 2;
    let strand = '+';

    temp_handle_id_alignment.dedup();
    temp_handle_id_alignment.reverse();
    let path: Vec<usize> = temp_handle_id_alignment
        .iter()
        .map(|id| *id as usize)
        .collect();
    let (path_len, path_start, path_end) = utils::get_path_len_start_end(
        &graph.nodes_id_pos,
        if i == 0 { i } else { i + 1 },
        ending_node,
        path_length,
    );

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = format!(
        "{}, best path: {}, score: {}\t{}",
        build_cigar(&temp_cigar),
        forw_best_path,
        best_score.0,
        path_sequence_string
    );

    GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_len,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    )
}

fn update_alignment(
    i: usize,
    j: usize,
    curr_score: i32,
    d: (i32, (usize, usize)),
    u: (i32, (usize, usize)),
    l: (i32, (usize, usize)),
    graph: &PathGraph,
    seq: &BString,
    cigar: &mut Vec<char>,
    path_length: usize,
    path_sequence: &mut Vec<u8>,
    handle_id_alignment: &mut Vec<u64>,
) -> (usize, usize, usize) {
    let mut new_path_length = path_length;
    let (new_i, new_j) = if curr_score == d.0 {
        if graph.lnz[i] != seq[j] {
            cigar.push('d');
        } else {
            cigar.push('D');
        }
        new_path_length += 1;
        path_sequence.push(graph.lnz[i]);
        handle_id_alignment.push(graph.nodes_id_pos[i]);
        (d.1 .0, d.1 .1)
    } else if curr_score == u.0 {
        cigar.push('U');
        new_path_length += 1;
        path_sequence.push(graph.lnz[i]);
        handle_id_alignment.push(graph.nodes_id_pos[i]);
        (u.1 .0, u.1 .1)
    } else {
        cigar.push('L');
        new_path_length += 1;
        path_sequence.push(b'-');
        (l.1 .0, l.1 .1)
    };
    (new_i, new_j, new_path_length)
}

fn ending_node(
    dpm: &DpMatrix,
    deltas: &DpDeltas,
    best_path: usize,
    paths_nodes: &Vec<BitVec>,
    seq_len: usize,
    graph: &PathGraph,
) -> usize {
    let mut best_score = None;
    let mut best_node = 0;
    for i in 1..graph.lnz.len() - 1 {
        if paths_nodes[i][best_path] {
            let score = get_abs_val(i, seq_len - 1, &graph.alphas, best_path, dpm, deltas);
            if best_score.is_none() || score > best_score.unwrap() {
                best_score = Some(score);
                best_node = i;
            }
        }
    }
    best_node
}
