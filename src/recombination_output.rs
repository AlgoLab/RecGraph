use bit_vec::BitVec;

use crate::{
    gaf_output::GAFStruct,
    pathwise_alignment_output::build_cigar,
    pathwise_alignment_recombination::{get_node_offset, get_rev_sequence},
    pathwise_graph::{PathGraph, PredHash},
    utils,
};
use std::collections::HashMap;

pub fn gaf_output_global_rec(
    dpm: &Vec<Vec<Vec<i32>>>,
    rev_dpm: &Vec<Vec<Vec<i32>>>,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    best_path: usize,
    rev_best_path: usize,
    forward_ending_node: usize,
    reverse_starting_node: usize,
    rec_col: usize,
    lnz: &Vec<char>,
    pred_hash: &PredHash,
    rev_pred_hash: &PredHash,
    nwp: &BitVec,
    rev_nwp: &BitVec,
    handles_nodes_id: &Vec<u64>,
) -> GAFStruct {
    let mut cigar = Vec::new();
    let mut rev_path_length: usize = 0;
    let mut i = reverse_starting_node;
    let mut j = rec_col;
    let mut handle_id_alignment = Vec::new();
    let mut rev_ending_node = 0;
    let ending_nodes = pred_hash.get_preds_and_paths(dpm.len() - 1);
    for (node, paths) in ending_nodes.iter() {
        if paths[rev_best_path] {
            rev_ending_node = *node
        }
    }

    // reverse alignment
    let r_seq = &get_rev_sequence(seq);
    while i > 0 && i < dpm.len() - 1 && j < dpm[0].len() - 1 {
        let mut predecessor = None;
        let (d, u, l) = if !rev_nwp[i] {
            (
                rev_dpm[i + 1][j + 1][rev_best_path] + scores.get(&(lnz[i], r_seq[j])).unwrap(),
                rev_dpm[i + 1][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap(),
                rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', r_seq[j])).unwrap(),
            )
        } else {
            let preds = rev_pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[rev_best_path] {
                    predecessor = Some(*pred);
                    d = rev_dpm[*pred][j + 1][rev_best_path]
                        + scores.get(&(lnz[i], r_seq[j])).unwrap();
                    u = rev_dpm[*pred][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    l = rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', r_seq[j])).unwrap();
                }
            }
            (d, u, l)
        };

        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if lnz[i] != r_seq[j] {
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
            rev_path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            rev_path_length += 1;
        } else {
            cigar.push('L');
            j += 1;
        }
    }
    while j < dpm[0].len() - 1 {
        cigar.push('L');
        j += 1;
    }
    while i < dpm.len() - 1 {
        let predecessor = if !rev_nwp[i] {
            i + 1
        } else {
            let preds = rev_pred_hash.get_preds_and_paths(i);
            let mut p = 0;
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    p = *pred;
                }
            }
            p
        };
        cigar.push('U');
        handle_id_alignment.push(handles_nodes_id[i]);
        i = predecessor;
        rev_path_length += 1;
    }

    let mut temp_cigar = Vec::new();
    let mut temp_handle_id_alignment = Vec::new();
    i = forward_ending_node;
    j = rec_col;
    let mut forw_path_length = 0;
    while i > 0 && j > 0 {
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
            if lnz[i] != seq[j] {
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
            forw_path_length += 1;
        } else if max == u {
            temp_cigar.push('U');
            temp_handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            forw_path_length += 1;
        } else {
            temp_cigar.push('L');
            j -= 1;
        }
    }
    while j > 0 {
        temp_cigar.push('L');
        j -= 1;
    }
    while i > 0 {
        let predecessor = if !nwp[i] {
            i - 1
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let mut p = 0;
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    p = *pred;
                }
            }
            p
        };
        temp_cigar.push('U');
        temp_handle_id_alignment.push(handles_nodes_id[i]);
        i = predecessor;
        forw_path_length += 1;
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
    let (path_len, path_start, path_end) = utils::get_rec_path_len_start_end(
        handles_nodes_id,
        forward_ending_node,
        reverse_starting_node,
        0,
        rev_ending_node,
        forw_path_length,
        rev_path_length,
    );
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let recombination = if best_path == rev_best_path {
        format!("No recombination, best path: {}", best_path)
    } else {
        let fen_offset = get_node_offset(handles_nodes_id, forward_ending_node);
        let rsn_offset = get_node_offset(handles_nodes_id, reverse_starting_node);
        format!(
            "recombination path {} {}, nodes {}[{}] {}[{}]",
            best_path,
            rev_best_path,
            handles_nodes_id[forward_ending_node],
            fen_offset,
            handles_nodes_id[reverse_starting_node],
            rsn_offset
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
        path_len,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}

pub fn gaf_output_global_no_rec(
    dpm: &Vec<Vec<Vec<i32>>>,
    graph: &PathGraph,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    best_path: usize,
) -> GAFStruct {
    let pred_hash = &graph.pred_hash;
    let nwp = &graph.nwp;
    let handles_nodes_id = &graph.nodes_id_pos;
    let lnz = &graph.lnz;

    let mut path_length: usize = 0;
    let mut cigar = Vec::new();
    let mut handle_id_alignment = Vec::new();
    let mut i = 0;
    let ending_nodes = pred_hash.get_preds_and_paths(dpm.len() - 1);
    for (node, paths) in ending_nodes.iter() {
        if paths[best_path] {
            i = *node;
        }
    }
    let ending_node = i;

    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
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
            if lnz[i] != seq[j] {
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
    while i > 0 {
        let predecessor = if !nwp[i] {
            i - 1
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let mut p = 0;
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    p = *pred;
                }
            }
            p
        };
        cigar.push('U');
        handle_id_alignment.push(handles_nodes_id[i]);
        i = predecessor;
        path_length += 1;
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
    let (path_len, path_start, path_end) =
        utils::get_path_len_start_end(&handles_nodes_id, 0, ending_node, path_length);

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
        path_len,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}

pub fn gaf_output_semiglobal_rec(
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
    best_score: (f32, i32),
) -> GAFStruct {
    let mut cigar = Vec::new();
    let mut rev_path_length: usize = 0;
    let mut i = reverse_starting_node;
    let mut j = rec_col;
    let mut handle_id_alignment = Vec::new();
    let mut rev_ending_node = i;
    let mut path_sequence = Vec::new();

    // reverse alignment
    let r_seq = &get_rev_sequence(seq);
    while i > 0 && i < dpm.len() - 1 && j < dpm[0].len() - 1 {
        let mut predecessor = None;
        let (d, u, l) = if !rev_nwp[i] {
            (
                rev_dpm[i + 1][j + 1][rev_best_path] + scores.get(&(lnz[i], r_seq[j])).unwrap(),
                rev_dpm[i + 1][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap(),
                rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', r_seq[j])).unwrap(),
            )
        } else {
            let preds = rev_pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[rev_best_path] {
                    predecessor = Some(*pred);
                    d = rev_dpm[*pred][j + 1][rev_best_path]
                        + scores.get(&(lnz[i], r_seq[j])).unwrap();
                    u = rev_dpm[*pred][j][rev_best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    l = rev_dpm[i][j + 1][rev_best_path] + scores.get(&('-', r_seq[j])).unwrap();
                }
            }
            (d, u, l)
        };

        let max = *[d, u, l].iter().max().unwrap();
        rev_ending_node = i;
        if max == d {
            if lnz[i] != r_seq[j] {
                cigar.push('d');
            } else {
                cigar.push('D');
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            path_sequence.push(lnz[i]);
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            j += 1;
            rev_path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            path_sequence.push(lnz[i]);
            i = if predecessor.is_none() {
                i + 1
            } else {
                predecessor.unwrap()
            };
            rev_path_length += 1;
        } else {
            cigar.push('L');
            j += 1;
        }
    }
    while j < dpm[0].len() - 1 {
        cigar.push('L');
        j += 1;
    }

    let mut path_length: usize = 0;
    let mut temp_cigar = Vec::new();
    let mut temp_handle_id_alignment = Vec::new();
    let mut temp_path_sequence = Vec::new();

    i = forward_ending_node;
    j = rec_col;

    while i > 0 && j > 0 {
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
            if lnz[i] != seq[j] {
                temp_cigar.push('d');
            } else {
                temp_cigar.push('D');
            }
            temp_handle_id_alignment.push(handles_nodes_id[i]);
            temp_path_sequence.push(lnz[i]);
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
            temp_path_sequence.push(lnz[i]);
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

    temp_path_sequence.reverse();
    temp_path_sequence.append(&mut path_sequence);

    let path_sequence_string: String = temp_path_sequence.into_iter().collect();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    let path: Vec<usize> = temp_handle_id_alignment
        .iter()
        .map(|id| *id as usize)
        .collect();
    let start = if i == 0 { i } else { i + 1 };
    let end = rev_ending_node;
    let (path_len, path_start, path_end) = utils::get_rec_path_len_start_end(
        handles_nodes_id,
        forward_ending_node,
        reverse_starting_node,
        start,
        end,
        path_length,
        rev_path_length,
    );
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let recombination = if best_path == rev_best_path {
        format!("No recombination, best path: {}", best_path)
    } else {
        let fen_offset = get_node_offset(handles_nodes_id, forward_ending_node);
        let rsn_offset = get_node_offset(handles_nodes_id, reverse_starting_node);
        format!(
            "recombination path {} {}, nodes {}[{}] {}[{}], score: {}, displacement: {}\t{}",
            best_path,
            rev_best_path,
            handles_nodes_id[forward_ending_node],
            fen_offset,
            handles_nodes_id[reverse_starting_node],
            rsn_offset,
            best_score.0,
            best_score.1,
            path_sequence_string
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
        path_len,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    gaf_output
}

pub fn gaf_output_semiglobal_no_rec(
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
    let mut path_sequence = Vec::new();

    let score = dpm[i][j][best_path];

    while i > 0 && j > 0 {
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
            if lnz[i] == seq[j] {
                cigar.push('D')
            } else {
                cigar.push('d')
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            path_sequence.push(lnz[i]);
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
            path_sequence.push(lnz[i]);

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
    path_sequence.reverse();
    let path_sequence_string: String = path_sequence.into_iter().collect();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    let path: Vec<usize> = handle_id_alignment.iter().map(|id| *id as usize).collect();
    let (path_len, path_start, path_end) = utils::get_path_len_start_end(
        &handles_nodes_id,
        if i == 0 { i } else { i + 1 },
        ending_node,
        path_length,
    );

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = format!(
        "{}, best path: {}, score: {}\t{}",
        build_cigar(&cigar),
        best_path,
        score,
        path_sequence_string
    );
    let gaf_output = GAFStruct::build_gaf_struct(
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
    );
    gaf_output
}
