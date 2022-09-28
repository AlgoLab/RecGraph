use bit_vec::BitVec;

use crate::{
    gaf_output::GAFStruct, pathwise_alignment_output::build_cigar,
    pathwise_alignment_recombination::get_node_offset, pathwise_graph::PredHash,
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
    let mut path_length: usize = 0;
    let mut i = reverse_starting_node;
    let mut j = rec_col + 1;
    let mut handle_id_alignment = Vec::new();
    // reverse alignment
    let mut rev_ending_node = i;
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
            if handles_nodes_id[i] != 0 {
                rev_ending_node = i;
            }
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
            if handles_nodes_id[i] != 0 {
                rev_ending_node = i;
            }
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
        path_length += 1;
        if handles_nodes_id[i] != 0 {
            rev_ending_node = i;
        }
    }

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
    let path_start = 0;
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
