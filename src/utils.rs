use std::{
    cmp::{self, Ordering},
    collections::HashMap,
    fs::{File, OpenOptions},
};

use bit_vec::BitVec;
use handlegraph::{handle::Handle, handlegraph::HandleGraph, hashgraph::HashGraph};

use crate::args_parser;
use std::io::{prelude::*, BufWriter};
use std::path::Path;

#[inline]
/// Needed for adaptive band settings, set the leftmost and rightmost position for each row of the dp matrix   
/// The algorithm used is the same as abPOA
pub fn set_ampl_for_row(
    i: usize,
    p_arr: &[usize],
    r_val: usize,
    best_scoring_pos: &[usize],
    seq_len: usize,
    bta: usize,
    simd_version: bool,
) -> (usize, usize) {
    let ms;
    let me;
    if i == 0 {
        ms = 0;
        me = 0;
    } else if p_arr.is_empty() {
        let pl = best_scoring_pos[i - 1];
        ms = pl + 1;
        me = pl + 1;
    } else {
        let mut pl = 0;
        let mut pr = 0;
        let mut first = true;
        for p in p_arr.iter() {
            let current_best = best_scoring_pos[*p];
            if first {
                pl = current_best;
                pr = current_best;
                first = false;
            }
            if current_best < pl {
                pl = current_best;
            }
            if current_best > pr {
                pr = current_best;
            }
        }
        ms = pl + 1;
        me = pr + 1;
    }
    let tmp_bs = cmp::min(ms as i32, (seq_len as i32 - r_val as i32) - bta as i32);
    let band_start = if tmp_bs < 0 {
        0
    } else {
        cmp::max(0, tmp_bs as usize)
    };
    let band_end = if seq_len > r_val {
        cmp::min(seq_len, cmp::max(me, seq_len - r_val) + bta)
    } else {
        cmp::min(seq_len, me + bta)
    };
    if simd_version {
        set_left_right_x64(band_start, band_end, seq_len)
    } else {
        (band_start, band_end)
    }
}

fn set_left_right_x64(left: usize, right: usize, seq_len: usize) -> (usize, usize) {
    let mut new_right = right;
    let mut new_left = left;
    while (new_right - new_left) % 8 != 0 {
        if (new_right - new_left) % 2 == 0 && new_right < seq_len {
            new_right += 1;
        } else if new_left > 0 {
            new_left -= 1;
        } else {
            break;
        }
    }
    if new_left == 0 {
        while (new_right - 1) % 8 != 0 && new_right < seq_len {
            new_right += 1;
        }
    }
    if new_right == seq_len {
        while (new_right - new_left) % 8 != 0 && new_left > 1 {
            new_left -= 1
        }
    }

    (new_left, new_right)
}

/// Set R score for each node of the graph, this is done before the dp algorithm.   
/// R represent the most likely distance of each node to the last node of the graph and is used
/// in order to compute the band size for this node in the DP matrix   
pub fn set_r_values(
    nwp: &bit_vec::BitVec,
    pred_hash: &HashMap<usize, Vec<usize>>,
    lnz_len: usize,
) -> Vec<usize> {
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[lnz_len - 1] = 0;
    for p in pred_hash.get(&(lnz_len - 1)).unwrap() {
        r_values[*p] = 0;
    }
    for i in (1..lnz_len - 1).rev() {
        if r_values[i] == -1 || r_values[i] > r_values[i + 1] + 1 {
            r_values[i] = r_values[i + 1] + 1;
        }
        if nwp[i] {
            for p in pred_hash.get(&i).unwrap() {
                if r_values[*p] == -1 || r_values[*p] > r_values[i] + 1 {
                    r_values[*p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values.iter().map(|x| *x as usize).collect()
}

#[inline]
pub fn get_max_d_u_l(d: i32, u: i32, l: i32) -> (i32, char) {
    match d.cmp(&u) {
        Ordering::Less => match u.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (u, 'U'),
        },
        _ => match d.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (d, 'D'),
        },
    }
}

/// Set for each node of the LnzGraph the handle id in the .gfa file,
/// this enable the creation of the gaf output with same nodes as the original .gfa file
pub fn create_handle_pos_in_lnz(
    nwp: &BitVec,
    file_path: &str,
    amb_mode: bool,
) -> HashMap<usize, String> {
    let sorted_handles = crate::graph::get_sorted_handles(file_path, amb_mode);
    let mut curr_handle_idx = 0;
    let mut handle_of_lnz_pos = HashMap::new();
    for i in 1..nwp.len() - 1 {
        if nwp[i] {
            curr_handle_idx += 1;
        }
        handle_of_lnz_pos.insert(
            i,
            sorted_handles[(curr_handle_idx - 1) as usize]
                .id()
                .to_string(),
        );
    }
    handle_of_lnz_pos.insert(0, String::from("-1"));
    handle_of_lnz_pos
}

/// Same as create_handle_pos_in_lnz, but works with an HashGraph and a LnzGraph instead of
/// a .gfa file
pub fn handle_pos_in_lnz_from_hashgraph(
    nwp: &BitVec,
    graph: &HashGraph,
    amb_mode: bool,
) -> HashMap<usize, String> {
    let mut sorted_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_handles.sort();
    if amb_mode {
        sorted_handles.reverse();
        sorted_handles = sorted_handles
            .iter()
            .map(|h| h.flip())
            .collect::<Vec<Handle>>();
    }
    let mut curr_handle_idx = 0;
    let mut handle_of_lnz_pos = HashMap::new();
    for i in 1..nwp.len() - 1 {
        if nwp[i] {
            curr_handle_idx += 1;
        }
        handle_of_lnz_pos.insert(
            i,
            sorted_handles[(curr_handle_idx - 1) as usize]
                .id()
                .to_string(),
        );
    }
    handle_of_lnz_pos.insert(0, String::from("-1"));
    handle_of_lnz_pos
}

pub fn write_gaf(gaf_out: &str, number: usize) {
    let out_file = args_parser::get_out_file();
    if out_file == "standard output" {
        println!("{}", gaf_out)
    } else {
        let file_name = Path::new(&out_file);
        let file = if file_name.exists() && number != 1 {
            OpenOptions::new()
                .write(true)
                .append(true)
                .open(file_name)
                .unwrap()
        } else {
            File::create(file_name).expect("unable to create file")
        };

        let f = &mut BufWriter::new(&file);
        writeln!(f, "{}", gaf_out).expect("error in writing");
    }
}

pub fn get_path_len_start_end(
    handles_nodes_id: &Vec<u64>,
    start: usize,
    end: usize,
    path_len: usize,
) -> (usize, usize, usize) {
    let mut path_start = 0;
    if start > 0 {
        let first_node_id = handles_nodes_id[start];
        let mut counter = start - 1;
        while counter > 0 && handles_nodes_id[counter] == first_node_id {
            counter -= 1;
            path_start += 1;
        }
    }
    let path_end = if path_len > 0 {
        path_start + path_len - 1
    } else {
        0
    };

    let mut end_offset = 0;
    if end > 0 {
        let last_node_id = handles_nodes_id[end];
        let mut counter = end + 1;
        while counter < handles_nodes_id.len() - 1 && handles_nodes_id[counter] == last_node_id {
            counter += 1;
            end_offset += 1;
        }
    }

    let path_len = path_end + end_offset + 1;
    (path_len, path_start, path_end)
}

pub fn get_rec_path_len_start_end(
    handles_nodes_id: &Vec<u64>,
    fen: usize,
    rsn: usize,
    start: usize,
    end: usize,
    forw_path_length: usize,
    rev_path_length: usize,
) -> (usize, usize, usize) {
    //forward path info
    let mut path_start = 0;
    if start > 0 {
        let first_node_id = handles_nodes_id[start];
        let mut counter = start - 1;
        while counter > 0 && handles_nodes_id[counter] == first_node_id {
            counter -= 1;
            path_start += 1;
        }
    }

    let forw_path_end = if forw_path_length > 0 {
        path_start + forw_path_length - 1
    } else {
        0
    };

    let mut forw_end_offset = 0;
    if fen > 0 {
        let last_node_id = handles_nodes_id[fen];
        let mut counter = fen + 1;
        while counter < handles_nodes_id.len() - 1 && handles_nodes_id[counter] == last_node_id {
            counter += 1;
            forw_end_offset += 1;
        }
    }
    let forw_path_len = forw_path_end + forw_end_offset + 1;

    //reverse path info
    let mut rev_path_start = 0;
    if rsn > 0 {
        let first_node_id = handles_nodes_id[rsn];
        let mut counter = rsn - 1;
        while counter > 0 && handles_nodes_id[counter] == first_node_id {
            counter -= 1;
            rev_path_start += 1;
        }
    }

    let rev_path_end = if rev_path_length > 0 {
        rev_path_start + rev_path_length - 1
    } else {
        0
    };
    let path_end = forw_path_len + rev_path_end;
    let mut end_offset = 0;
    if end > 0 {
        let last_node_id = handles_nodes_id[end];
        let mut counter = end + 1;
        while counter < handles_nodes_id.len() - 1 && handles_nodes_id[counter] == last_node_id {
            counter += 1;
            end_offset += 1;
        }
    }
    let rev_path_len = rev_path_end + end_offset + 1;
    let path_len = forw_path_len + rev_path_len;

    (path_len, path_start, path_end)
}
