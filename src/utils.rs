use std::
    fs::{File, OpenOptions}
;
use ahash::AHashMap as HashMap;

use bit_vec::BitVec;
use handlegraph::{handle::Handle, handlegraph::HandleGraph, hashgraph::HashGraph};

use std::io::{prelude::*, BufWriter};
use std::path::Path;

use crate::dp_matrix::{DpDeltas, DpMatrix};

const SIZE: usize = 6;

#[inline]
fn char_to_index(c: u8) -> usize {
    match c {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'-' => 4,
        b'N' => 5,
        _ => panic!("Carattere non valido {}", c as char),
    }
}

#[inline]
pub fn idx(c1: u8, c2: u8) -> usize {
    char_to_index(c1) * SIZE + char_to_index(c2)
}

/// Set for each node of the LnzGraph the handle id in the .gfa file,
/// this enable the creation of the gaf output with same nodes as the original .gfa file

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

pub fn write_gaf(gaf_out: &str, number: usize, out_file: &str) {
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
    handles_nodes_id: &[u64],
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
    handles_nodes_id: &[u64],
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

#[inline]
pub fn get_abs_val(
    i: usize,
    j: usize,
    alphas: &[usize],
    path: usize,
    dpm: &DpMatrix,
    deltas: &DpDeltas,
) -> i32 {
    let is_alpha = alphas[i] == path;
    if is_alpha {
        dpm.get(i, j)
    } else {
        let update = deltas.get_val(i, j, path);
        dpm.get(i, j) + update
    }
}
