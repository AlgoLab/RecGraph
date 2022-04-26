use crate::args_parser;
use crate::bitfield_path as bf;
use bit_vec::BitVec;
use bitvec::prelude::*;
use std::fs::OpenOptions;
use std::io::{prelude::*, BufWriter};
use std::path::Path;
use std::{collections::HashMap, fs::File};
fn create_handle_pos_in_lnz(
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

pub fn gaf_of_global_abpoa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: &str,
    //graph: &[char],  needed for path start and end?
    ampl_for_row: &[(usize, usize)],
    last_row: usize,
    last_col: usize,
    nwp: &BitVec,
    file_path: &str,
    amb_mode: bool,
) {
    let hofp = create_handle_pos_in_lnz(nwp, file_path, amb_mode);
    let mut col = last_col;
    let mut row = last_row;

    let mut handle_id_alignment = Vec::new();

    let mut cigars = Vec::new();
    let mut cigar = String::new();

    let mut count_m = 0;
    let mut count_i = 0;
    let mut count_d = 0;

    let mut curr_handle = "";
    let mut last_dir = ' ';
    let mut path_length = 0;
    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        let curr_bv = &path[row][col];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(curr_bv);

        if hofp.get(&row).unwrap() != curr_handle {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            cigars.insert(0, cigar);

            cigar = String::new();
            
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        curr_handle = hofp.get(&row).unwrap();
        if dir.to_ascii_uppercase() != last_dir.to_ascii_uppercase() {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        last_dir = dir;

        let p_left = ampl_for_row[pred].0;
        let j_pos = if ampl_for_row[row].0 < p_left {
            let delta = p_left - ampl_for_row[row].0;
            col - delta
        } else {
            let delta = ampl_for_row[row].0 - p_left;
            col + delta
        };

        match dir {
            'D' => {

                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col = j_pos - 1;
                count_m += 1;
                path_length += 1;
            }
            'd' => {

                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col = j_pos - 1;
                count_m += 1;
                path_length += 1;
            }
            'L' => {

                col -= 1;
                count_d += 1;
            }
            'U' => {
                handle_id_alignment.push(hofp.get(&row).unwrap());

                row = pred;
                col = j_pos;
                count_i += 1;
                path_length += 1;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
    cigars.insert(0, cigar);

    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    
    let seq_length = sequence.len() - 1; // $ doesn't count
    let query_start = col;
    let query_end = last_col + ampl_for_row.get(last_row).unwrap().0;
    let strand = if amb_mode { "-" } else { "+" };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(&hofp, row); // first letter used in first node of alignment
    let path_end = node_start(&hofp, last_row); // last letter used in last node of alignment
    let number_residue = "*"; // to set
    let align_block_length = path_length;
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_out = format!(
        "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        seq_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue,
        align_block_length,
        mapping_quality,
        comments
    );
    write_gaf(&gaf_out);
}
pub fn gaf_of_local_poa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: &str,
    //graph: &[char],  needed for path start and end?
    last_row: usize,
    last_col: usize,
    nwp: &BitVec,
    file_path: &str,
    amb_mode: bool,
) {
    let hofp = create_handle_pos_in_lnz(nwp, file_path, amb_mode);
    let mut col = last_col;
    let mut row = last_row;

    let mut handle_id_alignment = Vec::new();

    let mut cigars = Vec::new();
    let mut cigar = String::new();

    let mut count_m = 0;
    let mut count_i = 0;
    let mut count_d = 0;

    let mut curr_handle = "";
    let mut last_dir = ' ';
    let mut path_length = 0;
    let mut residue_matching = 0;
    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        let curr_bv = &path[row][col];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(curr_bv);

        if hofp.get(&row).unwrap() != curr_handle {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            cigars.insert(0, cigar);

            cigar = String::new();
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        curr_handle = hofp.get(&row).unwrap();
        if dir.to_ascii_uppercase() != last_dir.to_ascii_uppercase() {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        last_dir = dir;

        match dir {
            'D' => {
                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col -= 1;
                count_m += 1;
                path_length += 1;
                residue_matching += 1;
            }
            'd' => {
                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col -= 1;
                count_m += 1;
                path_length += 1;
                residue_matching += 1;
            }
            'L' => {
                col -= 1;
                count_d += 1;
            }
            'U' => {
                handle_id_alignment.push(hofp.get(&row).unwrap());

                row = pred;
                count_i += 1;
                path_length += 1;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
    cigars.insert(0, cigar);

    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    
    let seq_length = sequence.len() - 1; // $ doesn't count
    let query_start = col;
    let query_end = last_col;
    let strand = if amb_mode { "-" } else { "+" };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(&hofp, row); // first letter used in first node of alignment
    let path_end = node_start(&hofp, last_row); // last letter used in last node of alignment
    let number_residue_matching = residue_matching;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_out = format!(
        "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        seq_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue_matching,
        align_block_length,
        mapping_quality,
        comments
    );
    write_gaf(&gaf_out);
}
fn node_start(hofp: &HashMap<usize, String>, row: usize) -> usize {
    let handle_id = hofp.get(&row).unwrap();
    let mut i = row;
    while hofp.get(&i).unwrap() == handle_id && i > 0 {
        i -= 1;
    }
    row - i
}

fn write_gaf(gaf_out: &String) {
    let file_path = args_parser::get_graph_path();
    let file_name = Path::new(&file_path)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .split('.')
        .collect::<Vec<&str>>()[0];
    let file_name_out = String::from(file_name) + ".gaf";
    let path = project_root::get_project_root()
        .unwrap()
        .join(file_name_out);
    let file = if Path::new(&path).exists() {
        OpenOptions::new()
            .write(true)
            .append(true)
            .open(path)
            .unwrap()
    } else {
        File::create(path).expect("unable to create file")
    };

    let f = &mut BufWriter::new(&file);
    writeln!(f, "{}", gaf_out).expect("error in writing");
}
/*
fn write_alignment(
    ref_nodes: Vec<String>,
    read_nodes: Vec<String>,
    cigars: &Vec<String>,
    handle_align: &Vec<&String>,
    align_type: &str,
) {
    let file_name = String::from(align_type) + "_alignment.txt";

    let path = project_root::get_project_root().unwrap().join(file_name);
    let file = File::create(path).expect("unable to create file");
    let f = &mut BufWriter::new(&file);

    let handle_align: Vec<String> = handle_align
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect();
    let handle_ids = handle_align.join(",");

    writeln!(f, "{}", handle_ids).expect("unable to write");
    writeln!(f).expect("unable to write");

    for i in 0..ref_nodes.len() {
        writeln!(f, "{}", ref_nodes[i]).expect("unable to write");
        writeln!(f, "{}\t\t{}", read_nodes[i], cigars[i]).expect("unable to write");
        writeln!(f).expect("unable to write");
    }
}
 */
fn set_cigar_substring(count_m: i32, count_i: i32, count_d: i32, cs: String) -> String {
    let cigar;
    if count_m > 0 {
        cigar = format!("{}M{}", count_m, cs);
    } else if count_i > 0 {
        cigar = format!("{}I{}", count_i, cs);
    } else if count_d > 0 {
        cigar = format!("{}D{}", count_d, cs);
    } else {
        cigar = cs;
    };
    cigar
}
