use crate::bitfield_path as bf;
use bit_vec::BitVec;
use bitvec::prelude::*;
use std::io::{prelude::*, BufWriter};
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

pub fn gfa_of_abpoa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    graph: &[char],
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

    let mut ref_nodes = Vec::new();
    let mut ref_node = String::new();

    let mut read_nodes = Vec::new();
    let mut read_node = String::new();

    let mut count_M = 0;
    let mut count_I = 0;
    let mut count_D = 0;

    let mut curr_handle = "";
    let mut last_dir = ' ';
    let mut path_length = 0;
    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        let curr_bv = &path[row][col];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(curr_bv);

        if hofp.get(&row).unwrap() != curr_handle {
            cigar = set_cigar_substring(count_M, count_I, count_D, cigar);
            cigars.insert(0, cigar);

            cigar = String::new();
            read_node = format!("{}\t{}", curr_handle, read_node);
            ref_node = format!("{}\t{}", curr_handle, ref_node);
            read_nodes.insert(0, read_node);
            ref_nodes.insert(0, ref_node);

            read_node = String::new();
            ref_node = String::new();
            count_M = 0;
            count_I = 0;
            count_D = 0;
        }
        curr_handle = hofp.get(&row).unwrap();
        if dir.to_ascii_uppercase() != last_dir.to_ascii_uppercase() {
            cigar = set_cigar_substring(count_M, count_I, count_D, cigar);
            count_M = 0;
            count_I = 0;
            count_D = 0;
        }
        last_dir = dir;

        let left = ampl_for_row[row].0;
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
                ref_node.insert(0, graph[row]);
                read_node.insert(0, sequence[col + left]);

                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col = j_pos - 1;
                count_M += 1;
                path_length += 1;
            }
            'd' => {
                ref_node.insert(0, graph[row]);
                read_node.insert(0, sequence[col + left]);

                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col = j_pos - 1;
                count_M += 1;
                path_length += 1;
            }
            'L' => {
                read_node.insert(0, sequence[col + left]);

                col -= 1;
                count_I += 1;
            }
            'U' => {
                ref_node.insert(0, graph[row]);
                handle_id_alignment.push(hofp.get(&row).unwrap());

                row = pred;
                col = j_pos;
                count_D += 1;
                path_length += 1;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    cigar = set_cigar_substring(count_M, count_I, count_D, cigar);
    cigars.insert(0, cigar);

    read_node = format!("{}\t{}", curr_handle, read_node);
    read_nodes.insert(0, read_node);

    ref_node = format!("{}\t{}", curr_handle, ref_node);
    ref_nodes.insert(0, ref_node);

    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    write_alignment(
        ref_nodes,
        read_nodes,
        &cigars,
        &handle_id_alignment,
        "gfa_mk_poa",
    );
    let seq_name = "query_name"; // to set
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
    let path_start = "*"; // to set 
    let path_end = "*"; // to set
    let number_residue = "*"; // to set
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars.join("\t");
    let gaf_out = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
    let fasta_name = "fasta_name"; //to set
    let file_name = String::from(fasta_name) + "_gaf.txt";
    let path = project_root::get_project_root().unwrap().join(file_name);
    let file = File::create(path).expect("unable to create file");
    let f = &mut BufWriter::new(&file);
    writeln!(f,"{}", gaf_out).expect("error in writing");

}

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

fn set_cigar_substring(count_m: i32, count_i: i32, count_d: i32, cs: String) -> String {
    let cigar;
    if count_m > 0 {
        cigar = format!("{}M{}", count_m, cs);
    } else if count_i > 0 {
        cigar = format!("{}I{}", count_i, cs);
    } else if count_d > 0 {
        cigar = format!("{}D{}", count_d, cs);
    } else {
        cigar = format!("{}", cs);
    };
    cigar
}
