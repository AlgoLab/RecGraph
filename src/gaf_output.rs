use crate::args_parser;
use crate::bitfield_path as bf;
use bitvec::prelude::*;
use std::fs::OpenOptions;
use std::io::{prelude::*, BufWriter};
use std::path::Path;
use std::{collections::HashMap, fs::File};
// TODO: gaf_out diventa struct GAFStruct, con metodo to_string che ritorna una stringa = a gaf_out originale per write_output
// mettere flag per stabilire se scrivere output oppure ritornare GAFStruct, per farlo il tipo di ritorno deve essere
// Option(GafStruct), None se voglio scrivere output, Some se voglio struct.
pub struct GAFStruct {
    query_name: String,
    query_length: usize,
    query_start: usize,
    query_end: usize,
    strand: char,
    path: String,
    path_length: usize,
    path_start: usize,
    path_end: usize,
    residue_matches_number: usize,
    alignment_block_length: String,
    mapping_quality: String,
    comments: String,
}
impl GAFStruct {
    pub fn new(
        query_name: String,
        query_length: usize,
        query_start: usize,
        query_end: usize,
        strand: char,
        path: String,
        path_length: usize,
        path_start: usize,
        path_end: usize,
        residue_matches_number: usize,
        alignment_block_length: String,
        mapping_quality: String,
        comments: String,
    ) -> GAFStruct {
        GAFStruct {
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            path,
            path_length,
            path_start,
            path_end,
            residue_matches_number,
            alignment_block_length,
            mapping_quality,
            comments,
        }
    }
    pub fn to_string(self) -> String {
        let gaf_struct_to_string = format!(
            "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            self.path,
            self.path_length,
            self.path_start,
            self.path_end,
            self.residue_matches_number,
            self.alignment_block_length,
            self.mapping_quality,
            self.comments
        );
        gaf_struct_to_string
    }
}
pub fn gaf_of_gap_abpoa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_x: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_y: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: (&str, usize),
    ampl_for_row: &[(usize, usize)],
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
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
                residue_matching += 1;
            }
            'd' => {
                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col = j_pos - 1;
                count_m += 1;
                path_length += 1;
            }
            'L' => {
                if bf::dir_from_bitvec(&path_x[row][col]) == 'X' {
                    while bf::dir_from_bitvec(&path_x[row][col]) == 'X' {
                        count_d += 1;
                        col -= 1;
                    }
                } else {
                    count_d += 1;
                    col -= 1;
                }
            }
            'U' => {
                if bf::dir_from_bitvec(&path_y[row][col]) == 'Y' {
                    while bf::dir_from_bitvec(&path_y[row][col]) == 'Y' {
                        let left_row = ampl_for_row[row].0;
                        let p = bf::pred_from_bitvec(&path_y[row][col]);
                        let left_p = ampl_for_row[p].0;
                        let j_pos = if left_p < left_row {
                            col + (left_row - left_p)
                        } else {
                            col - (left_p - left_row)
                        };
                        handle_id_alignment.push(hofp.get(&row).unwrap());
                        count_i += 1;
                        path_length += 1;
                        col = j_pos;
                        row = p;
                    }
                } else {
                    handle_id_alignment.push(hofp.get(&row).unwrap());
                    count_i += 1;
                    path_length += 1;
                    row = pred;
                    col = j_pos;
                }
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
    let strand = if amb_mode { '-' } else { '+' };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(hofp, row); // first letter used in first node of alignment
    let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
    let number_residue = residue_matching;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");

    let gaf_struct = GAFStruct::new(
        String::from(seq_name.0),
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    let gaf_out = gaf_struct.to_string();
    write_gaf(&gaf_out, seq_name.1);
}
pub fn gaf_of_global_abpoa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: (&str, usize),
    ampl_for_row: &[(usize, usize)],
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
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
                residue_matching += 1;
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
    let strand = if amb_mode { '-' } else { '+' };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(hofp, row); // first letter used in first node of alignment
    let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
    let number_residue = residue_matching; // to set
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_struct = GAFStruct::new(
        String::from(seq_name.0),
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    let gaf_out = gaf_struct.to_string();
    write_gaf(&gaf_out, seq_name.1);
}
pub fn gaf_of_local_poa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: (&str, usize),
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
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
    let strand = if amb_mode { '-' } else { '+' };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(hofp, row); // first letter used in first node of alignment
    let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
    let number_residue_matching = residue_matching;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_struct = GAFStruct::new(
        String::from(seq_name.0),
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue_matching,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    let gaf_out = gaf_struct.to_string();
    write_gaf(&gaf_out, seq_name.1);
}

pub fn gaf_of_gap_local_poa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_x: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    path_y: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    seq_name: (&str, usize),
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
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
            }
            'L' => {
                if bf::dir_from_bitvec(&path_x[row][col]) == 'X' {
                    while bf::dir_from_bitvec(&path_x[row][col]) == 'X' {
                        count_d += 1;
                        col -= 1;
                    }
                } else {
                    count_d += 1;
                    col -= 1;
                }
            }
            'U' => {
                if bf::dir_from_bitvec(&path_y[row][col]) == 'Y' {
                    while bf::dir_from_bitvec(&path_y[row][col]) == 'Y' {
                        let p = bf::pred_from_bitvec(&path_y[row][col]);
                        handle_id_alignment.push(hofp.get(&row).unwrap());
                        row = p;
                        count_i += 1;
                        path_length += 1;
                    }
                } else {
                    handle_id_alignment.push(hofp.get(&row).unwrap());
                    count_i += 1;
                    path_length += 1;
                    row = pred;
                }
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
    let strand = if amb_mode { '-' } else { '+' };
    let path_matching: String = handle_id_alignment
        .iter()
        .map(|line| line.chars().collect::<Vec<char>>().into_iter().collect())
        .collect::<Vec<String>>()
        .join(">");
    //path_length obtained from iterating in path matrix
    let path_start = node_start(hofp, row); // first letter used in first node of alignment
    let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
    let number_residue_matching = residue_matching;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_struct = GAFStruct::new(
        String::from(seq_name.0),
        seq_length,
        query_start,
        query_end,
        strand,
        path_matching,
        path_length,
        path_start,
        path_end,
        number_residue_matching,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );

    let gaf_out = gaf_struct.to_string();
    write_gaf(&gaf_out, seq_name.1);
}

pub fn gaf_of_local_poa_simd(
    path: &[Vec<f32>],
    sequence: &[char],
    seq_name: (&str, usize),
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
    let mut col = last_col;
    let mut row = last_row;

    let mut handle_id_alignment = Vec::new();

    let mut cigars = Vec::new();
    let mut cigar = String::new();

    let mut count_m = 0;
    let mut count_i = 0;
    let mut count_d = 0;

    let mut curr_handle = "";
    let mut last_dir = -1;
    let mut path_length = 0;
    let mut residue_matching = 0;
    while path[row][col] != 0.0 {
        let val = path[row][col];

        let val_str = val.to_string();
        let pred_dir = val_str.split('.').collect::<Vec<&str>>();
        let pred = pred_dir[0].parse::<usize>().unwrap();
        let dir = pred_dir[1].parse::<i32>().unwrap();

        if hofp.get(&row).unwrap() != curr_handle {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            cigars.insert(0, cigar);

            cigar = String::new();
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        curr_handle = hofp.get(&row).unwrap();
        if dir != last_dir {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        last_dir = dir;

        match dir {
            1 => {
                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col -= 1;
                count_m += 1;
                path_length += 1;
                residue_matching += 1;
            }
            3 => {
                col -= 1;
                count_d += 1;
            }
            2 => {
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
    let path_start = node_start(hofp, row); // first letter used in first node of alignment
    let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
    let number_residue_matching = residue_matching;
    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = cigars[..cigars.len() - 1].join(",");
    let gaf_out = format!(
        "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        seq_name.0,
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
    write_gaf(&gaf_out, seq_name.1);
}

pub fn gaf_of_global_abpoa_simd(
    path: &[Vec<f32>],
    sequence: &[char],
    seq_name: (&str, usize),
    last_row: usize,
    last_col: usize,
    amb_mode: bool,
    hofp: &HashMap<usize, String>,
) {
    let mut col = last_col;
    let mut row = last_row;

    let mut handle_id_alignment = Vec::new();

    let mut cigars = Vec::new();
    let mut cigar = String::new();

    let mut count_m = 0;
    let mut count_i = 0;
    let mut count_d = 0;

    let mut curr_handle = "";
    let mut last_dir = -1;
    let mut path_length = 0;
    let mut residue_matching = 0;
    let mut out_ok = true;
    while path[row][col] != 0.0 {
        let val = path[row][col];
        if val == -1f32 {
            out_ok = false;
            break;
        }
        let val_str = val.to_string();
        let pred_dir = val_str.split('.').collect::<Vec<&str>>();
        let pred = pred_dir[0].parse::<usize>().unwrap();
        let dir = pred_dir[1].parse::<i32>().unwrap();

        if hofp.get(&row).unwrap() != curr_handle {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            cigars.insert(0, cigar);

            cigar = String::new();

            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        curr_handle = hofp.get(&row).unwrap();
        if dir != last_dir {
            cigar = set_cigar_substring(count_m, count_i, count_d, cigar);
            count_m = 0;
            count_i = 0;
            count_d = 0;
        }
        last_dir = dir;
        match dir {
            1 => {
                handle_id_alignment.push(hofp.get(&row).unwrap());
                row = pred;
                col -= 1;
                count_m += 1;
                path_length += 1;
                residue_matching += 1;
            }
            3 => {
                col -= 1;
                count_d += 1;
            }
            2 => {
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
    if out_ok {
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
        let path_start = node_start(hofp, row); // first letter used in first node of alignment
        let path_end = node_start(hofp, last_row); // last letter used in last node of alignment
        let number_residue = residue_matching; // to set
        let align_block_length = "*"; // to set
        let mapping_quality = "*"; // to set
        let comments = cigars[..cigars.len() - 1].join(",");
        let gaf_out = format!(
            "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            seq_name.0,
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
        write_gaf(&gaf_out, seq_name.1);
    } else {
        println!("band not enough for correct output");
    }
}

fn node_start(hofp: &HashMap<usize, String>, row: usize) -> usize {
    let handle_id = hofp.get(&row).unwrap();
    let mut i = row;
    while hofp.get(&i).unwrap() == handle_id && i > 0 {
        i -= 1;
    }
    row - i
}

fn write_gaf(gaf_out: &str, number: usize) {
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

fn set_cigar_substring(count_m: i32, count_i: i32, count_d: i32, cs: String) -> String {
    if (count_m * count_i) + (count_i * count_d) + (count_m * count_d) != 0 {
        panic!("wrong format in cigar string")
    }

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
