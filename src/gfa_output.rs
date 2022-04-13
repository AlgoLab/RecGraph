use std::{collections::HashMap, fs::File};
use bitvec::prelude::*;
use bit_vec::BitVec;
use crate::bitfield_path as bf;
use std::io::{prelude::*, BufWriter};

fn create_handle_pos_in_lnz(nwp: &BitVec, file_path: &str, amb_mode: bool) -> HashMap<usize, String> {
    let sorted_handles = crate::graph::get_sorted_handles(file_path, amb_mode);
    let mut curr_handle_idx = 0;
    let mut handle_of_lnz_pos = HashMap::new();
    for i in 1..nwp.len() - 1 {
        if nwp[i] && i > 1{
            curr_handle_idx += 1;
        }
        handle_of_lnz_pos.insert(i, sorted_handles[curr_handle_idx as usize].id().to_string());
    }
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
    amb_mode: bool
) {
    let hofp = create_handle_pos_in_lnz(nwp, file_path, amb_mode);

    let mut col = last_col;
    let mut row = last_row;
    
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();
    let mut handle_id_alignment = String::new();

    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        let curr_bv = &path[row][col];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(curr_bv);

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
                sequence_align.push(sequence[col + left]);
                graph_align.push(graph[row]);
                alignment_moves.push('|');
                handle_id_alignment.push(hofp.get(&row).unwrap().chars().nth(0).unwrap());
                row = pred;
                col = j_pos - 1;
            }
            'd' => {
                sequence_align.push(sequence[col + left]);
                graph_align.push(graph[row]);
                alignment_moves.push('.');
                handle_id_alignment.push(hofp.get(&row).unwrap().chars().nth(0).unwrap());
                row = pred;
                col = j_pos - 1;
            }
            'L' => {
                graph_align.push('-');
                sequence_align.push(sequence[col + left]);
                alignment_moves.push(' ');
                handle_id_alignment.push(' ');

                col -= 1;
            }
            'U' => {
                graph_align.push(graph[row]);
                sequence_align.push('-');
                alignment_moves.push(' ');
                handle_id_alignment.push(hofp.get(&row).unwrap().chars().nth(0).unwrap());

                row = pred;
                col = j_pos;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(graph_align, sequence_align, alignment_moves, handle_id_alignment, "gfa_mk_poa");
}

fn reverse_and_write(mut graph_al: String, mut seq_al: String, mut al_moves: String, handle_align: String,align_type: &str) {
    graph_al = graph_al.chars().rev().collect();
    al_moves = al_moves.chars().rev().collect();
    seq_al = seq_al.chars().rev().collect();
    let file_name = String::from(align_type) + "_alignment.txt";

    let path = project_root::get_project_root().unwrap().join(file_name);
    let file = File::create(path).expect("unable to create file");
    let f = &mut BufWriter::new(&file);
    let mut i = 0;
    while i < graph_al.len() {
        if i + 80 < graph_al.len() {
            write!(f, "{: >80}", "").expect("unable to write");
            writeln!(f, "[{}-{}]", i, i + 80).expect("unable to write");

            write!(f, "{}", &handle_align[i..i + 80]).expect("unable to write");
            writeln!(f, "\thandle_id").expect("unable to write");

            write!(f, "{}", &graph_al[i..i + 80]).expect("unable to write");
            writeln!(f, "\tgraph").expect("unable to write");

            write!(f, "{}", &al_moves[i..i + 80]).expect("unable to write");
            writeln!(f, "\tmatc/mis").expect("unable to write");

            write!(f, "{}", &seq_al[i..i + 80]).expect("unable to write");
            writeln!(f, "\tseq").expect("unable to write");

            writeln!(f).expect("unable to write");
        } else {
            write!(f, "{: >80}", "").expect("unable to write");
            writeln!(f, "[{}-{}]", i, graph_al.len()).expect("unable to write");

            write!(f, "{}", &handle_align[i..]).expect("unable to write");
            writeln!(f, "\thandle_id").expect("unable to write");

            write!(f, "{}", &graph_al[i..]).expect("unable to write");
            writeln!(f, "\tgraph").expect("unable to write");

            write!(f, "{}", &al_moves[i..]).expect("unable to write");
            writeln!(f, "\tmatc/mis").expect("unable to write");

            write!(f, "{}", &seq_al[i..]).expect("unable to write");
            writeln!(f, "\tseq").expect("unable to write");

            writeln!(f).expect("unable to write");
        }
        i += 80;
    }
}
