use crate::bitfield_path as bf;
use bitvec::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufWriter};

pub fn write_align_banded_poa(
    path: &[Vec<bitvec::prelude::BitVec<u16, Msb0>>],
    sequence: &[char],
    graph: &[char],
    ampl_for_row: &[(usize, usize)],
    last_row: usize,
    last_col: usize,
) {
    let mut col = last_col;
    let mut row = last_row;
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();

    while bf::dir_from_bitvec(&path[row][col]) != 'O' {
        let curr_bv = &path[row][col];
        let pred = bf::pred_from_bitvec(curr_bv);
        let dir = bf::dir_from_bitvec(&curr_bv);

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
                row = pred;
                col = j_pos - 1;
            }
            'd' => {
                sequence_align.push(sequence[col + left]);
                graph_align.push(graph[row]);
                alignment_moves.push('.');
                row = pred;
                col = j_pos - 1;
            }
            'L' => {
                graph_align.push('-');
                sequence_align.push(sequence[col + left]);
                alignment_moves.push(' ');
                col -= 1;
            }
            'U' => {
                graph_align.push(graph[row]);
                sequence_align.push('-');
                alignment_moves.push(' ');
                row = pred;
                col = j_pos;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(graph_align, sequence_align, alignment_moves, "mk_poa");
}

pub fn write_align_local_poa(
    path: &[Vec<(char, usize)>],
    sequence: &[char],
    graph: &[char],
    best_row: usize,
    best_col: usize,
) {
    let mut col = best_col;
    let mut row = best_row;
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();
    while path[row][col] != ('O', 0) {
        match path[row][col] {
            ('D', _) => {
                sequence_align.push(sequence[col]);
                graph_align.push(graph[row]);
                alignment_moves.push('|');
                row = path[row][col].1;
                col -= 1;
            }
            ('d', _) => {
                sequence_align.push(sequence[col]);
                graph_align.push(graph[row]);
                alignment_moves.push('.');
                row = path[row][col].1;
                col -= 1;
            }
            ('L', _) => {
                graph_align.push('-');
                sequence_align.push(sequence[col]);
                alignment_moves.push(' ');
                col -= 1;
            }
            ('U', _) => {
                graph_align.push(graph[row]);
                sequence_align.push('-');
                alignment_moves.push(' ');
                row = path[row][col].1;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(graph_align, sequence_align, alignment_moves, "local_poa");
}

pub fn write_align_gap_local_poa(
    path: &[Vec<(char, usize)>],
    path_x: &[Vec<(char, usize)>],
    path_y: &[Vec<(char, usize)>],
    sequence: &[char],
    graph: &[char],
    best_row: usize,
    best_col: usize,
) {
    let mut col = best_col;
    let mut row = best_row;
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();
    while path[row][col] != ('O', 0) {
        match path[row][col] {
            ('D', _) => {
                sequence_align.push(sequence[col]);
                graph_align.push(graph[row]);
                alignment_moves.push('|');
                row = path[row][col].1;
                col -= 1;
            }
            ('d', _) => {
                sequence_align.push(sequence[col]);
                graph_align.push(graph[row]);
                alignment_moves.push('.');
                row = path[row][col].1;
                col -= 1;
            }
            ('L', _) => {
                if path_x[row][col].0 == 'X' {
                    while path_x[row][col].0 == 'X' {
                        graph_align.push('-');
                        sequence_align.push(sequence[col]);
                        alignment_moves.push(' ');
                        col -= 1;
                    }
                } else {
                    graph_align.push('-');
                    sequence_align.push(sequence[col]);
                    alignment_moves.push(' ');
                    col -= 1;
                }
            }
            ('U', _) => {
                if path_y[row][col].0 == 'Y' {
                    while path_y[row][col].0 == 'Y' {
                        graph_align.push(graph[row]);
                        sequence_align.push('-');
                        alignment_moves.push(' ');
                        row = path_y[row][col].1;
                    }
                } else {
                    graph_align.push(graph[row]);
                    sequence_align.push('-');
                    alignment_moves.push(' ');
                    row = path_y[row][col].1;
                }
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(
        graph_align,
        sequence_align,
        alignment_moves,
        "gap_local_poa",
    );
}

pub fn write_align_gap_mk_abpoa(
    path: &[Vec<(char, usize)>],
    path_x: &[Vec<(char, usize)>],
    path_y: &[Vec<(char, usize)>],
    ampl_for_row: &[(usize, usize)],
    sequence: &[char],
    graph: &[char],
    best_row: usize,
    best_col: usize,
    pred_hash: &HashMap<usize, Vec<usize>>,
) {
    let mut col = best_col;
    let mut row = best_row;
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();
    while path[row][col] != ('O', 0) {
        let left = ampl_for_row[row].0;
        match path[row][col] {
            ('D', _) => {
                sequence_align.push(sequence[col + left]);
                graph_align.push(graph[row]);
                alignment_moves.push('|');
                let p = path[row][col].1;
                let left_p = ampl_for_row[p].0;
                let j_pos = if left_p < left {
                    col + (left - left_p)
                } else {
                    col - (left_p - left)
                };
                col = j_pos - 1;
                row = p;
            }
            ('d', _) => {
                sequence_align.push(sequence[col + left]);
                graph_align.push(graph[row]);
                alignment_moves.push('.');
                let p = path[row][col].1;
                let left_p = ampl_for_row[p].0;
                let j_pos = if left_p < left {
                    col + (left - left_p)
                } else {
                    col - (left_p - left)
                };
                col = j_pos - 1;
                row = p;
            }
            ('L', _) => {
                if path_x[row][col].0 == 'X' {
                    while path_x[row][col].0 == 'X' && col > 0 {
                        graph_align.push('-');
                        sequence_align.push(sequence[col + left]);
                        alignment_moves.push(' ');
                        col -= 1
                    }
                } else {
                    graph_align.push('-');
                    sequence_align.push(sequence[col + left]);
                    alignment_moves.push(' ');
                    col -= 1
                }
            }
            ('U', _) => {
                if path_y[row][col].0 == 'Y' {
                    while path_y[row][col].0 == 'Y' {
                        let left_row = ampl_for_row[row].0;
                        graph_align.push(graph[row]);
                        sequence_align.push('-');
                        alignment_moves.push(' ');
                        let p = path_y[row][col].1;
                        let left_p = ampl_for_row[p].0;
                        let j_pos = if left_p < left_row {
                            col + (left_row - left_p)
                        } else {
                            col - (left_p - left_row)
                        };
                        col = j_pos;
                        row = p;
                    }
                } else {
                    graph_align.push(graph[row]);
                    sequence_align.push('-');
                    alignment_moves.push(' ');
                    let p = path[row][col].1;
                    let left_p = ampl_for_row[p].0;
                    let j_pos = if left_p < left {
                        col + (left - left_p)
                    } else {
                        col - (left_p - left)
                    };
                    col = j_pos;
                    row = p;
                }
            }
            ('u', _) => {
                println!("u");
                col = 0;
                while row > 0 {
                    graph_align.push(graph[row]);
                    sequence_align.push('-');
                    alignment_moves.push(' ');
                    match pred_hash.get(&row) {
                        Some(p_arr) => {
                            row = *p_arr.iter().min().unwrap();
                        }
                        _ => {
                            row -= 1;
                        }
                    }
                }
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(graph_align, sequence_align, alignment_moves, "gap_mk_abpoa");
}

fn reverse_and_write(mut s1_al: String, mut s2_al: String, mut al_moves: String, align_type: &str) {
    s1_al = s1_al.chars().rev().collect();
    al_moves = al_moves.chars().rev().collect();
    s2_al = s2_al.chars().rev().collect();
    let file_name = String::from(align_type) + "_alignment.txt";

    let path = project_root::get_project_root().unwrap().join(file_name);
    let file = File::create(path).expect("unable to create file");
    let f = &mut BufWriter::new(&file);
    let mut i = 0;
    while i < s1_al.len() {
        if i + 80 < s1_al.len() {
            write!(f, "{: >80}", "").expect("unable to write");
            writeln!(f, "[{}-{}]", i, i + 80).expect("unable to write");

            write!(f, "{}", &s1_al[i..i + 80]).expect("unable to write");
            writeln!(f, "\tgraph").expect("unable to write");

            write!(f, "{}", &al_moves[i..i + 80]).expect("unable to write");
            writeln!(f, "\tmatc/mis").expect("unable to write");

            write!(f, "{}", &s2_al[i..i + 80]).expect("unable to write");
            writeln!(f, "\tseq").expect("unable to write");

            writeln!(f).expect("unable to write");
        } else {
            write!(f, "{: >80}", "").expect("unable to write");
            writeln!(f, "[{}-{}]", i, s1_al.len()).expect("unable to write");

            write!(f, "{}", &s1_al[i..]).expect("unable to write");
            writeln!(f, "\tgraph").expect("unable to write");

            write!(f, "{}", &al_moves[i..]).expect("unable to write");
            writeln!(f, "\tmatc/mis").expect("unable to write");

            write!(f, "{}", &s2_al[i..]).expect("unable to write");
            writeln!(f, "\tseq").expect("unable to write");

            writeln!(f).expect("unable to write");
        }
        i += 80;
    }
    drop(f);
}
