mod ab_gap_global_alignment;
mod ab_global_alignment;
mod ab_mk_edit_distance;
mod args_parser;
mod basic_output;
mod global_alignment_affine_gap;
mod local_alignment;
mod matrix;

use std::cmp;
use std::fs::File;
use std::io::{prelude::*, BufReader};

fn main() {
    let sequences = get_sequences();

    let mut s1: Vec<char> = sequences[0].chars().collect();
    let mut s2: Vec<char> = sequences[1].chars().collect();
    s1.insert(0, '$');
    s2.insert(0, '$');

    let score_matrix = matrix::create_score_matrix();
    let ampl = match s1.len() < s2.len() {
        true => s2.len() - s1.len(),
        _ => s1.len() - s2.len(),
    };

    let align_mode = args_parser::get_align_mode();
    match align_mode {
        0 => ab_global_alignment::exec(&s1, &s2, &score_matrix, cmp::max(ampl * 2 + 1, 3)),
        1 => {
            /*
            if ampl * 2 + 1 > cmp::min(s1.len(), s2.len()) {
                local_alignment::exec(&s1, &s2, &score_matrix)
            } else {
                ab_local_alignment::exec(&s1, &s2, &score_matrix, cmp::max(ampl * 2 + 1, 3))
            }
            */
            local_alignment::exec(&s1, &s2, &score_matrix)
        }
        2 => ab_mk_edit_distance::exec(&s1, &s2, cmp::max((ampl * 2 + 1) as i32, 3)),
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            global_alignment_affine_gap::exec(&s1, &s2, &score_matrix, g_open, g_ext);
            ab_gap_global_alignment::exec(
                &s1,
                &s2,
                &score_matrix,
                cmp::max(ampl * 2 + 1, 3),
                g_open,
                g_ext,
            );
        }
        _ => panic!("alignment mode must be 0, 1, 2 or 3"),
    }
}

fn get_sequences() -> Vec<String> {
    let mut sequences: Vec<String> = Vec::new();
    let file_path = project_root::get_project_root()
        .unwrap()
        .join("sequences.txt");

    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    for line in reader.lines().flatten() {
        sequences.push(line);
    }
    sequences
}
