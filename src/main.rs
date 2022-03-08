mod ab_mk_edit_distance;
use ab_local_alignment::ab_loc_alignment;
use ab_mk_edit_distance as abmked;

mod ab_global_alignment;
mod ab_local_alignment;
mod args_parser;
mod basic_output;
mod matrix;

use std::cmp;
use std::fs::File;
use std::io::{prelude::*, BufReader};

fn main() {
    let sequences = get_sequences();

    let mut s1: Vec<char> = sequences[3].chars().collect();
    let mut s2: Vec<char> = sequences[4].chars().collect();
    s1.insert(0, '$');
    s2.insert(0, '$');

    let ampl = match s1.len() < s2.len() {
        true => s2.len() - s1.len(),
        _ => s1.len() - s2.len(),
    };

    let score_matrix = matrix::create_score_matrix();
    // TODO gestione tipologia alignment da linea di comando
    //edit distance con banda su matrice m*k
    abmked::ab_ed_km(&s1, &s2, cmp::max((ampl * 2 + 1) as i32, 3));

    //glob alignment
    ab_global_alignment::ab_glob_alignement(
        &s1,
        &s2,
        &score_matrix,
        cmp::max((ampl * 2 + 1) as i32, 3),
    );

    //local alignment
    ab_loc_alignment(&s1, &s2, &score_matrix, cmp::max((ampl * 2 + 1) as i32, 3))
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
