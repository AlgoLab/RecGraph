mod ab_mk_edit_distance;
use ab_mk_edit_distance as abmked;

mod ab_global_alignement;
mod args_parser;
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

    let ampl = match s1.len() < s2.len() {
        true => s2.len() - s1.len(),
        _ => s1.len() - s2.len(),
    };

    let score_matrix = matrix::create_score_matrix();

    //edit distance con banda su matrice m*k
    abmked::ab_ed_km(&s1, &s2, cmp::max((ampl * 2 + 1) as i32, 3));

    //glob alignment
    ab_global_alignement::ab_glob_alignement(
        &s1,
        &s2,
        &score_matrix,
        cmp::max((ampl * 2 + 1) as i32, 3),
    );
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
