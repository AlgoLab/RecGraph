use crate::args_parser;
use std::{
    collections::HashMap,
    fs::File,
    io::{prelude::*, BufReader},
};

pub fn create_score_matrix() -> HashMap<(char, char), i32> {
    let matrix_type = args_parser::get_matrix_type();
    match matrix_type.as_str() {
        "HOXD70.mtx" => create_score_matrix_hoxd70(),
        "none" => {
            let (match_score, mismatch_score) = args_parser::get_match_mismatch();
            create_score_matrix_match_mis(match_score, mismatch_score)
        }
        _ => {
            panic!("wrong matrix type")
        }
    }
}

fn create_score_matrix_match_mis(m: i32, x: i32) -> HashMap<(char, char), i32> {
    let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
    for i in ['A', 'C', 'G', 'T', '-'].iter() {
        for j in ['A', 'C', 'G', 'T', '-'].iter() {
            if i == j {
                score_matrix.insert((*i, *j), m);
            } else {
                score_matrix.insert((*i, *j), x);
            }
        }
    }
    score_matrix
}

fn create_score_matrix_hoxd70() -> HashMap<(char, char), i32> {
    let mut matrix: Vec<Vec<String>> = Vec::new();
    let file_path = project_root::get_project_root().unwrap().join("HOXD70.mtx");

    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    for line in reader.lines().flatten() {
        let mut splitted_line: Vec<String> = Vec::new();
        for elem in line.split(' ') {
            splitted_line.push(String::from(elem));
        }
        splitted_line.retain(|x| !x.is_empty());
        matrix.push(splitted_line);
    }
    matrix[0].insert(0, String::from("X"));

    let mut matrix_score: HashMap<(char, char), i32> = HashMap::new();
    for i in 1..matrix.len() {
        for j in 1..matrix[0].len() {
            let c1 = matrix[i][0]
                .chars()
                .next()
                .expect("failed to create HOXD70");
            let c2 = matrix[0][j]
                .chars()
                .next()
                .expect("failed to create HOXD70");

            matrix_score.insert((c1, c2), matrix[i][j].parse().unwrap());
        }
    }
    for ch in ['A', 'C', 'G', 'T', 'N'].iter() {
        matrix_score.insert((*ch, '-'), -200);
        matrix_score.insert(('-', *ch), -200);
    }
    matrix_score
}
