use crate::args_parser;
use std::{
    collections::HashMap,
    fs::File,
    io::{prelude::*, BufReader},
};

/// Returns the same matrix of create_score_matrix with every score transformed into an f32
/// this is done in order to works with f32 values needed in alignment with simd instruction.
pub fn create_f32_scores_matrix() -> HashMap<(char, char), f32> {
    let matrix = create_score_matrix();
    let mut f32_matrix: HashMap<(char, char), f32> = HashMap::new();
    for (k, v) in matrix.iter() {
        f32_matrix.insert(*k, *v as f32);
    }
    f32_matrix
}

/// Returned score matrix can be set by match/mismatch score or by a .mtx file (currently only HOXD70 and HOXD55).
/// This function is meant to be used by recgraph directly, if you want to create a score matrix use the functions defined in api.rs
pub fn create_score_matrix() -> HashMap<(char, char), i32> {
    let matrix_type = args_parser::get_matrix_type();
    match matrix_type.as_str() {
        "HOXD70.mtx" | "HOXD70" => create_score_matrix_from_matrix_file("HOXD70.mtx"),
        "HOXD55.mtx" | "HOXD55" => create_score_matrix_from_matrix_file("HOXD55.mtx"),
        "none" => {
            let (match_score, mismatch_score) = args_parser::get_match_mismatch();
            create_score_matrix_match_mis(match_score, mismatch_score)
        }
        _ => {
            panic!("wrong matrix type")
        }
    }
}
pub fn create_score_matrix_match_mis(m: i32, x: i32) -> HashMap<(char, char), i32> {
    let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
    for i in ['A', 'C', 'G', 'T', 'N', '-'].iter() {
        for j in ['A', 'C', 'G', 'T', 'N', '-'].iter() {
            if i == j {
                score_matrix.insert((*i, *j), m);
            } else if *i == '-' || *j == '-' {
                score_matrix.insert((*i, *j), x * 2);
            } else {
                score_matrix.insert((*i, *j), x);
            }
        }
    }
    score_matrix.insert(('N', 'N'), x);
    score_matrix.remove(&('-', '-'));
    score_matrix
}
pub fn create_score_matrix_match_mis_f32(m: f32, x: f32) -> HashMap<(char, char), f32> {
    let mut score_matrix: HashMap<(char, char), f32> = HashMap::new();
    for i in ['A', 'C', 'G', 'T', 'N', '-'].iter() {
        for j in ['A', 'C', 'G', 'T', 'N', '-'].iter() {
            if i == j {
                score_matrix.insert((*i, *j), m);
            } else {
                score_matrix.insert((*i, *j), x);
            }
        }
    }
    score_matrix.insert(('N', 'N'), x);
    score_matrix.remove(&('-', '-'));
    score_matrix
}
pub fn create_score_matrix_from_matrix_file(matrix_file: &str) -> HashMap<(char, char), i32> {
    let mut matrix: Vec<Vec<String>> = Vec::new();
    let file_path = project_root::get_project_root().unwrap().join(matrix_file);

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
    matrix_score.remove(&('-', '-'));
    matrix_score
}

#[cfg(test)]
mod tests {
    #[test]
    fn match_miss_matrix_correct() {
        let score_matrix = super::create_score_matrix_match_mis(10, -10);
        assert_eq!(*score_matrix.get(&('A', 'A')).unwrap(), 10);
        assert_eq!(*score_matrix.get(&('A', 'C')).unwrap(), -10);
        assert_eq!(*score_matrix.get(&('N', 'N')).unwrap(), -10);
        assert_eq!(score_matrix.get(&('-', '-')), None);
    }
    #[test]
    fn hoxd_correct() {
        let score_matrix_d70 = super::create_score_matrix_from_matrix_file("HOXD70.mtx");
        let score_matrix_d55 = super::create_score_matrix_from_matrix_file("HOXD55.mtx");

        assert_eq!(*score_matrix_d70.get(&('A', 'A')).unwrap(), 91);
        assert_eq!(*score_matrix_d70.get(&('T', 'G')).unwrap(), -144);

        assert_eq!(*score_matrix_d55.get(&('A', 'A')).unwrap(), 91);
        assert_eq!(*score_matrix_d55.get(&('T', 'G')).unwrap(), -90);

        assert_eq!(score_matrix_d70.get(&('-', '-')), None);
        assert_eq!(score_matrix_d55.get(&('-', '-')), None);
    }
}
