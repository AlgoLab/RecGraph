use crate::utils::idx;

const SIZE: usize = 6;

pub fn create_score_matrix(m: i32, x: i32, i: i32) -> Vec<i32> {
    let mut scores_matrix: Vec<i32> = vec![0; SIZE * SIZE];

    for first_char in "ACGT-N".chars() {
        for second_char in "ACGT-N".chars() {
            let index = idx(first_char as u8, second_char as u8);
            match (first_char, second_char) {
                ('N', _) | (_, 'N') => scores_matrix[index] = x,
                ('-', _) | (_, '-') => scores_matrix[index] = i,
                _ => scores_matrix[index] = if first_char == second_char { m } else { x },
            }
        }
    }
    scores_matrix
}
