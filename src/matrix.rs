use std::collections::HashMap;

pub fn create_score_matrix(m: i32, x: i32) -> HashMap<(char, char), i32> {
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
