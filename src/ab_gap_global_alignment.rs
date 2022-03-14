use crate::basic_output;
use std::{cmp, collections::HashMap};

pub fn exec(
    s1: &[char],
    s2: &[char],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: usize,
    o: i32,
    e: i32,
) -> i32 {
    let s1_len = s1.len();
    let s2_len = s2.len();

    let mut m = vec![vec![0; ampl]; s1_len];
    let mut x = vec![vec![0; ampl]; s1_len];
    let mut y = vec![vec![0; ampl]; s1_len];

    let mut path = vec![vec!['x'; ampl]; s1_len];
    let mut path_x = vec![vec!['M'; ampl]; s1_len];
    let mut path_y = vec![vec!['M'; ampl]; s1_len];

    m[0][ampl / 2] = 0;
    path[0][ampl / 2] = 'O';

    for j in ampl / 2 + 1..ampl {
        // prima riga i = 0
        if j - (ampl / 2) < s2_len {
            y[0][j] = o + e * (j - ampl / 2) as i32;
            m[0][j] = y[0][j];
            path[0][j] = 'L';
        }
    }

    for i in 1..s1_len {
        for j in 0..ampl {
            if i + j < s2_len + ampl / 2 && i + j >= ampl / 2 {
                // else out of band
                if i + j == ampl / 2 {
                    // j = 0 in mn
                    x[i][j] = o + e * i as i32;
                    m[i][j] = x[i][j];

                    path[i][j] = 'U';
                } else if j == 0 {
                    // bordo inf banda
                    // set y
                    y[i][j] = cmp::max(y[i - 1][j + 1] + e, m[i - 1][j + 1] + o + e);
                    if y[i][j] != m[i - 1][j + 1] + o + e {
                        path_y[i - 1][j + 1] = 'Y'
                    }

                    //set m
                    let index_of_s2 = i + j - ampl / 2;

                    let d = m[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = y[i][j];
                    match d.cmp(&u) {
                        cmp::Ordering::Less => {
                            m[i][j] = u;
                            path[i][j] = 'U'
                        }
                        _ => {
                            m[i][j] = d;
                            if s1[i] == s2[index_of_s2] {
                                path[i][j] = 'D'
                            } else {
                                path[i][j] = 'd'
                            }
                        }
                    }
                } else if j == ampl - 1 {
                    // bordo sup banda
                    // set x
                    x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);
                    if x[i][j] != m[i][j - 1] + o + e {
                        path_x[i][j - 1] = 'X'
                    }

                    // set m
                    let index_of_s2 = i + j - ampl / 2;

                    let d = m[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let l = x[i][j];
                    match d.cmp(&l) {
                        cmp::Ordering::Less => {
                            m[i][j] = l;
                            path[i][j] = 'L'
                        }
                        _ => {
                            m[i][j] = d;
                            if s1[i] == s2[index_of_s2] {
                                path[i][j] = 'D'
                            } else {
                                path[i][j] = 'd'
                            }
                        }
                    }
                } else {
                    // celle interne banda
                    // set x
                    x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);
                    if x[i][j] != m[i][j - 1] + o + e {
                        path_x[i][j - 1] = 'X'
                    }

                    // set y
                    y[i][j] = cmp::max(y[i - 1][j + 1] + e, m[i - 1][j + 1] + o + e);
                    if y[i][j] != m[i - 1][j + 1] + o + e {
                        path_y[i - 1][j + 1] = 'Y'
                    }

                    // set m
                    let index_of_s2 = i + j - ampl / 2;

                    let d = m[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = y[i][j];
                    let l = x[i][j];

                    match d.cmp(&l) {
                        cmp::Ordering::Less => match l.cmp(&u) {
                            cmp::Ordering::Less => {
                                m[i][j] = u;
                                path[i][j] = 'U'
                            }
                            _ => {
                                m[i][j] = l;
                                path[i][j] = 'L'
                            }
                        },
                        _ => match d.cmp(&u) {
                            cmp::Ordering::Less => {
                                m[i][j] = u;
                                path[i][j] = 'U';
                            }
                            _ => {
                                m[i][j] = d;
                                if s1[i] == s2[index_of_s2] {
                                    path[i][j] = 'D'
                                } else {
                                    path[i][j] = 'd'
                                }
                            }
                        },
                    }
                }
            }
        }
    }

    match ampl_is_enough_iterative(
        &path,
        &path_x,
        &path_y,
        s2_len - 1 + (ampl / 2) - (s1_len - 1),
    ) {
        true => {
            println!(
                "Ab Gap Alignement: {}",
                m[s1_len - 1][s2_len - 1 + (ampl / 2) - (s1_len - 1)]
            );
            basic_output::write_alignment_ab_gap(&path, &path_x, &path_y, s1, s2, "ab_gap");
            m[s1_len - 1][s2_len - 1 + (ampl / 2) - (s1_len - 1)]
        }
        false => exec(s1, s2, score_matrix, ampl * 2 + 1, o, e),
    }
}

fn ampl_is_enough_iterative(
    path: &[Vec<char>],
    path_x: &[Vec<char>],
    path_y: &[Vec<char>],
    start_col: usize,
) -> bool {
    let mut row = path.len() - 1;
    let mut col = start_col;
    let col_number = path[0].len();
    while path[row][col] != 'O' {
        if col == 0 || col == col_number - 1 {
            if path[row][col] == 'D' {
                row -= 1; // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            match path[row][col] {
                'D' | 'd' => {
                    row -= 1;
                }
                'U' => {
                    if path_y[row][col] == 'M' {
                        row -= 1;
                        col += 1;
                    } else {
                        while path_y[row][col] == 'Y' {
                            if col == 0 || col == col_number - 1 {
                                return false;
                            }
                            row -= 1;
                            col += 1;
                        }
                    }
                }
                'L' => {
                    if path_x[row][col] == 'M' {
                        col -= 1;
                    } else {
                        while path_x[row][col] == 'X' {
                            if col == 0 || col == col_number - 1 {
                                return false;
                            }
                            col -= 1;
                        }
                    }
                }
                _ => panic!("ampl_is_enough panic "),
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    #[test]
    fn correct_score_from_ab_gap_align() {
        let s1: Vec<char> = "$AATTAACTTTCGC".chars().collect();
        let s2: Vec<char> = "$AACCC".chars().collect();

        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                score_matrix.insert((*c1, *c2), 1);
                } else {
                score_matrix.insert((*c1, *c2), -1);
                }
            }
        }
        let align_score = super::exec(&s1, &s2, &score_matrix, (s1.len()-s2.len())*2+1, -10, -1);
        assert_eq!(align_score, -15);
    }
    #[test]
    fn same_string_align_best_possible_score() {
        let s1: Vec<char> = "$AATTAACTTTCGC".chars().collect();
        
        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                score_matrix.insert((*c1, *c2), 1);
                } else {
                score_matrix.insert((*c1, *c2), -1);
                }
            }
        }
        let align_score = super::exec(&s1, &s1, &score_matrix, 3, -10, -1);
        assert_eq!(align_score, (s1.len()-1) as i32);
    }
    #[test]
    fn empty_string_align_zero() {
        let s1: Vec<char> = "$".chars().collect();
        
        let mut score_matrix: HashMap<(char, char), i32> = HashMap::new();
        for c1 in ['A', 'C', 'G', 'T'].iter() {
            for c2 in ['A', 'C', 'G', 'T'].iter() {
                if c1 == c2 {
                score_matrix.insert((*c1, *c2), 1);
                } else {
                score_matrix.insert((*c1, *c2), -1);
                }
            }
        }
        let align_score = super::exec(&s1, &s1, &score_matrix, 3, -10, -1);
        assert_eq!(align_score, 0);
    }
}