use crate::basic_output;
use std::{cmp, collections::HashMap};

pub fn exec(s1: &[char], s2: &[char], score_matrix: &HashMap<(char, char), i32>, ampl: usize) {
    let s1_len = s1.len();
    let s2_len = s2.len();
    let mut a = vec![vec![0; ampl as usize]; s1_len];
    let mut path = vec![vec!['x'; ampl as usize]; s1_len];

    a[0][(ampl / 2) as usize] = 0;
    path[0][(ampl / 2) as usize] = 'O';

    for j in ampl / 2 + 1..ampl {
        // prima riga
        if j - (ampl / 2) < s2_len {
            a[0][j] = a[0][j - 1] + score_matrix.get(&('-', s2[j - ampl / 2])).unwrap();
            path[0][j] = 'L';
        }
    }

    for i in 1..s1_len {
        for j in 0..ampl {
            if i + j < s2_len + ampl / 2 && i + j >= ampl / 2 {
                if i + j == ampl / 2 {
                    // elementi con j = 0 matrice mn
                    a[i][j] = a[i - 1][j + 1] + score_matrix.get(&(s1[i], '-')).unwrap();
                    path[i][j] = 'U';
                } else if j == 0 {
                    // bordo inf banda
                    let index_of_s2 = i + j - ampl / 2;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = a[i - 1][j + 1] + score_matrix.get(&(s1[i], '-')).unwrap();
                    match d.cmp(&u) {
                        cmp::Ordering::Less => {
                            a[i][j] = u;
                            path[i][j] = 'U'
                        }
                        _ => {
                            a[i][j] = d;
                            if s1[i] == s2[index_of_s2] {
                                path[i][j] = 'D'
                            } else {
                                path[i][j] = 'd'
                            }
                        }
                    }
                } else if j == ampl - 1 {
                    // bordo sup banda
                    let index_of_s2 = i + j - ampl / 2;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let l = a[i][j - 1] + score_matrix.get(&(s2[index_of_s2], '-')).unwrap();
                    match d.cmp(&l) {
                        cmp::Ordering::Less => {
                            a[i][j] = l;
                            path[i][j] = 'L'
                        }
                        _ => {
                            a[i][j] = d;
                            if s1[i] == s2[index_of_s2] {
                                path[i][j] = 'D'
                            } else {
                                path[i][j] = 'd'
                            }
                        }
                    }
                } else {
                    // celle interne banda
                    let index_of_s2 = i + j - ampl / 2;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = a[i - 1][j + 1] + score_matrix.get(&(s1[i], '-')).unwrap();
                    let l = a[i][j - 1] + score_matrix.get(&(s2[index_of_s2], '-')).unwrap();

                    match d.cmp(&l) {
                        cmp::Ordering::Less => match l.cmp(&u) {
                            cmp::Ordering::Less => {
                                a[i][j] = u;
                                path[i][j] = 'U'
                            }
                            _ => {
                                a[i][j] = l;
                                path[i][j] = 'L'
                            }
                        },
                        _ => match d.cmp(&u) {
                            cmp::Ordering::Less => {
                                a[i][j] = u;
                                path[i][j] = 'U';
                            }
                            _ => {
                                a[i][j] = d;
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

    match ampl_is_enough_iterative(&path, s2_len - 1 + (ampl / 2) - (s1_len - 1)) {
        true => {
            println!(
                "Global Alignement: {}",
                a[s1_len - 1][s2_len - 1 + (ampl / 2) - (s1_len - 1)]
            );
            basic_output::write_alignment(
                &path,
                s1_len - 1,
                s2_len - 1 + (ampl / 2) - (s1_len - 1),
                s1,
                s2,
                "global",
            )
        }
        false => exec(s1, s2, score_matrix, ampl * 2 + 1),
    }
}

fn ampl_is_enough_iterative(path: &[Vec<char>], start_col: usize) -> bool {
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
                    row -= 1;
                    col += 1;
                }
                'L' => {
                    col -= 1;
                }
                _ => panic!("ampl_is_enough panic"),
            }
        }
    }
    true
}
