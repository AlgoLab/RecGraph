use crate::basic_output;
use std::{cmp, collections::HashMap};

pub fn ab_loc_alignment(
    s1: &[char],
    s2: &[char],
    score_matrix: &HashMap<(char, char), i32>,
    ampl: i32,
) {
    
    let s1_len = s1.len();
    let s2_len = s2.len();
    
    let mut a = vec![vec![0; ampl as usize]; s1_len];
    let mut path = vec![vec!['x'; ampl as usize]; s1_len];

    a[0][(ampl / 2) as usize] = 0;
    path[0][(ampl / 2) as usize] = 'O';

    for j in (ampl / 2 + 1) as usize..ampl as usize {
        // prima riga
        a[0][j] = 0;
        path[0][j] = 'O';
    }

    for i in 1..s1_len {
        for j in 0..ampl as usize {
            let i_32 = i as i32;
            let j_32 = j as i32;

            if i_32 + j_32 - ampl / 2 < s2_len as i32 && i_32 + j_32 >= ampl / 2 {
                if i_32 + j_32 == ampl / 2 {
                    // elementi con j = 0 matrice mn
                    a[i][j] = 0;
                    path[i][j] = 'O';
                } else if j == 0 {
                    // bordo inf banda
                    let index_of_s2 = i + j - (ampl / 2) as usize;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = a[i - 1][j + 1] + score_matrix.get(&(s1[i], '-')).unwrap();

                    if d < 0 && u < 0 {
                        a[i][j] = 0;
                        path[i][j] = 'O';
                    } else {
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
                    }
                } else if j_32 == ampl - 1 {
                    // bordo sup banda
                    let index_of_s2 = i + j - (ampl / 2) as usize;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let l = a[i][j - 1] + score_matrix.get(&(s2[index_of_s2], '-')).unwrap();

                    if d < 0 && l < 0 {
                        a[i][j] = 0;
                        path[i][j] = 'O';
                    } else {
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
                    }
                } else {
                    // celle interne banda
                    let index_of_s2 = i + j - (ampl / 2) as usize;

                    let d = a[i - 1][j] + score_matrix.get(&(s2[index_of_s2], s1[i])).unwrap();
                    let u = a[i - 1][j + 1] + score_matrix.get(&(s1[i], '-')).unwrap();
                    let l = a[i][j - 1] + score_matrix.get(&(s2[index_of_s2], '-')).unwrap();

                    if d < 0 && u < 0 && l < 0 {
                        a[i][j] = 0;
                        path[i][j] = 'O';
                    } else {
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
    }
    // ampl is enough deve andare da valore max a primo 'o'
    let (row, col) = get_max(&a);
    match ampl_is_enough_iterative(&path, row, col) {
        true => {
            println!("Local Alignement: {}", a[row][col]);

            basic_output::write_alignment(&path, row, col, s1, s2, "local_ab")
        }
        false => ab_loc_alignment(s1, s2, score_matrix, ampl * 2 + 1),
    }
}

fn ampl_is_enough_iterative(path: &[Vec<char>], mut row: usize, mut col: usize) -> bool {
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

fn get_max(a: &[Vec<i32>]) -> (usize, usize) {
    let mut max_row = 0;
    let mut max_col = 0;

    for i in 0..a.len() {
        for j in 0..a[0].len() {
            if a[i][j] >= a[max_row][max_col] {
                max_row = i;
                max_col = j;
            }
        }
    }
    (max_row, max_col)
}
