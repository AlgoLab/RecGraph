use crate::basic_output;
use std::{cmp, collections::HashMap};
pub fn exec(s1: &[char], s2: &[char], matrix: &HashMap<(char, char), i32>) {
    let mut a = vec![vec![0; s2.len()]; s1.len()];
    let mut path = vec![vec!['x'; s2.len()]; s1.len()];
    let mut max_row = 0;
    let mut max_col = 0;

    for (row, char_s1) in s1.iter().enumerate() {
        for (col, char_s2) in s2.iter().enumerate() {
            match (row, col) {
                (_, 0) | (0, _) => {
                    a[row][col] = 0;
                    path[row][col] = 'O';
                }
                _ => {
                    let d = a[row - 1][col - 1] + matrix.get(&(*char_s1, *char_s2)).unwrap();
                    let l = a[row - 1][col] + matrix.get(&(*char_s1, '-')).unwrap();
                    let u = a[row][col - 1] + matrix.get(&('-', *char_s2)).unwrap();

                    if d < 0 && l < 0 && u < 0 {
                        a[row][col] = 0;
                        path[row][col] = 'O';
                    }
                    match d.cmp(&l) {
                        cmp::Ordering::Less => match l.cmp(&u) {
                            cmp::Ordering::Less => {
                                a[row][col] = u;
                                path[row][col] = 'U'
                            }
                            _ => {
                                a[row][col] = l;
                                path[row][col] = 'L'
                            }
                        },
                        _ => match d.cmp(&u) {
                            cmp::Ordering::Less => {
                                a[row][col] = u;
                                path[row][col] = 'U'
                            }
                            _ => {
                                a[row][col] = d;
                                if s1[row] == s2[col] {
                                    path[row][col] = 'D'
                                } else {
                                    path[row][col] = 'd'
                                }
                            }
                        },
                    }
                }
            }
            if a[row][col] >= a[max_row][max_col] {
                max_row = row;
                max_col = col;
            }
        }
    }
    println!("Local Alignement: {}", a[max_row][max_col]);
    basic_output::write_alignment_no_ab(&path, max_row, max_col, s1, s2, "local")
}
