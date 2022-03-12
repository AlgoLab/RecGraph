use crate::basic_output;
use std::{cmp, collections::HashMap};

// guarda pag 244
pub fn exec(s1: &[char], s2: &[char], matrix: &HashMap<(char, char), i32>, o: i32, e: i32) {
    // usato se differenza lunghezze sequenze > lunghezza sequenza più corta
    let mut m = vec![vec![0; s2.len()]; s1.len()]; // best alignment
    let mut x = vec![vec![0; s2.len()]; s1.len()]; // best alignment ending with gap in s1
    let mut y = vec![vec![0; s2.len()]; s1.len()]; // best alignment ending with gap in s2

    let mut path = vec![vec!['x'; s2.len()]; s1.len()];
    let mut path_x = vec![vec!['M'; s2.len()]; s1.len()];
    let mut path_y = vec![vec!['M'; s2.len()]; s1.len()];

    for (i, char_1) in s1.iter().enumerate() {
        for (j, char_2) in s2.iter().enumerate() {
            match (i, j) {
                (0, 0) => {
                    path[i][j] = 'O';
                }
                (_, 0) => {
                    x[i][j] = o + e * i as i32;
                    m[i][j] = x[i][j];

                    path[i][j] = 'U';
                }
                (0, _) => {
                    y[i][j] = o + e * j as i32;
                    m[i][j] = y[i][j];

                    path[i][j] = 'L';
                }
                _ => {
                    // set x
                    x[i][j] = cmp::max(x[i][j - 1] + e, m[i][j - 1] + o + e);
                    if x[i][j] != m[i][j - 1] + o + e {
                        path_x[i][j - 1] = 'X'
                    }

                    // set y
                    y[i][j] = cmp::max(y[i - 1][j] + e, m[i - 1][j] + o + e);
                    if y[i][j] != m[i - 1][j] + o + e {
                        path_y[i - 1][j] = 'Y'
                    }
                    // set m
                    let d = m[i - 1][j - 1] + matrix.get(&(*char_1, *char_2)).unwrap();
                    let u = y[i][j];
                    let l = x[i][j];

                    match d.cmp(&u) {
                        cmp::Ordering::Less => match u.cmp(&l) {
                            cmp::Ordering::Less => {
                                m[i][j] = l;
                                path[i][j] = 'L';
                            }
                            _ => {
                                m[i][j] = u;
                                path[i][j] = 'U';
                            }
                        },
                        _ => match d.cmp(&l) {
                            cmp::Ordering::Less => {
                                m[i][j] = l;
                                path[i][j] = 'L';
                            }
                            _ => {
                                m[i][j] = d;
                                if char_1 == char_2 {
                                    path[i][j] = 'D';
                                } else {
                                    path[i][j] = 'd';
                                }
                            }
                        },
                    }
                }
            }
        }
    }
    println!("gap alignment: {}", m[s1.len() - 1][s2.len() - 1]);
    basic_output::write_alignment_gap(
        &path,
        &path_x,
        &path_y,
        s1.len() - 1,
        s2.len() - 1,
        s1,
        s2,
        "gap",
    )
}
