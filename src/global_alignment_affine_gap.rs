use std::{cmp, collections::HashMap};
use std::fs::File;
use std::io::{prelude::*, BufWriter};
// guarda pag 244
pub fn exec(s1: &[char], s2: &[char], matrix: &HashMap<(char, char), i32>, o: i32, e: i32) {
    // usato se differenza lunghezze sequenze > lunghezza sequenza più corta
    let mut m = vec![vec![0; s2.len()]; s1.len()]; // best alignment
    let mut x = vec![vec![0; s2.len()]; s1.len()]; // best alignment ending with gap in s1
    let mut y = vec![vec![0; s2.len()]; s1.len()]; // best alignment ending with gap in s2

    let mut path = vec![vec!['x'; s2.len()]; s1.len()];

    for (i, char_1) in s1.iter().enumerate() {
        for (j, char_2) in s2.iter().enumerate() {
            match (i, j) {
                (0, 0) => {
                    path[i][j] = 'O';
                },
                (_, 0) => {
                    x[i][j] = o + e*i as i32;
                    m[i][j] = x[i][j];

                    path[i][j] = 'U';
                },
                (0, _) => {
                    y[i][j] = o + e*j as i32;
                    m[i][j] = y[i][j];

                    path[i][j] = 'L';
                },
                _ => {
                    // set x
                    x[i][j] = cmp::max(
                        x[i][j-1]+e,
                        m[i][j-1]+o+e
                    );
                    
                    // set y
                    y[i][j] = cmp::max(
                        y[i-1][j]+e,
                        m[i-1][j]+o+e
                    );

                    // set m
                    let d = m[i-1][j-1] + matrix.get(&(*char_1, *char_2)).unwrap();
                    let u = y[i][j];
                    let l = x[i][j];

                    match d.cmp(&u) {
                        cmp::Ordering::Less => {
                            match u.cmp(&l) {
                                cmp::Ordering::Less => {
                                    m[i][j] = l;
                                    path[i][j] = 'L';
                                },
                                _ => {
                                    m[i][j] = u;
                                    path[i][j] = 'U';
                                },
                            }
                        },
                        _ => {
                            match d.cmp(&l) {
                                cmp::Ordering::Less => {
                                    m[i][j] = l;
                                    path[i][j] = 'L';
                                },
                                _ => {
                                    m[i][j] = d;
                                    if char_1 == char_2 {
                                        path[i][j] = 'D';
                                    } else {
                                        path[i][j] = 'd';
                                    }
                                }
                            }
                        },
                    }
                }
            }
        }
    }
    println!("gap alignment: {}", m[s1.len()-1][s2.len()-1]); 
    write_alignment(&path, s1.len()-1, s2.len()-1, s1, s2)
}

fn write_alignment(path: &[Vec<char>], mut row: usize, mut col: usize, s1: &[char], s2: &[char]) {
    let mut s1_align = String::new();
    let mut s2_align = String::new();
    let mut alignment_moves = String::new();

    while row > 0 || col >0 {
        match path[row][col] {
            'D' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[col]);
                alignment_moves.push('|');
                row -= 1;
                col -= 1;
            }
            'd' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[col]);
                alignment_moves.push('.');
                col -= 1;
                row -= 1;
            }
            'U' => {
                s1_align.push(s1[row]);
                s2_align.push('-');
                alignment_moves.push(' ');
                row -= 1;
            }
            'L' => {
                s1_align.push('-');
                s2_align.push(s2[col]);
                alignment_moves.push(' ');
                col -= 1;
            }
            _ => panic!("print panic"),
        }
    }

    s1_align = s1_align.chars().rev().collect();
    alignment_moves = alignment_moves.chars().rev().collect();
    s2_align = s2_align.chars().rev().collect();
    let file_name = "gap_alignment.txt";

    let path = project_root::get_project_root().unwrap().join(file_name);
    let f = File::create(path).expect("unable to create file");
    let mut f = BufWriter::new(f);

    writeln!(f, "{}", s1_align).expect("unable to write");
    writeln!(f, "{}", alignment_moves).expect("unable to write");
    writeln!(f, "{}", s2_align).expect("unable to write");
}
