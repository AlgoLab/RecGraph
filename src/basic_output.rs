use std::fs::File;
use std::io::{prelude::*, BufWriter};

pub fn write_alignment_ab(
    path: &[Vec<char>],
    mut row: usize,
    mut col: usize,
    s1: &[char],
    s2: &[char],
    align_type: &str,
) {
    let col_number = path[0].len();

    let mut s1_align = String::new();
    let mut s2_align = String::new();
    let mut alignment_moves = String::new();

    while path[row][col] != 'O' {
        let index_of_s2 = row + col - col_number / 2;
        match path[row][col] {
            'D' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[index_of_s2]);
                alignment_moves.push('|');
                row -= 1;
            }
            'd' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[index_of_s2]);
                alignment_moves.push('.');
                row -= 1;
            }
            'U' => {
                s1_align.push(s1[row]);
                s2_align.push('-');
                alignment_moves.push(' ');
                row -= 1;
                col += 1;
            }
            'L' => {
                s1_align.push('-');
                s2_align.push(s2[index_of_s2]);
                alignment_moves.push(' ');
                col -= 1;
            }
            _ => panic!("ampl_is_enough panic"),
        }
    }
    reverse_and_write(s1_align, s2_align, alignment_moves, align_type);
}

pub fn write_alignment_no_ab(
    path: &[Vec<char>],
    mut row: usize,
    mut col: usize,
    s1: &[char],
    s2: &[char],
    align_type: &str,
) {
    let mut s1_align = String::new();
    let mut s2_align = String::new();
    let mut alignment_moves = String::new();

    while path[row][col] != 'O' {
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
            _ => panic!("ampl_is_enough panic"),
        }
    }
    reverse_and_write(s1_align, s2_align, alignment_moves, align_type);
}

pub fn write_alignment_ab_gap(
    path: &[Vec<char>],
    path_x: &[Vec<char>],
    path_y: &[Vec<char>],
    s1: &[char],
    s2: &[char],
    align_type: &str,
) {
    let col_number = path[0].len();

    let mut s1_align = String::new();
    let mut s2_align = String::new();
    let mut alignment_moves = String::new();

    let mut row = s1.len() - 1;
    let mut col = s2.len() - 1 + (col_number / 2) - (s1.len() - 1);

    while path[row][col] != 'O' {
        let index_of_s2 = row + col - col_number / 2;
        match path[row][col] {
            'D' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[index_of_s2]);
                alignment_moves.push('|');
                row -= 1;
            }
            'd' => {
                s1_align.push(s1[row]);
                s2_align.push(s2[index_of_s2]);
                alignment_moves.push('.');
                row -= 1;
            }
            'U' => {
                if path_y[row][col] == 'M' {
                    s1_align.push(s1[row]);
                    s2_align.push('-');
                    alignment_moves.push(' ');
                    row -= 1;
                    col += 1;
                } else {
                    while path_y[row][col] == 'Y' {
                        s1_align.push(s1[row]);
                        s2_align.push('-');
                        alignment_moves.push(' ');
                        row -= 1;
                        col += 1;
                    }
                }
            }
            'L' => {
                if path_x[row][col] == 'M' {
                    s1_align.push('-');
                    s2_align.push(s2[row + col - col_number / 2]);
                    alignment_moves.push(' ');
                    col -= 1;
                } else {
                    while path_x[row][col] == 'X' {
                        s1_align.push('-');
                        s2_align.push(s2[row + col - col_number / 2]);
                        alignment_moves.push(' ');
                        col -= 1;
                    }
                }
            }
            _ => panic!("ampl_is_enough panic"),
        }
    }
    reverse_and_write(s1_align, s2_align, alignment_moves, align_type);
}

pub fn write_alignment_gap(
    path: &[Vec<char>],
    path_x: &[Vec<char>],
    path_y: &[Vec<char>],
    s1: &[char],
    s2: &[char],
    align_type: &str,
) {
    let mut row = s1.len() - 1;
    let mut col = s2.len() - 1;
    let mut s1_align = String::new();
    let mut s2_align = String::new();
    let mut alignment_moves = String::new();

    while path[row][col] != 'O' {
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
                if path_y[row][col] == 'M' {
                    s1_align.push(s1[row]);
                    s2_align.push('-');
                    alignment_moves.push(' ');
                    row -= 1;
                } else {
                    while path_y[row][col] == 'Y' {
                        s1_align.push(s1[row]);
                        s2_align.push('-');
                        alignment_moves.push(' ');
                        row -= 1;
                    }
                }
            }
            'L' => {
                if path_x[row][col] == 'M' {
                    s1_align.push('-');
                    s2_align.push(s2[col]);
                    alignment_moves.push(' ');
                    col -= 1;
                } else {
                    while path_x[row][col] == 'X' {
                        s1_align.push('-');
                        s2_align.push(s2[col]);
                        alignment_moves.push(' ');
                        col -= 1;
                    }
                }
            }
            _ => panic!("ampl_is_enough panic"),
        }
    }
    reverse_and_write(s1_align, s2_align, alignment_moves, align_type);
}

pub fn write_align_poa(path: &[Vec<(char, i32)>], sequence: &[char], graph: &[(char, Vec<usize>)]) {
    let mut row = path.len() - 1;
    let mut col = path[path.len() - 1][path[0].len() - 1].1 as usize;
    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();

    while path[row][col] != ('O', 0) {
        match path[row][col] {
            ('D', _) => {
                sequence_align.push(sequence[row]);
                graph_align.push(graph[col].0);
                alignment_moves.push('|');
                row -= 1;
                col = path[row][col].1 as usize;
            }
            ('d', _) => {
                sequence_align.push(sequence[row]);
                graph_align.push(graph[col].0);
                alignment_moves.push('.');
                row -= 1;
                col = path[row][col].1 as usize;
            }
            ('L', _) => {
                sequence_align.push('-');
                graph_align.push(graph[col].0);
                alignment_moves.push(' ');
                col = path[row][col].1 as usize;
            }
            ('U', _) => {
                sequence_align.push(sequence[row]);
                graph_align.push('-');
                alignment_moves.push(' ');
                row -= 1;
            }
            _ => {
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(sequence_align, graph_align, alignment_moves, "poa");
}

pub fn write_align_ab_poa(
    path: &[Vec<(char, i32)>],
    sequence: &[char],
    graph: &[(char, Vec<usize>)],
) {
    let mut row = path[graph.len() - 1]
        [(sequence.len() - 1) + (path[0].len() / 2) - (graph.len() - 1)]
        .1 as usize;
    let mut col = (sequence.len() - 1) + (path[0].len() / 2) - (graph.len() - 1);

    let mut sequence_align = String::new();
    let mut graph_align = String::new();
    let mut alignment_moves = String::new();

    while path[row][col].0 != 'O' {
        let index_of_seq = row + col - path[0].len() / 2;

        match path[row][col] {
            ('D', _) => {
                sequence_align.push(sequence[index_of_seq]);
                graph_align.push(graph[row].0);
                alignment_moves.push('|');
                row = path[row][col].1 as usize;
            }
            ('d', _) => {
                sequence_align.push(sequence[index_of_seq]);
                graph_align.push(graph[row].0);
                alignment_moves.push('.');
                row = path[row][col].1 as usize;
            }
            ('L', _) => {
                sequence_align.push(sequence[index_of_seq]);
                graph_align.push('-');
                alignment_moves.push(' ');

                col -= 1;
            }
            ('U', _) => {
                sequence_align.push('-');
                graph_align.push(graph[row].0);
                alignment_moves.push(' ');

                let delta = row - path[row][col].1 as usize;
                row = path[row][col].1 as usize;
                col += delta;
            }
            _ => {
                println!("{} {} {:?}", row, col, path[row][col]);
                panic!("impossible value in poa path")
            }
        }
    }
    reverse_and_write(sequence_align, graph_align, alignment_moves, "ab_poa");
}

fn reverse_and_write(mut s1_al: String, mut s2_al: String, mut al_moves: String, align_type: &str) {
    s1_al = s1_al.chars().rev().collect();
    al_moves = al_moves.chars().rev().collect();
    s2_al = s2_al.chars().rev().collect();
    let file_name = String::from(align_type) + "_alignment.txt";

    let path = project_root::get_project_root().unwrap().join(file_name);
    let f = File::create(path).expect("unable to create file");
    let mut f = BufWriter::new(f);

    writeln!(f, "{}", s1_al).expect("unable to write");
    writeln!(f, "{}", al_moves).expect("unable to write");
    writeln!(f, "{}", s2_al).expect("unable to write");
}
