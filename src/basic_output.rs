use std::fs::File;
use std::io::{prelude::*, BufWriter};

pub fn write_alignment(path: &[Vec<char>], start_col: usize, s1: &[char], s2: &[char]) {
    let mut row = path.len() - 1;
    let mut col = start_col;
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

    s1_align = s1_align.chars().rev().collect();
    alignment_moves = alignment_moves.chars().rev().collect();
    s2_align = s2_align.chars().rev().collect();

    let path = project_root::get_project_root()
        .unwrap()
        .join("alignment.txt");
    let f = File::create(path).expect("unable to create file");
    let mut f = BufWriter::new(f);

    writeln!(f, "{}", s1_align).expect("unable to write");
    writeln!(f, "{}", alignment_moves).expect("unable to write");
    writeln!(f, "{}", s2_align).expect("unable to write");
}
