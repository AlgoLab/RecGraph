use crate::args_parser;
use std::fs::File;
use std::io::{prelude::*, BufReader};

pub fn get_sequences() -> Vec<Vec<char>> {
    let file_path = args_parser::get_sequence_path();
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    let mut sequences: Vec<Vec<char>> = Vec::new();
    for line in reader.lines().flatten() {
        if !line.starts_with('>') && !line.is_empty() {
            let mut line: Vec<char> = line
                .chars()
                .map(|c| {
                    if c == '-' {
                        'N'
                    } else {
                        c.to_ascii_uppercase()
                    }
                })
                .collect::<Vec<char>>();
            line.insert(0, '$');
            sequences.push(line);
        }
    }
    sequences
}
