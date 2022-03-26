use std::fs::File;
use std::io::{prelude::*, BufReader};
use crate::args_parser;

pub fn get_sequences() -> Vec<String> {
    let file_path = args_parser::get_sequence_path();
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    let mut sequences: Vec<String> = Vec::new();
    for line in reader.lines().flatten() {
        if !line.starts_with('>') && !line.is_empty(){
            let line: String = line
                .chars()
                .map(|c| if c == '-' { 'N' } else { c })
                .collect::<String>()
                .to_uppercase();
            sequences.push(line);
        }
    }
    sequences
}
