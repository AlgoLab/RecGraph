use crate::args_parser;
use std::fs::File;
use std::io::{prelude::*, BufReader};
pub fn get_sequences() -> (Vec<Vec<char>>, Vec<String>) {
    let file_path = args_parser::get_sequence_path();
    let file = File::open(file_path).unwrap();
    let reader = BufReader::new(file);

    let mut sequences: Vec<Vec<char>> = Vec::new();
    let mut sequences_name: Vec<String> = Vec::new();

    let mut sequence: Vec<char> = Vec::new();
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
            sequence.append(&mut line);
        } else if line.starts_with('>') {
            sequences_name.push(line);
            if !sequence.is_empty() {
                sequence.insert(0, '$');
                sequences.push(sequence);
            }
            sequence = Vec::new();
        }
    }
    if !sequence.is_empty() {
        sequences.push(sequence);
    }

    if sequences.len() != sequences_name.len() {
        panic!("wrong fasta file format");
    }
    (sequences, sequences_name) //update with also sequences_name
}
//TODO: verify score from rev&cmpl string equal score from rev&cmpl graph
pub fn rev_and_compl(seq: &[char]) -> Vec<char> {
    let mut rev_seq = seq[1..]
        .iter()
        .map(|c| match *c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'N' => 'N',
            _ => {
                panic!("wrong char: {}, unable to rev&compl", c)
            }
        })
        .collect::<Vec<char>>();
    rev_seq.reverse();
    rev_seq.insert(0, '$');
    rev_seq
}
#[cfg(test)]
mod tests {
    #[test]
    fn rev_and_compl_of_seq_correct() {
        let s1 = ['$', 'A', 'A', 'T'];
        let s1_rc = super::rev_and_compl(&s1);
        for i in 0..s1_rc.len() {
            assert_eq!(['$', 'A', 'T', 'T'][i], s1_rc[i]);
        }
    }
    #[test]
    fn rev_and_compl_with_every_symbol() {
        let s1 = ['$', 'A', 'T', 'C', 'G', 'N'];
        let s1_rc = super::rev_and_compl(&s1);
        for i in 0..s1_rc.len() {
            assert_eq!(['$', 'N', 'C', 'G', 'A', 'T'][i], s1_rc[i]);
        }
    }
}
