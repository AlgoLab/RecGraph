mod basic_edit_distance;
use basic_edit_distance as bed;

mod ab_edit_distance;
use ab_edit_distance as abed;

mod ab_mk_edit_distance;
use ab_mk_edit_distance as abmked;

use std::{cmp, time::Instant};
use std::fs::File;
use std::io::{prelude::*, BufReader};

fn main() {
    // fare esecuzioni in loop e media dei tempi?
   
    let mut sequences: Vec<String> = Vec::new();
    let file = File::open("./sequences.txt").unwrap();
    let reader = BufReader::new(file);
    
    for line in reader.lines () {
        if let Ok(seq) = line {
            sequences.push(seq);
        }
    }

    let mut s1: Vec<char> = sequences[0].chars().collect();
    let mut s2: Vec<char> = sequences[1].chars().collect();
    s1.insert(0, '$');
    s2.insert(0, '$');
    

    //prima versione
    println!("Basic edit distance\n");
    let start = Instant::now();
    bed::exec_ed(&s1, &s2);
    let basic_duration = start.elapsed();

    //con banda su matrice m*n
    let ampl;
    match s1.len() < s2.len() {
        true => ampl = s2.len()-s1.len(),
        _ => ampl = s1.len()-s2.len()
    };
    println!("");
    println!("Band edit distance (m*n)\n");
    let start = Instant::now();
    let mut a =vec![vec![0;s2.len()];s1.len()];
    let mut path =vec![vec!['x';s2.len()];s1.len()];
    abed::ab_ed(&s1, &s2, cmp::max((ampl * 2 + 1) as i32, 3), &mut a, &mut path);
    println!("");
    let ab_duration = start.elapsed();

    //con banda su matrice m*k
    println!("");
    println!("Band edit distance (m*k)\n");
    let start = Instant::now();
    abmked::ab_ed_km(&s1, &s2, cmp::max((ampl * 2 + 1) as i32, 3));
    let ab_mk_duration = start.elapsed();

    println!("");
    println!("BASIC: {:?}", basic_duration);
    println!("BAND: {:?}", ab_duration);
    println!("MK BAND: {:?}", ab_mk_duration);

}