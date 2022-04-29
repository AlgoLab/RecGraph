use std::collections::{HashMap, HashSet};
use bitvec::prelude::*;
use crate::{graph::LnzGraph, bitfield_path as bf};

pub fn exec(
    sequence: &[char],
    graph: &LnzGraph,
    path_node: &[Vec<usize>],
    score_matrix: &HashMap<(char, char), i32>,
    path_number: usize
) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;

    let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let alpha = 0;
    let mut path = vec![vec![bitvec![u16, Msb0; 0; 32]; sequence.len()]; lnz.len()];

    for i in 0..lnz.len()-1{
        for j in 0..sequence.len(){
            match (i, j) {
                (0,0) => {
                    dpm[i][j] = vec![0; path_number];
                    path[i][j] = bf::set_path_cell(0, 'O');
                },
                (_, 0) => {
                    if nodes_with_pred[i] {
                        let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                        let pred_paths = path_node[i-1].iter().collect::<HashSet<_>>();
                        let x = curr_node_paths.intersection(&pred_paths).collect::<HashSet<_>>();
                        
                        if x.contains(&&alpha) {

                        }
                    } else {

                    }
                },
                (0, _) => {
                    dpm[i][j][alpha] = dpm[i][j-1][alpha] + score_matrix.get(&(sequence[j], '-')).unwrap();
                    for k in alpha+1..path_number{
                        dpm[i][j][k] = dpm[i][j-1][k];
                    }
                    path[i][j] = bf::set_path_cell(i, 'L');
                },
                _ => {}
            }
        }
    }
}