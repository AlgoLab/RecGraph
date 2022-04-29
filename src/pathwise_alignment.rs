use crate::{bitfield_path as bf, graph::LnzGraph};
use bitvec::prelude::*;
use std::collections::{HashMap, HashSet};

pub fn exec(
    sequence: &[char],
    graph: &LnzGraph,
    path_node: &[Vec<usize>],
    score_matrix: &HashMap<(char, char), i32>,
    path_number: usize,
) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;

    let mut dpm = vec![vec![vec![0; path_number]; sequence.len()]; lnz.len()];
    let mut path = vec![vec![bitvec![u16, Msb0; 0; 32]; sequence.len()]; lnz.len()];

    let mut alphas = vec![path_number + 1; lnz.len()];

    for i in 0..lnz.len() - 1 {
        for j in 0..sequence.len() {
            match (i, j) {
                (0, 0) => {
                    alphas[0] = 0;
                    dpm[i][j] = vec![0; path_number];
                    path[i][j] = bf::set_path_cell(0, 'O');
                }
                (_, 0) => {
                    if !nodes_with_pred[i] {
                        let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                        let pred_paths = path_node[i - 1].iter().collect::<HashSet<_>>();
                        let x = curr_node_paths
                            .intersection(&pred_paths)
                            .collect::<Vec<_>>();

                        if x.contains(&&&alphas[i - 1]) {
                            for k in x.iter() {
                                if ***k == alphas[i - 1] {
                                    dpm[i][j][***k] = dpm[i - 1][j][***k]
                                        + score_matrix.get(&(lnz[i], '-')).unwrap();
                                }
                                dpm[i][j][***k] = dpm[i - 1][j][***k];
                            }
                            alphas[i] = alphas[i - 1];
                        } else {
                            // set new alpha for this node, update delta
                            alphas[i] = ***x.iter().min().unwrap();
                            for k in x.iter() {
                                if ***k == alphas[i] {
                                    dpm[i][j][***k] = dpm[i - 1][j][alphas[i - 1]]
                                        + dpm[i - 1][j][***k]
                                        + score_matrix.get(&(lnz[i], '-')).unwrap();
                                }
                                dpm[i][j][***k] =
                                    dpm[i - 1][j][alphas[i - 1]] + dpm[i - 1][j][***k];
                            }
                        }
                        path[i][j] = bf::set_path_cell(i - 1, 'U');
                    } else {
                        // TODO: consider also multiple alpha
                        for p in pred_hash.get(&i).unwrap().iter() {
                            let curr_node_paths = path_node[i].iter().collect::<HashSet<&usize>>();
                            let pred_paths = path_node[*p].iter().collect::<HashSet<_>>();
                            let x = curr_node_paths
                                .intersection(&pred_paths)
                                .collect::<Vec<_>>();

                            if x.contains(&&&alphas[*p]) {
                                for k in x.iter() {
                                    if ***k == alphas[i - 1] {
                                        dpm[i][j][***k] = dpm[*p][j][***k]
                                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                                    }
                                    dpm[i][j][***k] = dpm[*p][j][***k];
                                }
                                alphas[i] = alphas[*p];
                            } else {
                                // set new alpha for this node, update delta
                                alphas[i] = ***x.iter().min().unwrap();
                                for k in x.iter() {
                                    if ***k == alphas[i] {
                                        dpm[i][j][***k] = dpm[*p][j][alphas[*p]]
                                            + dpm[*p][j][***k]
                                            + score_matrix.get(&(lnz[i], '-')).unwrap();
                                    }
                                    dpm[i][j][***k] =
                                        dpm[*p][j][alphas[*p]] + dpm[*p][j][***k];
                                }
                            }
                        }
                    }
                }
                (0, _) => {
                    dpm[i][j][alphas[0]] =
                        dpm[i][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
                    for k in alphas[0] + 1..path_number {
                        dpm[i][j][k] = dpm[i][j - 1][k];
                    }
                    path[i][j] = bf::set_path_cell(i, 'L');
                }
                _ => {}
            }
        }
    }
}
