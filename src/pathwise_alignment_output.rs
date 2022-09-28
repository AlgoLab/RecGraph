use std::collections::HashMap;

use bit_vec::BitVec;

use crate::{
    gaf_output::GAFStruct,
    pathwise_graph::{PathGraph, PredHash},
};

pub fn build_alignment(
    dpm: &Vec<Vec<Vec<i32>>>,
    graph: &PathGraph,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    best_path: usize,
) -> GAFStruct {
    let pred_hash = &graph.pred_hash;
    let alphas = &graph.alphas;
    let nwp = &graph.nwp;
    let handles_nodes_id = &graph.nodes_id_pos;
    let lnz = &graph.lnz;

    let mut path_length: usize = 0;
    let mut cigar = Vec::new();
    let mut handle_id_alignment = Vec::new();
    let mut i = 0;
    let ending_nodes = pred_hash.get_preds_and_paths(dpm.len() - 1);
    for (node, paths) in ending_nodes.iter() {
        if paths[best_path] {
            i = *node;
        }
    }
    let ending_node = i;

    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
        let mut predecessor = None;
        let curr_score = if alphas[i] == best_path {
            dpm[i][j][best_path]
        } else {
            dpm[i][j][best_path] + dpm[i][j][alphas[i]]
        };
        let (d, u, l) = if !nwp[i] {
            (
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j - 1][best_path]
                } else {
                    dpm[i - 1][j - 1][best_path] + dpm[i - 1][j - 1][alphas[i - 1]]
                } + scores.get(&(lnz[i], seq[j])).unwrap(),
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j][best_path]
                } else {
                    dpm[i - 1][j][best_path] + dpm[i - 1][j][alphas[i - 1]]
                } + scores.get(&(lnz[i], '-')).unwrap(),
                if alphas[i] == best_path {
                    dpm[i][j - 1][best_path]
                } else {
                    dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]]
                } + scores.get(&('-', seq[j])).unwrap(),
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    if alphas[*pred] == best_path {
                        d = dpm[*pred][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap();
                        u = dpm[*pred][j][best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    } else {
                        d = dpm[*pred][j - 1][best_path]
                            + dpm[*pred][j - 1][alphas[*pred]]
                            + scores.get(&(lnz[i], seq[j])).unwrap();
                        u = dpm[*pred][j][best_path]
                            + dpm[*pred][j][alphas[*pred]]
                            + scores.get(&(lnz[i], '-')).unwrap();
                    }
                    if alphas[i] == best_path {
                        l = dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap();
                    } else {
                        l = dpm[i][j - 1][best_path]
                            + dpm[i][j - 1][alphas[i]]
                            + scores.get(&('-', seq[j])).unwrap();
                    }
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
            path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        } else {
            cigar.push('L');
            j -= 1;
        }
    }
    while j > 0 {
        cigar.push('L');
        j -= 1;
    }
    while i > 0 {
        let predecessor = if !nwp[i] {
            i - 1
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let mut p = 0;
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    p = *pred;
                }
            }
            p
        };
        cigar.push('U');
        handle_id_alignment.push(handles_nodes_id[i]);
        i = predecessor;
        path_length += 1;
    }
    cigar.reverse();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    let path: Vec<usize> = handle_id_alignment.iter().map(|id| *id as usize).collect();
    //path length already set
    let path_start = 0;
    let path_end = get_node_offset(handles_nodes_id, ending_node) as usize; // last letter used in last node of alignment

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = format!("{}, best path: {}", build_cigar(&cigar), best_path);
    let gaf_output = GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_length,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}

pub fn build_alignment_semiglobal(
    dpm: &Vec<Vec<Vec<i32>>>,
    lnz: &Vec<char>,
    seq: &[char],
    scores: &HashMap<(char, char), i32>,
    alphas: &Vec<usize>,
    best_path: usize,
    pred_hash: &PredHash,
    nwp: &BitVec,
    handles_nodes_id: &Vec<u64>,
    ending_node: usize,
) -> GAFStruct {
    let mut cigar = Vec::new();
    let mut path_length: usize = 0;
    let mut i = ending_node;
    let mut j = dpm[i].len() - 1;
    let mut handle_id_alignment = Vec::new();

    while i > 0 && j > 0 {
        let curr_score = if alphas[i] == best_path {
            dpm[i][j][best_path]
        } else {
            dpm[i][j][best_path] + dpm[i][j][alphas[i]]
        };
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j - 1][best_path]
                } else {
                    dpm[i - 1][j - 1][best_path] + dpm[i - 1][j - 1][alphas[i - 1]]
                } + scores.get(&(lnz[i], seq[j])).unwrap(),
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j][best_path]
                } else {
                    dpm[i - 1][j][best_path] + dpm[i - 1][j][alphas[i - 1]]
                } + scores.get(&(lnz[i], '-')).unwrap(),
                if alphas[i] == best_path {
                    dpm[i][j - 1][best_path]
                } else {
                    dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]]
                } + scores.get(&('-', seq[j])).unwrap(),
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    if alphas[*pred] == best_path {
                        d = dpm[*pred][j - 1][best_path] + scores.get(&(lnz[i], seq[j])).unwrap();
                        u = dpm[*pred][j][best_path] + scores.get(&(lnz[i], '-')).unwrap();
                    } else {
                        d = dpm[*pred][j - 1][best_path]
                            + dpm[*pred][j - 1][alphas[*pred]]
                            + scores.get(&(lnz[i], seq[j])).unwrap();
                        u = dpm[*pred][j][best_path]
                            + dpm[*pred][j][alphas[*pred]]
                            + scores.get(&(lnz[i], '-')).unwrap();
                    }
                    if alphas[i] == best_path {
                        l = dpm[i][j - 1][best_path] + scores.get(&('-', seq[j])).unwrap();
                    } else {
                        l = dpm[i][j - 1][best_path]
                            + dpm[i][j - 1][alphas[i]]
                            + scores.get(&('-', seq[j])).unwrap();
                    }
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
            path_length += 1;
        } else if max == u {
            cigar.push('U');
            handle_id_alignment.push(handles_nodes_id[i]);
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            path_length += 1;
        } else {
            cigar.push('L');
            j -= 1;
        }
    }
    while j > 0 {
        cigar.push('L');
        j -= 1;
    }

    cigar.reverse();

    let query_name = String::from("Temp");
    let seq_length = dpm[0].len() - 1;
    let query_start = 0;
    let query_end = dpm[0].len() - 2;
    let strand = '+';
    handle_id_alignment.dedup();
    handle_id_alignment.reverse();
    let path: Vec<usize> = handle_id_alignment.iter().map(|id| *id as usize).collect();
    //path length already set
    let path_start = get_node_offset(handles_nodes_id, if i == 0 { i } else { i + 1 }) as usize; // first letter used in first node of alignment
    let path_end = get_node_offset(handles_nodes_id, ending_node) as usize; // last letter used in last node of alignment

    let align_block_length = "*"; // to set
    let mapping_quality = "*"; // to set
    let comments = format!("{}, best path: {}", build_cigar(&cigar), best_path);
    let gaf_output = GAFStruct::build_gaf_struct(
        query_name,
        seq_length,
        query_start,
        query_end,
        strand,
        path,
        path_length,
        path_start,
        path_end,
        0,
        String::from(align_block_length),
        String::from(mapping_quality),
        comments,
    );
    gaf_output
}

pub fn build_alignment_gap(
    dpm: &Vec<Vec<Vec<i32>>>,
    x: &Vec<Vec<Vec<i32>>>,
    y: &Vec<Vec<Vec<i32>>>,
    alphas: &Vec<usize>,
    best_path: usize,
    pred_hash: &PredHash,
    nwp: &BitVec,
) -> String {
    let mut cigar = Vec::new();
    let mut i = 0;
    let ending_nodes = pred_hash.get_preds_and_paths(dpm.len() - 1);
    for (node, paths) in ending_nodes.iter() {
        if paths[best_path] {
            i = *node;
        }
    }
    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
        let curr_score = if alphas[i] == best_path {
            dpm[i][j][best_path]
        } else {
            dpm[i][j][best_path] + dpm[i][j][alphas[i]]
        };
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j - 1][best_path]
                } else {
                    dpm[i - 1][j - 1][best_path] + dpm[i - 1][j - 1][alphas[i - 1]]
                },
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j][best_path]
                } else {
                    dpm[i - 1][j][best_path] + dpm[i - 1][j][alphas[i - 1]]
                },
                if alphas[i] == best_path {
                    dpm[i][j - 1][best_path]
                } else {
                    dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]]
                },
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    if alphas[*pred] == best_path {
                        d = dpm[*pred][j - 1][best_path];
                        u = dpm[*pred][j][best_path];
                    } else {
                        d = dpm[*pred][j - 1][best_path] + dpm[*pred][j - 1][alphas[*pred]];
                        u = dpm[*pred][j][best_path] + dpm[*pred][j][alphas[*pred]];
                    }
                    if alphas[i] == best_path {
                        l = dpm[i][j - 1][best_path];
                    } else {
                        l = dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]];
                    }
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }

            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
        } else if max == u {
            cigar.push('U');
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            while dpm[i][j][best_path] < y[i][j][best_path] {
                cigar.push('U');
                if nwp[i] {
                    let preds = pred_hash.get_preds_and_paths(i);
                    for (pred, paths) in preds.iter() {
                        if paths[best_path] {
                            predecessor = Some(*pred);
                        }
                    }
                } else {
                    predecessor = Some(i - 1);
                }
                i = predecessor.unwrap();
            }
        } else {
            cigar.push('L');
            j -= 1;
            while dpm[i][j][best_path] < x[i][j][best_path] {
                cigar.push('L');
                j -= 1;
            }
        }
    }
    while j > 0 {
        cigar.push('L');
        j -= 1;
    }
    while i > 0 {
        cigar.push('U');
        i -= 1;
    }
    cigar.reverse();
    cigar.pop();
    build_cigar(&cigar)
}

pub fn build_alignment_semiglobal_gap(
    dpm: &Vec<Vec<Vec<i32>>>,
    x: &Vec<Vec<Vec<i32>>>,
    y: &Vec<Vec<Vec<i32>>>,
    alphas: &Vec<usize>,
    best_path: usize,
    pred_hash: &PredHash,
    nwp: &BitVec,
    ending_node: usize,
) -> String {
    let mut cigar = Vec::new();
    let mut i = ending_node;
    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
        let curr_score = if alphas[i] == best_path {
            dpm[i][j][best_path]
        } else {
            dpm[i][j][best_path] + dpm[i][j][alphas[i]]
        };
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j - 1][best_path]
                } else {
                    dpm[i - 1][j - 1][best_path] + dpm[i - 1][j - 1][alphas[i - 1]]
                },
                if alphas[i - 1] == best_path {
                    dpm[i - 1][j][best_path]
                } else {
                    dpm[i - 1][j][best_path] + dpm[i - 1][j][alphas[i - 1]]
                },
                if alphas[i] == best_path {
                    dpm[i][j - 1][best_path]
                } else {
                    dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]]
                },
            )
        } else {
            let preds = pred_hash.get_preds_and_paths(i);
            let (mut d, mut u, mut l) = (0, 0, 0);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    predecessor = Some(*pred);
                    if alphas[*pred] == best_path {
                        d = dpm[*pred][j - 1][best_path];
                        u = dpm[*pred][j][best_path];
                    } else {
                        d = dpm[*pred][j - 1][best_path] + dpm[*pred][j - 1][alphas[*pred]];
                        u = dpm[*pred][j][best_path] + dpm[*pred][j][alphas[*pred]];
                    }
                    if alphas[i] == best_path {
                        l = dpm[i][j - 1][best_path];
                    } else {
                        l = dpm[i][j - 1][best_path] + dpm[i][j - 1][alphas[i]];
                    }
                }
            }
            (d, u, l)
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            if curr_score < d {
                cigar.push('d');
            } else {
                cigar.push('D');
            }

            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            j -= 1;
        } else if max == u {
            cigar.push('U');
            i = if predecessor.is_none() {
                i - 1
            } else {
                predecessor.unwrap()
            };
            while dpm[i][j][best_path] < y[i][j][best_path] {
                cigar.push('U');
                if nwp[i] {
                    let preds = pred_hash.get_preds_and_paths(i);
                    for (pred, paths) in preds.iter() {
                        if paths[best_path] {
                            predecessor = Some(*pred);
                        }
                    }
                } else {
                    predecessor = Some(i - 1);
                }
                i = predecessor.unwrap();
            }
        } else {
            cigar.push('L');
            j -= 1;
            while dpm[i][j][best_path] < x[i][j][best_path] {
                cigar.push('L');
                j -= 1;
            }
        }
    }
    while j > 0 {
        cigar.push('L');
        j -= 1;
    }

    cigar.reverse();
    let mut starting_node = 0;
    while i > 0 {
        if nwp[i] {
            let preds = pred_hash.get_preds_and_paths(i);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    i = *pred;
                }
            }
        } else {
            i -= 1
        }
        starting_node += 1;
    }

    let mut final_node = 0;
    i = ending_node;
    while i > 0 {
        if nwp[i] {
            let preds = pred_hash.get_preds_and_paths(i);
            for (pred, paths) in preds.iter() {
                if paths[best_path] {
                    i = *pred;
                }
            }
        } else {
            i -= 1
        }
        final_node += 1;
    }
    let cigar = build_cigar(&cigar);
    let cigar_output = format!("{}\t({} {})", cigar, starting_node, final_node);
    cigar_output
}

pub fn extract_best_path_matrix(
    best_path: usize,
    dpm: &Vec<Vec<Vec<i32>>>,
    alphas: &Vec<usize>,
) -> Vec<Vec<i32>> {
    let mut bpm: Vec<Vec<i32>> = vec![vec![0; dpm[0].len()]; dpm.len()];
    for i in 0..dpm.len() {
        for j in 0..dpm[0].len() {
            bpm[i][j] = if alphas[i] == best_path {
                dpm[i][j][best_path]
            } else {
                dpm[i][j][best_path] + dpm[i][j][alphas[i]]
            }
        }
    }
    bpm
}

pub fn build_cigar(cigar: &Vec<char>) -> String {
    let mut output_string = String::new();

    let mut d_count = 0;
    let mut u_count = 0;
    let mut l_count = 0;
    let mut mm_count = 0;
    for ch in cigar.iter() {
        match ch {
            'D' => {
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                d_count += 1;
            }

            'U' => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                u_count += 1;
            }
            'd' => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                mm_count += 1;
            }
            _ => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0;
                }
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                l_count += 1;
            }
        }
    }
    if d_count != 0 {
        output_string = format!("{}{}M", output_string, d_count);
    }
    if u_count != 0 {
        output_string = format!("{}{}I", output_string, u_count);
    }
    if l_count != 0 {
        output_string = format!("{}{}D", output_string, l_count);
    }
    if mm_count != 0 {
        output_string = format!("{}{}X", output_string, mm_count);
    }
    output_string
}

fn get_node_offset(nodes_handles: &Vec<u64>, curr_node: usize) -> i32 {
    let handle = nodes_handles[curr_node];
    if handle == 0 {
        0
    } else {
        let mut counter = curr_node;
        let mut offset = 0;
        while nodes_handles[counter - 1] == handle {
            counter -= 1;
            offset += 1;
        }
        offset
    }
}
