use bit_vec::BitVec;

use crate::pathwise_graph::PredHash;

pub fn build_alignment(
    dpm: &Vec<Vec<Vec<i32>>>,
    alphas: &Vec<usize>,
    best_path: usize,
    pred_hash: &PredHash,
    nwp: &BitVec,
) -> String {
    let mut cigar = Vec::new();
    let mut i = dpm.len() - 1;
    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
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
            cigar.push('D');
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
        cigar.push('U');
        i -= 1;
    }
    cigar.reverse();
    cigar.pop();
    build_cigar(&cigar)
}

pub fn build_alignment_semiglobal(
    dpm: &Vec<Vec<Vec<i32>>>,
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
            cigar.push('D');
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

    build_cigar(&cigar)
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
    let mut i = dpm.len() - 1;
    let mut j = dpm[i].len() - 1;
    while i != 0 && j != 0 {
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
            cigar.push('D');
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
/*
pub fn build_best_path_alignment(
    bpm: &Vec<Vec<i32>>,
    pred_hash: &PredHash,
    nwp: &BitVec,
    best_path: usize,
) -> String {
    let mut cigar = Vec::new();
    let mut i = bpm.len() - 1;
    let mut j = bpm[i].len() - 1;

    while i != 0 && j != 0 {
        let mut predecessor = None;
        let (d, u, l) = if !nwp[i] {
            (bpm[i - 1][j - 1], bpm[i - 1][j], bpm[i][j - 1])
        } else {
            let preds = pred_hash.get_preds_and_paths(i);

            let mut pred = 0;
            for p in preds.iter() {
                if p.1[best_path] {
                    pred = p.0
                }
            }
            predecessor = Some(pred);
            (bpm[pred][j - 1], bpm[pred][j], bpm[i][j - 1])
        };
        let max = *[d, u, l].iter().max().unwrap();
        if max == d {
            cigar.push('D');
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
        cigar.push('U');
        i = if !nwp[i] {
            i - 1
        } else {
            let preds = pred_hash.get_preds_and_paths(i);

            let mut pred = 0;
            for p in preds.iter() {
                if p.1[best_path] {
                    pred = p.0
                }
            }
            pred
        };
    }
    cigar.reverse();

    build_cigar(&cigar)
}
*/
fn build_cigar(cigar: &Vec<char>) -> String {
    let mut output_string = String::new();

    let mut d_count = 0;
    let mut u_count = 0;
    let mut l_count = 0;

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
                u_count += 1;
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
    output_string
}
