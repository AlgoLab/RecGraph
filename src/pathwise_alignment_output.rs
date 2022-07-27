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
            i -= 1;
            j -= 1;
        } else if max == u {
            cigar.push('U');
            i -= 1;
        } else {
            cigar.push('L');
            j -= 1;
        }
    }
    cigar.reverse();
    build_cigar(&cigar)
}

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
