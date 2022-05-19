use std::{
    cmp::{self, Ordering},
    collections::HashMap,
};

use bit_vec::BitVec;
use handlegraph::{handle::Handle, handlegraph::HandleGraph, hashgraph::HashGraph};
#[inline]
pub fn set_ampl_for_row(
    i: usize,
    p_arr: &[usize],
    r_val: usize,
    best_scoring_pos: &[usize],
    seq_len: usize,
    bta: usize,
    simd_version: bool,
) -> (usize, usize) {
    let ms;
    let me;
    if i == 0 {
        ms = 0;
        me = 0;
    } else if p_arr.is_empty() {
        let pl = best_scoring_pos[i - 1];
        ms = pl + 1;
        me = pl + 1;
    } else {
        let mut pl = 0;
        let mut pr = 0;
        let mut first = true;
        for p in p_arr.iter() {
            let current_best = best_scoring_pos[*p];
            if first {
                pl = current_best;
                pr = current_best;
                first = false;
            }
            if current_best < pl {
                pl = current_best;
            }
            if current_best > pr {
                pr = current_best;
            }
        }
        ms = pl + 1;
        me = pr + 1;
    }
    let tmp_bs = cmp::min(ms as i32, (seq_len as i32 - r_val as i32) - bta as i32);
    let band_start = if tmp_bs < 0 {
        0
    } else {
        cmp::max(0, tmp_bs as usize)
    };
    let band_end = if seq_len > r_val {
        cmp::min(seq_len, cmp::max(me, seq_len - r_val) + bta)
    } else {
        cmp::min(seq_len, me + bta)
    };
    if simd_version {
        set_left_right_x64(band_start, band_end, seq_len)
    } else {
        (band_start, band_end)
    }
}

fn set_left_right_x64(left: usize, right: usize, seq_len: usize) -> (usize, usize) {
    let mut new_right = right;
    let mut new_left = left;
    while (new_right - new_left) % 8 != 0 {
        if (new_right - new_left) % 2 == 0 && new_right < seq_len {
            new_right += 1;
        } else if new_left > 0 {
            new_left -= 1;
        } else {
            break;
        }
    }
    if new_left == 0 {
        while (new_right - 1) % 8 != 0 && new_right < seq_len {
            new_right += 1;
        }
    }
    if new_right == seq_len {
        while (new_right - new_left) % 8 != 0 && new_left > 1 {
            new_left -= 1
        }
    }

    (new_left, new_right)
}

pub fn set_r_values(
    nwp: &bit_vec::BitVec,
    pred_hash: &HashMap<usize, Vec<usize>>,
    lnz_len: usize,
) -> Vec<usize> {
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[lnz_len - 1] = 0;
    for p in pred_hash.get(&(lnz_len - 1)).unwrap() {
        r_values[*p] = 0;
    }
    for i in (1..lnz_len - 1).rev() {
        if r_values[i] == -1 || r_values[i] > r_values[i + 1] + 1 {
            r_values[i] = r_values[i + 1] + 1;
        }
        if nwp[i] {
            for p in pred_hash.get(&i).unwrap() {
                if r_values[*p] == -1 || r_values[*p] > r_values[i] + 1 {
                    r_values[*p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values.iter().map(|x| *x as usize).collect()
}

#[inline]
pub fn get_max_d_u_l(d: i32, u: i32, l: i32) -> (i32, char) {
    match d.cmp(&u) {
        Ordering::Less => match u.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (u, 'U'),
        },
        _ => match d.cmp(&l) {
            Ordering::Less => (l, 'L'),
            _ => (d, 'D'),
        },
    }
}

pub fn output_creation(path: &[Vec<f32>], graph_lnz: &Vec<char>, read_number: usize) {
    let mut row = path.len() - 2;
    let mut col = path[row].len() - 1;
    let mut cigar = String::new();
    let mut output_ok = true;
    let mut graph_sequence = String::new();

    while path[row][col] != 0.0 {
        let val = path[row][col];
        if val == -1f32 {
            println!("larger band needed for a correct output");
            output_ok = false;
            break;
        }
        let val_str = val.to_string();
        let pred_dir = val_str.split('.').collect::<Vec<&str>>();
        let pred = pred_dir[0].parse::<usize>().unwrap();
        let dir = pred_dir[1].parse::<i32>().unwrap();
        match dir {
            1 => {
                cigar.push('M');
                col -= 1;
                graph_sequence.push(graph_lnz[row]);
            }
            2 => {
                cigar.push('I');
                graph_sequence.push(graph_lnz[row]);
            }
            3 => {
                cigar.push('D');
                col -= 1;
            }
            _ => {
                panic!();
            }
        };
        row = pred;
    }
    if output_ok {
        cigar = cigar.chars().rev().collect::<String>();
        graph_sequence = graph_sequence.chars().rev().collect();
        let mut count_m = 0;
        let mut count_i = 0;
        let mut count_d = 0;
        let mut output = String::new();
        for c in cigar.chars() {
            match c {
                'M' => {
                    if count_i > 0 {
                        output = format!("{}{}I", output, count_i);
                        count_i = 0;
                    }
                    if count_d > 0 {
                        output = format!("{}{}D", output, count_d);
                        count_d = 0;
                    }
                    count_m += 1;
                }
                'I' => {
                    if count_m > 0 {
                        output = format!("{}{}M", output, count_m);
                        count_m = 0;
                    }
                    if count_d > 0 {
                        output = format!("{}{}D", output, count_d);
                        count_d = 0;
                    }
                    count_i += 1;
                }
                'D' => {
                    if count_m > 0 {
                        output = format!("{}{}M", output, count_m);
                        count_m = 0;
                    }
                    if count_i > 0 {
                        output = format!("{}{}I", output, count_i);
                        count_i = 0;
                    }
                    count_d += 1;
                }
                _ => {
                    panic!()
                }
            }
        }
        if count_m > 0 {
            output = format!("{}{}M", output, count_m);
        }
        if count_i > 0 {
            output = format!("{}{}I", output, count_i);
        }
        if count_d > 0 {
            output = format!("{}{}D", output, count_d);
        }
        println!("graph:\t{}", graph_sequence);
        println!(">{}\t{}", read_number, output);
        println!();
    }
}
pub fn create_handle_pos_in_lnz(
    nwp: &BitVec,
    file_path: &str,
    amb_mode: bool,
) -> HashMap<usize, String> {
    let sorted_handles = crate::graph::get_sorted_handles(file_path, amb_mode);
    let mut curr_handle_idx = 0;
    let mut handle_of_lnz_pos = HashMap::new();
    for i in 1..nwp.len() - 1 {
        if nwp[i] {
            curr_handle_idx += 1;
        }
        handle_of_lnz_pos.insert(
            i,
            sorted_handles[(curr_handle_idx - 1) as usize]
                .id()
                .to_string(),
        );
    }
    handle_of_lnz_pos.insert(0, String::from("-1"));
    handle_of_lnz_pos
}

pub fn handle_pos_in_lnz_from_hashgraph(
    nwp: &BitVec,
    graph: &HashGraph,
    amb_mode: bool,
) -> HashMap<usize, String> {
    let mut sorted_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_handles.sort();
    if amb_mode {
        sorted_handles.reverse();
        sorted_handles = sorted_handles
            .iter()
            .map(|h| h.flip())
            .collect::<Vec<Handle>>();
    }
    let mut curr_handle_idx = 0;
    let mut handle_of_lnz_pos = HashMap::new();
    for i in 1..nwp.len() - 1 {
        if nwp[i] {
            curr_handle_idx += 1;
        }
        handle_of_lnz_pos.insert(
            i,
            sorted_handles[(curr_handle_idx - 1) as usize]
                .id()
                .to_string(),
        );
    }
    handle_of_lnz_pos.insert(0, String::from("-1"));
    handle_of_lnz_pos
}
