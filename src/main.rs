mod ab_gap_global_alignment;
mod ab_global_alignment;
mod args_parser;
mod banded_mk_poa;
mod basic_output;
mod gap_partial_order_alignment;
mod global_alignment_affine_gap;
mod graph;
mod local_alignment;
mod matrix;
mod partial_order_alignment_global;
mod sequences;
use std::cmp;

fn main() {
    let sequences = sequences::get_sequences();
    let mut s1: Vec<char> = sequences[10].chars().collect();
    let mut s2: Vec<char> = sequences[11].chars().collect();
    s1.insert(0, '$');
    s2.insert(0, '$');

    let score_matrix = matrix::create_score_matrix();
    let align_mode = args_parser::get_align_mode();

    // set band ampl
    let (b, f) = args_parser::get_b_f();
    let bases_to_add = (b + f * s1.len() as f32) as usize;
    let ampl = match s1.len() < s2.len() {
        true => s2.len() - s1.len() + bases_to_add,
        _ => s1.len() - s2.len() + bases_to_add,
    };

    match align_mode {
        0 => ab_global_alignment::exec(&s1, &s2, &score_matrix, cmp::max(ampl * 2 + 1, 3)),
        1 => {
            /*
            if ampl * 2 + 1 > cmp::min(s1.len(), s2.len()) {
                local_alignment::exec(&s1, &s2, &score_matrix)
            } else {
                ab_local_alignment::exec(&s1, &s2, &score_matrix, cmp::max(ampl * 2 + 1, 3))
            }
            */
            local_alignment::exec(&s1, &s2, &score_matrix)
        }
        2 => println!("edit distance removed"),
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            global_alignment_affine_gap::exec(&s1, &s2, &score_matrix, g_open, g_ext);
            ab_gap_global_alignment::exec(
                &s1,
                &s2,
                &score_matrix,
                cmp::max(ampl * 2 + 1, 3),
                g_open,
                g_ext,
            );
        }
        4 => {
            let mut sequence: Vec<char> = sequences[13].chars().collect();
            sequence.insert(0, '$');
            let graph_path = args_parser::get_graph_path();
            let linearization = graph::get_linearization(&graph_path);

            partial_order_alignment_global::exec(&sequence, &linearization, &score_matrix);
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            gap_partial_order_alignment::exec(
                &sequence,
                &linearization,
                &score_matrix,
                g_open,
                g_ext,
            );
            banded_mk_poa::exec(&sequence, &linearization, &score_matrix, 5);
        }
        _ => panic!("alignment mode must be 0, 1, 2 or 3"),
    }
}
