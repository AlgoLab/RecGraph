use std::collections::HashMap;
use std::time::Instant;

use rspoa::args_parser;
use rspoa::gaf_output;
use rspoa::gap_local_poa;
use rspoa::gap_mk_abpoa;
use rspoa::global_mk_abpoa;
use rspoa::graph;
use rspoa::local_poa;
use rspoa::matrix;
use rspoa::pathwise_alignment;
use rspoa::sequences;
use rspoa::simd_poa_ed;
fn main() {
    // get sequence
    let (sequences, seq_names) = sequences::get_sequences();

    //get graph
    let graph_path = args_parser::get_graph_path();
    let graph_struct = graph::read_graph(&graph_path, false);

    //get score matrix
    let score_matrix = matrix::create_score_matrix();

    //get alignment option
    let align_mode = args_parser::get_align_mode();
    let amb_strand = args_parser::get_amb_strand_mode();
    let (b, f) = args_parser::get_b_f();

    //get handle position for output
    let hofp_forward = gaf_output::create_handle_pos_in_lnz(&graph_struct.nwp, &graph_path, false);
    let mut hofp_reverse = HashMap::new();
    unsafe {
        let start1 = Instant::now();
        let simd_res = simd_poa_ed::exec(
            &sequences[0].iter().map(|c| *c as u8).collect::<Vec<u8>>(),
            &graph_struct,
        );
        let i1 = start1.elapsed();
        let start2 = Instant::now();
        let not_simd_res = simd_poa_ed::exec_no_simd(
            &sequences[0].iter().map(|c| *c as u8).collect::<Vec<u8>>(),
            &graph_struct,
        );
        let i2 = start2.elapsed();
        println!("Simd score: {simd_res} No simd score: {not_simd_res}");
        println!("Simd: {} No Simd: {}", i1.as_micros(), i2.as_micros());
    }

    match align_mode {
        //global alignment
        0 => {
            for (i, seq) in sequences.iter().enumerate() {
                let bases_to_add = (b + f * seq.len() as f32) as usize;
                let align_score = global_mk_abpoa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    bases_to_add,
                    false,
                    &hofp_forward,
                );
                if amb_strand && align_score < 0 {
                    if hofp_reverse.is_empty() {
                        hofp_reverse = gaf_output::create_handle_pos_in_lnz(
                            &graph_struct.nwp,
                            &graph_path,
                            true,
                        );
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    global_mk_abpoa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        bases_to_add,
                        true,
                        &hofp_reverse,
                    );
                }
            }
        }
        //local alignment
        1 => {
            for (i, seq) in sequences.iter().enumerate() {
                let align_score = local_poa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    false,
                    &hofp_forward,
                );
                if align_score < 0 && amb_strand {
                    if hofp_reverse.is_empty() {
                        hofp_reverse = gaf_output::create_handle_pos_in_lnz(
                            &graph_struct.nwp,
                            &graph_path,
                            true,
                        );
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    local_poa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        true,
                        &hofp_reverse,
                    );
                }
            }
        }
        //affine gap global alignment
        2 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();

            for (i, seq) in sequences.iter().enumerate() {
                let bases_to_add = (b + f * seq.len() as f32) as usize;
                let align_score = gap_mk_abpoa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    bases_to_add,
                    false,
                    &hofp_forward,
                );

                if amb_strand && align_score < 0 {
                    if hofp_reverse.is_empty() {
                        hofp_reverse = gaf_output::create_handle_pos_in_lnz(
                            &graph_struct.nwp,
                            &graph_path,
                            true,
                        );
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    gap_mk_abpoa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        bases_to_add,
                        true,
                        &hofp_reverse,
                    );
                }
            }
        }
        //affine gap local alignment
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let align_score = gap_local_poa::exec(
                    seq,
                    (&seq_names[i], i + 1),
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    false,
                    &hofp_forward,
                );
                if amb_strand && align_score < 0 {
                    if hofp_reverse.is_empty() {
                        hofp_reverse = gaf_output::create_handle_pos_in_lnz(
                            &graph_struct.nwp,
                            &graph_path,
                            true,
                        );
                    }
                    let rev_seq = sequences::rev_and_compl(seq);
                    gap_local_poa::exec(
                        &rev_seq,
                        (&seq_names[i], i + 1),
                        &graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        false,
                        &hofp_reverse,
                    );
                }
            }
        }
        4 => {
            let path_node = graph::create_nodes_paths(&graph_path);
            pathwise_alignment::exec(&sequences[4], &graph_struct, &path_node, &score_matrix, 3);
        }
        _ => {
            panic!("alignment mode must be 0, 1, 2 or 3");
        }
    }
}
