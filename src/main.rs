mod args_parser;
mod bitfield_path;
mod gaf_output;
mod gap_local_poa;
mod gap_mk_abpoa;
mod global_mk_abpoa;
mod graph;
mod local_poa;
mod matrix;
mod sequences;
fn main() {
    // get sequence
    let (sequences, seq_names) = sequences::get_sequences();
    let seq: &Vec<char> = &sequences[0];

    //get graph
    let graph_path = args_parser::get_graph_path();
    let graph_struct = graph::read_graph(&graph_path, false);

    //get score matrix
    let score_matrix = matrix::create_score_matrix();

    //get alignment option
    let align_mode = args_parser::get_align_mode();
    let amb_strand = args_parser::get_amb_strand_mode();
    let (b, f) = args_parser::get_b_f();
    let bases_to_add = (b + f * seq.len() as f32) as usize;

    match align_mode {
        //global alignment
        0 => {
            for (i, seq) in sequences.iter().enumerate() {
                let align_score = global_mk_abpoa::exec(
                    seq,
                    &seq_names[i],
                    &graph_struct,
                    &score_matrix,
                    bases_to_add,
                    &graph_path,
                    false,
                );
                if amb_strand && align_score < 0 {
                    let rev_graph_struct = graph::read_graph(&graph_path, true);
                    global_mk_abpoa::exec(
                        seq,
                        &seq_names[i],
                        &rev_graph_struct,
                        &score_matrix,
                        bases_to_add,
                        &graph_path,
                        true,
                    );
                }
            }
        }
        //local alignment
        1 => {
            for (i, seq) in sequences.iter().enumerate() {
                local_poa::exec(
                    seq,
                    &seq_names[i],
                    &graph_struct,
                    &score_matrix,
                    &graph_path,
                    false,
                );
            }
        }
        //affine gap global alignment
        2 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let align_score = gap_mk_abpoa::exec(
                    seq,
                    &seq_names[i],
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    bases_to_add,
                    &graph_path,
                    false,
                );

                if amb_strand && align_score < 0 {
                    let rev_graph_struct = graph::read_graph(&graph_path, true);
                    gap_mk_abpoa::exec(
                        seq,
                        &seq_names[i],
                        &rev_graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        bases_to_add,
                        &graph_path,
                        true,
                    );
                }
            }
        }
        //affine gap global alignment
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            for (i, seq) in sequences.iter().enumerate() {
                let align_score = gap_local_poa::exec(
                    seq,
                    &seq_names[i],
                    &graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    &graph_path,
                    false,
                );
                if amb_strand && align_score < 0 {
                    let rev_graph_struct = graph::read_graph(&graph_path, true);
                    gap_local_poa::exec(
                        seq,
                        &seq_names[i],
                        &rev_graph_struct,
                        &score_matrix,
                        g_open,
                        g_ext,
                        &graph_path,
                        true,
                    );
                }
            }
        }
        _ => {
            panic!("alignment mode must be 0, 1, 2 or 3");
        }
    }
}
