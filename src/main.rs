mod args_parser;
mod banded_mk_poa;
mod basic_output;
mod bitfield_path;
mod gap_local_poa;
mod gap_mk_abpoa;
mod graph;
mod local_poa;
mod matrix;
mod sequences;
fn main() {
    bitfield_path::test1();
    bitfield_path::example_path();
    // get sequence
    let sequences = sequences::get_sequences();
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
            let ampl = match seq.len() < graph_struct.lnz.len() {
                true => graph_struct.lnz.len() - 1 - seq.len(),
                _ => seq.len() - graph_struct.lnz.len() + 1,
            };
            let align_score = banded_mk_poa::exec(seq, &graph_struct, &score_matrix, ampl * 2);
            if amb_strand && align_score < 0 {
                let rev_graph_struct = graph::read_graph(&graph_path, true);
                banded_mk_poa::exec(seq, &rev_graph_struct, &score_matrix, ampl * 2);
            }
        }
        //local alignment
        1 => {
            local_poa::exec(seq, &graph_struct, &score_matrix);
        }
        //affine gap global alignment
        2 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            let align_score = gap_mk_abpoa::exec(
                seq,
                &graph_struct,
                &score_matrix,
                g_open,
                g_ext,
                bases_to_add,
            );

            if amb_strand && align_score < 0 {
                let rev_graph_struct = graph::read_graph(&graph_path, true);
                gap_mk_abpoa::exec(
                    seq,
                    &rev_graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    bases_to_add,
                );
            }
        }
        //affine gap global alignment
        3 => {
            let (g_open, g_ext) = args_parser::get_gap_open_gap_ext();
            let align_score = gap_local_poa::exec(seq, &graph_struct, &score_matrix, g_open, g_ext);

            if amb_strand && align_score < 0 {
                let rev_graph_struct = graph::read_graph(&graph_path, true);
                gap_mk_abpoa::exec(
                    seq,
                    &rev_graph_struct,
                    &score_matrix,
                    g_open,
                    g_ext,
                    bases_to_add,
                );
            }
        }
        _ => {
            panic!("alignment mode must be 0, 1, 2 or 3");
        }
    }
}
