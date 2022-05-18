use std::collections::HashMap;

use crate::{
    gaf_output, gap_local_poa, gap_mk_abpoa, global_mk_abpoa, graph, local_poa, matrix, sequences,
    utils,
};
use handlegraph::hashgraph::HashGraph;

/*
TODO: alignment algorithm as
    handlegraph
    read as string
    alignment mode as i32
    Option(other par)
*/
pub fn align_global_no_gap(
    read: &String,
    graph: &HashGraph,
    sequence_name: Option<(&str, usize)>,
    score_matrix: Option<HashMap<(char, char), f32>>,
    bta: Option<usize>,
) {
    let read_for_alignment = sequences::build_align_string(read);
    let lnz_graph = graph::create_graph_struct(graph, false);
    let score_matrix_f32 =
        score_matrix.unwrap_or(matrix::create_score_matrix_match_mis_f32(2f32, -4f32));
    let bases_to_add = bta.unwrap_or(read_for_alignment.len() * 0.1 as usize);

    let r_values = utils::set_r_values(&lnz_graph.nwp, &lnz_graph.pred_hash, lnz_graph.lnz.len());
    //TODO: modify create_handle_pos in order to work with lnz graph
    //let hofp_forward = gaf_output::create_handle_pos_in_lnz(&lnz_graph.nwp, &graph_path, false);
    /*
    unsafe {
        global_mk_abpoa::exec_simd(
            &read_for_alignment,
            (&seq_names[i], i + 1),
            &lnz_graph,
            &score_matrix_f32,
            bases_to_add,
            false,
            &hofp_forward,
            &r_values,
        );
    }
    */
}
