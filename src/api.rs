use std::collections::HashMap;

use crate::{
    gap_local_poa, gap_mk_abpoa, global_mk_abpoa, graph, local_poa, matrix, sequences, utils,
};
use handlegraph::hashgraph::HashGraph;

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
    let hofp = utils::handle_pos_in_lnz_from_hashgraph(&lnz_graph.nwp, &graph, false);

    unsafe {
        global_mk_abpoa::exec_simd(
            &read_for_alignment,
            sequence_name.unwrap_or(("no_name", 1)),
            &lnz_graph,
            &score_matrix_f32,
            bases_to_add,
            false,
            &hofp,
            &r_values,
        );
    }
}

pub fn align_global_gap(
    read: &String,
    graph: &HashGraph,
    sequence_name: Option<(&str, usize)>,
    score_matrix: Option<HashMap<(char, char), i32>>,
    bta: Option<usize>,
    o: Option<i32>,
    e: Option<i32>,
) {
    let read_for_alignment = sequences::build_align_string(read);
    let lnz_graph = graph::create_graph_struct(graph, false);
    let score_matrix_i32 = score_matrix.unwrap_or(matrix::create_score_matrix_match_mis(2, -4));
    let bases_to_add = bta.unwrap_or(read_for_alignment.len() * 0.1 as usize);

    let hofp = utils::handle_pos_in_lnz_from_hashgraph(&lnz_graph.nwp, &graph, false);

    gap_mk_abpoa::exec(
        &read_for_alignment,
        sequence_name.unwrap_or(("no_name", 1)),
        &lnz_graph,
        &score_matrix_i32,
        o.unwrap_or(-10),
        e.unwrap_or(-6),
        bases_to_add,
        false,
        &hofp,
    );
}

pub fn align_local_no_gap(
    read: &String,
    graph: &HashGraph,
    sequence_name: Option<(&str, usize)>,
    score_matrix: Option<HashMap<(char, char), f32>>,
) {
    let read_for_alignment = sequences::build_align_string(read);
    let lnz_graph = graph::create_graph_struct(graph, false);
    let score_matrix_f32 =
        score_matrix.unwrap_or(matrix::create_score_matrix_match_mis_f32(2f32, -4f32));
    let hofp = utils::handle_pos_in_lnz_from_hashgraph(&lnz_graph.nwp, &graph, false);

    unsafe {
        local_poa::exec_simd(
            &read_for_alignment,
            sequence_name.unwrap_or(("no_name", 1)),
            &lnz_graph,
            &score_matrix_f32,
            false,
            &hofp,
        );
    }
}

pub fn align_local_gap(
    read: &String,
    graph: &HashGraph,
    sequence_name: Option<(&str, usize)>,
    score_matrix: Option<HashMap<(char, char), i32>>,
    o: Option<i32>,
    e: Option<i32>,
) {
    let read_for_alignment = sequences::build_align_string(read);
    let lnz_graph = graph::create_graph_struct(graph, false);
    let score_matrix_i32 = score_matrix.unwrap_or(matrix::create_score_matrix_match_mis(2, -4));

    let hofp = utils::handle_pos_in_lnz_from_hashgraph(&lnz_graph.nwp, &graph, false);

    gap_local_poa::exec(
        &read_for_alignment,
        sequence_name.unwrap_or(("no_name", 1)),
        &lnz_graph,
        &score_matrix_i32,
        o.unwrap_or(-10),
        e.unwrap_or(-6),
        false,
        &hofp,
    );
}
