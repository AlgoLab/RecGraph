use std::time::Duration;

use crate::{
    a_star::{
        a_star_output::{build_gaf, Gaf},
        *,
    },
    args_parser::EstimateFunction,
    new_path_graph::path_graph::PathGraph,
};
use ahash::AHashMap;
use bstr::BString;
use handlegraph::hashgraph;

pub fn astar_align(
    read: &BString,
    graph: &PathGraph,
    indexes: &Vec<(
        lt_fm_index::LtFmIndex<u32, lt_fm_index::blocks::Block3<u128>>,
        Vec<u32>,
    )>,
    est_function: i32,
    is_local: bool,
    rec_cost: i32,
    max_rec: i32,
    seed_length: usize,
    sequence_name: Option<&String>,
    gap_open: u16,
) -> String {
    let mut heuristic = match est_function {
        0 => chain_heur::build_heuristic(
            &indexes,
            read,
            seed_length,
            rec_cost as usize,
            graph,
            max_rec > 0,
        ),
        1 => seed_heuristic::build_heuristic(
            &indexes,
            read,
            seed_length,
            rec_cost as usize,
            graph,
            max_rec > 0,
        ),
        _ => fast_heuristic::build_heuristic(
            &indexes,
            read,
            seed_length,
            rec_cost as usize,
            graph,
            max_rec > 0,
        ),
    };
    let (end_pos, mut alignment_graph) = a_star_visit::exec(
        read,
        &mut heuristic,
        &graph,
        is_local,
        rec_cost as u16,
        max_rec as u32,
        match est_function {
            0 => &EstimateFunction::Chaining,
            1 => &EstimateFunction::Seeding,
            _ => &EstimateFunction::Fast,
        },
        gap_open as u16,
    );
    build_gaf(
        &mut alignment_graph,
        &end_pos,
        graph,
        read,
        is_local,
        &Vec::new(),
        0,
        indexes,
        &sequence_name
            .map(|s| BString::from(s.as_str()))
            .unwrap_or_else(|| BString::from("")),
        Duration::new(0, 0),
        &AHashMap::new(),
        0.0,
    )
    .to_string()
}

pub fn convert_hash_graph_to_path_graph(graph: &hashgraph::HashGraph) -> PathGraph {
    PathGraph::from_hash_graph(graph)
}

pub fn build_indexes(
    graph: &PathGraph,
) -> Vec<(
    lt_fm_index::LtFmIndex<u32, lt_fm_index::blocks::Block3<u128>>,
    Vec<u32>,
)> {
    graph.get_indexes()
}

pub fn convert_gaf_to_bstring(gaf: &Gaf) -> BString {
    let gaf_str = gaf.to_string();
    BString::from(gaf_str)
}
