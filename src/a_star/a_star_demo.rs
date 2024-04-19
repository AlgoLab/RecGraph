use std::time::Instant;

use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;

use crate::a_star::seeding_heurisitc;
use crate::args_parser::ClArgs;
use crate::new_path_graph::path_graph::PathGraph;
use crate::sequences;

use super::a_star_output::build_gaf;
use super::a_star_visit;
use super::seeding_heurisitc::get_chaining_sh;

pub fn a_star_demo_chain() {
    let start = Instant::now();
    let args = ClArgs::parse();

    let file_path = args.graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let (sequences, _) = sequences::get_sequences(args.sequence_path);

    let chunk_size = ClArgs::parse().seed_len;
    let indexes = seeding_heurisitc::get_fm_index(&graph);

    sequences.iter().for_each(|seq| {
        let crumbs = get_chaining_sh(
            seq,
            &graph,
            chunk_size as usize,
            &indexes,
            args.base_rec_cost as usize,
        );
        let (end_pos, mut alignment_graph) = a_star_visit::exec(seq, crumbs, &path_graph);

        build_gaf(&mut alignment_graph, &end_pos, &path_graph, seq);
    });
    println!("chain time: {:?}", start.elapsed());
}
