use std::time::Instant;

use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;

use crate::args_parser::ClArgs;
use crate::pathwise_graph::create_path_graph;
use crate::sequences;

use super::a_star_output::build_gaf;
use super::a_star_visit;
use super::matches;
use super::matches::get_base_sh;
use super::seeding_heurisitc::get_chaining_sh;

pub fn a_star_demo() {
    let start = Instant::now();
    let args = ClArgs::parse();

    let file_path = args.graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = create_path_graph(&graph, false);

    let (sequences, _) = sequences::get_sequences(args.sequence_path);

    let chunk_size = ClArgs::parse().seed_len;
    let indexes = matches::get_fm_index(&graph);

    sequences.iter().for_each(|seq| {
        let crumbs = get_base_sh(seq, &graph, chunk_size as usize, &indexes);
        let (end_pos, mut alignment_graph) = a_star_visit::exec(seq, crumbs, &path_graph);

        build_gaf(&mut alignment_graph, &end_pos, &path_graph, seq);
    });
    println!("base time: {:?}", start.elapsed());
}

pub fn a_star_demo_chain() {
    let start = Instant::now();
    let args = ClArgs::parse();

    let file_path = args.graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = create_path_graph(&graph, false);

    let (sequences, _) = sequences::get_sequences(args.sequence_path);

    let chunk_size = ClArgs::parse().seed_len;
    let indexes = matches::get_fm_index(&graph);

    sequences.iter().for_each(|seq| {
        let crumbs = get_chaining_sh(seq, &graph, chunk_size as usize, &indexes);
        crumbs.iter().enumerate().for_each(|(path, crumb)| {
            println!("{path} {:?}", crumb);
        });
        let (end_pos, mut alignment_graph) = a_star_visit::exec(seq, crumbs, &path_graph);

        build_gaf(&mut alignment_graph, &end_pos, &path_graph, seq);
    });
    println!("chain time: {:?}", start.elapsed());
}
