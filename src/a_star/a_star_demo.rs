use std::time::{Duration, Instant};

use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;

use crate::a_star::{a_star_visit, approx_matching};
use crate::args_parser::ClArgs;
use crate::new_path_graph::path_graph::PathGraph;
use crate::sequences;

use super::a_star_output::build_gaf;

pub fn a_star_demo_chain() {
    let args = ClArgs::parse();

    let file_path = args.graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let (sequences, _) = sequences::get_sequences(args.sequence_path);

    let chunk_size = ClArgs::parse().seed_len;
    let start = Instant::now();
    let linearized_paths = approx_matching::get_linearized_paths_and_handles(&graph);
    let mut outs = Vec::new();
    let mut count = Duration::new(0, 0);
    sequences.iter().for_each(|seq| {
        let istant = Instant::now();
        let crumbs = approx_matching::build_heuristic(&linearized_paths, seq, chunk_size as usize);
        count += Instant::now() - istant;
        let (end_pos, mut alignment_graph) = a_star_visit::exec(seq, &crumbs, &path_graph);

        outs.push(build_gaf(&mut alignment_graph, &end_pos, &path_graph, seq));
    });
    outs.iter().for_each(|out| println!("{}", out));
    println!("Heu: {:?}", count);
    println!("Approx time: {:?}", start.elapsed());
}
