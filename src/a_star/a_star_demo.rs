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
    let linearized_paths = approx_matching::get_linearized_paths(&graph);
    let mut outs = Vec::new();
    let mut duration = Duration::new(0, 0);
    sequences.iter().for_each(|seq| {
        let istant = Instant::now();
        let crumbs = approx_matching::build_heuristic(
            &linearized_paths,
            seq,
            chunk_size as usize,
            args.base_rec_cost as usize,
            args.mex_err_seed,
        );
        duration += istant.elapsed();
        let (end_pos, mut alignment_graph) = a_star_visit::exec(
            seq,
            &crumbs,
            &path_graph,
            args.alignment_mode,
            args.base_rec_cost as u32,
        );

        outs.push((build_gaf(
            &mut alignment_graph,
            &end_pos,
            &path_graph,
            seq,
            args.alignment_mode,
        ), istant.elapsed(), seq.len()));
    });
    println!("seeding: {:?}", duration);
    outs.iter().for_each(|out| println!("{}\t{}\t{}", out.0, out.1.as_millis(), out.2));
    println!("Approx time: {:?}", start.elapsed());
}
