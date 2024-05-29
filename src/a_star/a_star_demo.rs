use bit_vec::BitVec;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use std::time::Instant;

use crate::a_star::{a_star_visit, build_heuristic as new_heuristic, check_ed};
use crate::args_parser::ClArgs;
use crate::new_path_graph::path_graph::{remove_duplicate_paths, PathGraph};
use crate::sequences;

use super::a_star_output::build_gaf;

pub fn a_star_demo_chain() {
    let args = ClArgs::parse();
    let file_path = args.graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let mut graph: HashGraph = HashGraph::from_gfa(&gfa);
    //remove_duplicate_paths(&mut graph);
    let path_graph = PathGraph::from_hash_graph(&graph);
    /*
    println!("LNZ: {:?}", path_graph.lnz.len());
    println!("Paths: {:?}", path_graph.get_node_path(36066));
    path_graph.succ_hash.get_all_paths().iter().for_each(|(handle, paths)| {
        if paths[9] {
            println!("{:?}", handle);
        }
    });
    let mut i = 29897;
    let mut paths = BitVec::from_elem(path_graph.succ_hash.paths_number as usize, false);
    paths.set(path_graph.succ_hash.paths_number as usize - 1, true);
    while i < path_graph.lnz.len() {
        if &paths !=  path_graph.get_node_path(i as u32) {
            println!("{i}\t{:?}", path_graph.get_node_path(i as u32));
        }
        i += 1;
    }
    */

    let (sequences, _) = sequences::get_sequences(args.sequence_path);
    let chunk_size = ClArgs::parse().seed_len;
    let start = Instant::now();
    let indexes = path_graph.get_indexes();
    let mut outs = Vec::new();
    let init = peak_mem_usage().unwrap();
    sequences.iter().for_each(|seq| {
        let istant = Instant::now();
        let mut heuristic = new_heuristic::build_heuristic(
            &indexes,
            seq,
            chunk_size as usize,
            args.base_rec_cost as usize,
        );
        let (end_pos, mut alignment_graph) = a_star_visit::exec(
            seq,
            &mut heuristic,
            &path_graph,
            args.alignment_mode,
            args.base_rec_cost as u32,
        );

        outs.push((
            build_gaf(
                &mut alignment_graph,
                &end_pos,
                &path_graph,
                seq,
                args.alignment_mode,
            ),
            istant.elapsed(),
            seq.len(),
        ));
    });
    outs.iter()
        .for_each(|out| println!("{}\t{}\t{}", out.0, out.1.as_millis(), out.2));
    println!("Approx time: {:?}", start.elapsed());
    let mem = peak_mem_usage().unwrap();
    println!("Init memory usage: {} B", init);
    println!("Peak memory usage: {} B", mem);
    check_ed::test(&path_graph, &sequences);
}

#[cfg(target_os = "linux")]
fn peak_mem_usage() -> Result<usize, &'static str> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize * 1024),
            _ => Err("Error getting memory usage"),
        }
    }
}
