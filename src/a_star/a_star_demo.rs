use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use rayon::prelude::*;
use std::time::{Duration, Instant};

use crate::a_star::{a_star_output, a_star_visit, build_heuristic as new_heuristic, check_ed};
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
    let original_path_ids = remove_duplicate_paths(&mut graph);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let (sequences, names) = sequences::get_sequences(args.sequence_path);
    let chunk_size = ClArgs::parse().seed_len;
    let indexes = path_graph.get_indexes();
    let mut outs = Vec::new();
    //let mut explored_pos_vec = Vec::new();
    let mut heur_tot_time = Duration::new(0, 0);
    let mut explore_tot_time = Duration::new(0, 0);
    sequences.iter().zip(names).for_each(|(seq, name)| {
        let istant = Instant::now();
        let mut heuristic = new_heuristic::build_heuristic(
            &indexes,
            seq,
            chunk_size as usize,
            args.base_rec_cost as usize,
            &path_graph,
        );
        let matches_in_path: Vec<_> = heuristic
            .1
            .par_iter()
            .map(|m| a_star_output::get_matches_end_in_path(m))
            .collect();
        let heu_time = istant.elapsed();
        heur_tot_time += heu_time;

        let explore_start = Instant::now();
        let (end_pos, mut alignment_graph /* , explored_pos*/) = a_star_visit::exec(
            seq,
            &mut heuristic,
            &path_graph,
            args.alignment_mode,
            args.base_rec_cost as u16,
            args.max_rec as u32,
        );
        let explore_time = explore_start.elapsed();
        explore_tot_time += explore_time;

        //explored_pos_vec.push(explored_pos);
        let graph_size = (path_graph.lnz.len() * seq.len()) as f32;
        let explored_cells = (alignment_graph.len() as f32 / graph_size) * 100.0;

        outs.push(build_gaf(
            &mut alignment_graph,
            &end_pos,
            &path_graph,
            seq,
            args.alignment_mode,
            &matches_in_path,
            chunk_size as usize,
            &indexes,
            &name,
            explore_time + heu_time,
            &original_path_ids,
            explored_cells,
        ));
    });
    //a_star_output::save_coords(&args.out_file, &explored_pos_vec);
    outs.iter().for_each(|out| println!("{}", out.to_string()));
    let mem = peak_mem_usage().unwrap();

    eprintln!("Peak memory (Byte)\t{}", mem);
    eprintln!("Heuristic tot time\t{:?}", heur_tot_time);
    eprintln!("Explore tot time\t{:?}", explore_tot_time);
    /*
    if args.alignment_mode {
        check_ed::semiglobal_test(&path_graph, &sequences)
    } else {
        check_ed::test(&path_graph, &sequences)
    }
    */
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
