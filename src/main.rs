use recgraph::args_parser::ClArgs;
use recgraph::pathwise_alignment_recombination;
use recgraph::pathwise_graph;
use recgraph::pathwise_graph::nodes_displacement_matrix;
use recgraph::score_matrix;
use recgraph::sequences;
use recgraph::utils;

use std::time::SystemTime;

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let now = SystemTime::now();

    let args = ClArgs::parse();

    let (sequences, ids) = sequences::get_sequences(args.sequence_path);
    let graph = pathwise_graph::read_graph_w_path(&args.graph_path, false);
    let rev_graph = pathwise_graph::create_reverse_path_graph(&graph);

    let displ_matrix = nodes_displacement_matrix(&graph, &rev_graph);

    let score_matrix =
        score_matrix::create_score_matrix(args.match_score, args.mismatch_score, args.gap_ext);

    let (base_rec_cost, multi_rec_cost) = (args.base_rec_cost, args.multi_rec_cost);
    let is_local = args.alignment_mode;
    let mut gafs = Vec::new();

    sequences.iter().enumerate().for_each(|(i, seq)| {
        let mut gaf = pathwise_alignment_recombination::exec(
            is_local,
            seq,
            &graph,
            &rev_graph,
            &score_matrix,
            base_rec_cost,
            multi_rec_cost,
            &displ_matrix,
            args.rec_number,
        );
        gaf.query_name = ids[i].to_string();
        gafs.push(gaf.to_string());
    });

    for (i, gaf) in gafs.iter().enumerate() {
        utils::write_gaf(gaf, i, args.out_file.as_str());
    }
    match now.elapsed() {
        Ok(elapsed) => {
            // it prints '2'
            eprintln!("Done in {:?}.", elapsed);
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {e:?}");
        }
    }
}
