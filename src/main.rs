use bstr::BStr;
use bstr::BString;

use recgraph::a_star::a_star_demo;
use recgraph::args_parser::ClArgs;
use recgraph::node_displacement::DisplacementMatrix;
use recgraph::pathwise_alignment_recombination;
use recgraph::pathwise_graph;

use recgraph::score_matrix;
use recgraph::sequences;
use recgraph::utils;

use std::time::Instant;
use std::time::SystemTime;

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let now = SystemTime::now();

    let args = ClArgs::parse();

    let (sequences, ids) = sequences::get_sequences(args.sequence_path);
    let graph = pathwise_graph::read_graph_w_path(&args.graph_path, false);

    let displ_matrix = DisplacementMatrix::new(&graph);

    let score_matrix =
        score_matrix::create_score_matrix(args.match_score, args.mismatch_score, args.gap_ext);

    let (base_rec_cost, multi_rec_cost) = (args.base_rec_cost, args.multi_rec_cost);
    let is_local = args.alignment_mode;
    let start = Instant::now();
    sequences.iter().enumerate().for_each(|(i, seq)| {
        a_star_demo::a_star_demo(seq);
    });
    println!("basic time: {:?}", start.elapsed());

    let start = Instant::now();
    sequences.iter().enumerate().for_each(|(i, seq)| {
        a_star_demo::a_star_demo_chain(seq);
    });
    println!("chain time: {:?}", start.elapsed());

    /*
    let mut gafs = Vec::new();

    let mut ress = Vec::new();

    sequences.iter().enumerate().for_each(|(i, seq)| {
        let res = a_star_demo::a_star_demo(&seq);
        let mut gaf = pathwise_alignment_recombination::exec(
            is_local,
            seq,
            &graph,
            &score_matrix,
            base_rec_cost,
            multi_rec_cost,
            &displ_matrix,
            args.rec_number,
        );
        gaf.query_name = ids[i].to_string();
        gafs.push(gaf.to_string());
        ress.push(res);
    });

    for (i, gaf) in gafs.iter().enumerate() {
        utils::write_gaf(gaf, i, args.out_file.as_str());
        println!("{:?}",ress[i]);
        println!();
        println!()
    }

     */
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
