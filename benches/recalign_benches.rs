use std::env;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use recalign::{a_star::a_star_demo, new_path_graph::path_graph::{remove_duplicate_paths, PathGraph}, sequences};

pub fn criterion_benchmark(c: &mut Criterion) {
    let seqs_path = String::from("example/reads.fa");
    let (sequences, _) = sequences::get_sequences(seqs_path);
    let graph_path = "example/graph.gfa";
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(graph_path).unwrap();
    let mut graph: HashGraph = HashGraph::from_gfa(&gfa);
    remove_duplicate_paths(&mut graph);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let indexes = path_graph.get_indexes();
    c.bench_function("a_star_demo_chain", |b| b.iter(|| a_star_demo::alignment_bench(black_box(&sequences), black_box(&path_graph), 10, 4, true, 2,black_box( &indexes))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);