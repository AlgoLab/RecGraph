use criterion::{black_box, criterion_group, criterion_main, Criterion};
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use recgraph::{
    a_star::{a_star_visit, approx_matching},
    new_path_graph::path_graph::PathGraph,
    sequences,
};

fn bench_global_alignment(c: &mut Criterion) {
    let file_path = "example/tests/L-3139.sort.gfa";
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let (sequences, _) = sequences::get_sequences(String::from("example/modified/L-3139_10.fa"));

    c.bench_function("a_star_visit", |b| {
        b.iter(|| {
            let linearized_paths = approx_matching::get_linearized_paths(&graph);
            let crumbs = approx_matching::build_heuristic(
                (&linearized_paths.0, &linearized_paths.1),
                &sequences[0],
                12,
                4,
                1,
            );
            let _ = a_star_visit::exec(
                black_box(&sequences[0]),
                black_box(&crumbs),
                black_box(&path_graph),
                false,
                4,
                12,
            );
        })
    });
}

criterion_group!(benches, bench_global_alignment);
criterion_main!(benches);
