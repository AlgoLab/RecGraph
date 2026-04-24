use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::hashgraph::HashGraph;
use recalign::{
    a_star::a_star_demo,
    new_path_graph::path_graph::{remove_duplicate_paths, PathGraph},
    sequences,
};

pub fn compare_est_functions(c: &mut Criterion) {
    let mut group = c.benchmark_group("est_functions");
    for gene in ["B-3136", "C-3137"].iter() {
        let seqs_path = format!("example/{}.fa", gene);
        let (sequences, _) = sequences::get_sequences(seqs_path, false);
        let graph_path = format!("example/{}.gfa", gene);
        let parser = GFAParser::new();
        let gfa: GFA<usize, ()> = parser.parse_file(graph_path).unwrap();
        let mut graph: HashGraph = HashGraph::from_gfa(&gfa);
        remove_duplicate_paths(&mut graph);
        let path_graph = PathGraph::from_hash_graph(&graph);
        let indexes = path_graph.get_indexes();
        group.bench_with_input(
            BenchmarkId::new("Chaining", gene),
            gene,
            |b, _gene: &&str| {
                b.iter(|| {
                    a_star_demo::alignment_bench(
                        black_box(&sequences),
                        black_box(&path_graph),
                        10,
                        4,
                        true,
                        2,
                        black_box(&indexes),
                        recalign::args_parser::EstimateFunction::Chaining,
                    )
                })
            },
        );
        group.bench_with_input(
            BenchmarkId::new("Seeding", gene),
            gene,
            |b, _gene: &&str| {
                b.iter(|| {
                    a_star_demo::alignment_bench(
                        black_box(&sequences),
                        black_box(&path_graph),
                        10,
                        4,
                        true,
                        2,
                        black_box(&indexes),
                        recalign::args_parser::EstimateFunction::Seeding,
                    )
                })
            },
        );
        group.bench_with_input(BenchmarkId::new("Fast", gene), gene, |b, _gene: &&str| {
            b.iter(|| {
                a_star_demo::alignment_bench(
                    black_box(&sequences),
                    black_box(&path_graph),
                    10,
                    4,
                    true,
                    2,
                    black_box(&indexes),
                    recalign::args_parser::EstimateFunction::Fast,
                )
            })
        });
    }
}

criterion_group!(benches, compare_est_functions);
criterion_main!(benches);
