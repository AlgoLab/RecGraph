use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rspoa::{global_mk_abpoa, graph, matrix};

fn criterion_benchmark(c: &mut Criterion) {
    let seq = "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG";
    let mut sequence = seq.chars().collect::<Vec<char>>();
    sequence.insert(0,'$');

    let graph_struct = graph::read_graph(&"prova.gfa", false);

    let score_matrix = matrix::create_score_matrix_match_mis(2, -4);

    let (b, f) = (1.0, 0.3);
    let bta = (b + f * sequence.len() as f32) as usize;

    c.bench_function("global_abpoa", |b| b.iter(|| global_mk_abpoa::exec(black_box(&sequence),("bench_sequence", 1),black_box(&graph_struct), black_box(&score_matrix), bta, "prova.gfa", false)));
}

criterion_group!(benches,  criterion_benchmark);
criterion_main!(benches);