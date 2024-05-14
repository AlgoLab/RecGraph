use bio::alignment::distance::simd::*;
use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
pub fn test(graph: &HashGraph, sequences: &Vec<BString>) {
    let mut paths_iter = graph.paths.iter().collect::<Vec<_>>();
    paths_iter.sort_by(|a, b| a.0.cmp(b.0));
    let paths = paths_iter
        .iter()
        .map(|(_, path)| {
            let path_str = path
                .nodes
                .iter()
                .flat_map(|node| graph.sequence(*node))
                .collect::<Vec<_>>();
            BString::from(path_str)
        })
        .collect::<Vec<BString>>();
    sequences.iter().for_each(|seq| {
        let best_ed = paths
            .iter()
            .enumerate()
            .map(|(idx, path)| {
                let ed = levenshtein(&seq[1..], path);
                (ed, idx)
            })
            .min()
            .unwrap();
        println!("ED\t{}\t{}", best_ed.0, best_ed.1);
    });
}
