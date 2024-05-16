use super::super::new_path_graph::path_graph::PathGraph;
use bio::alignment::distance::simd::*;
use bstr::BString;
pub fn test(graph: &PathGraph, sequences: &Vec<BString>) {
    let paths = (0..graph.succ_hash.paths_number)
        .into_iter()
        .map(|idx| graph.extract_path(idx as usize))
        .collect::<Vec<_>>();

    sequences.iter().for_each(|seq| {
        let best_ed = paths
            .iter()
            .enumerate()
            .map(|(idx, (path, _))| {
                let ed = levenshtein(&seq[1..], path);
                (ed, idx)
            })
            .min()
            .unwrap();
        println!("ED\t{}\t{}", best_ed.0, best_ed.1);
    });
}
