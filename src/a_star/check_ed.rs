use super::super::new_path_graph::path_graph::PathGraph;
use bio::alignment::distance::simd::*;
use bio::alignment::pairwise::*;
use bstr::BString;
use rayon::prelude::*;

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
                let ed = levenshtein(&seq[1..&seq.len() - 1], path);
                (ed, idx)
            })
            .min()
            .unwrap();
        eprintln!("ED\t{}\t{}", best_ed.0, best_ed.1);
    });
}

pub fn semiglobal_test(graph: &PathGraph, sequences: &Vec<BString>) {
    let paths = (0..graph.succ_hash.paths_number)
        .into_iter()
        .map(|idx| graph.extract_path(idx as usize))
        .collect::<Vec<_>>();

    sequences.iter().for_each(|seq| {
        let best_ed = paths
            .par_iter()
            .enumerate()
            .map(|(idx, (path, _))| {
                let score = |a: u8, b: u8| if a == b { 0i32 } else { -1i32 };
                let scoring = Scoring::new(0, -1, &score);

                let mut aligner = banded::Aligner::with_scoring(scoring, 8, 10);
                let alignment = aligner.semiglobal(&seq[1..&seq.len() - 1], path);
                let ed = -alignment.score;
                (ed, idx)
            })
            .min()
            .unwrap();
        eprintln!("ED\t{}\t{}", best_ed.0, best_ed.1);
    });
}
