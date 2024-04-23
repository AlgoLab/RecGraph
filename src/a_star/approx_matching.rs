use bio::pattern_matching::{myers::Myers, ukkonen};
use rayon::prelude::*;

use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};

use super::seeding_heurisitc::Match;

pub fn test(seq: &Vec<u8>, pattern: &Vec<u8>) {
    let start = std::time::Instant::now();
    let mut uk = ukkonen::Ukkonen::with_capacity(10, ukkonen::unit_cost);
    let occ: Vec<(usize, usize)> = uk.find_all_end(pattern, seq, 1).collect();
    println!("Ukkonen {:?}", start.elapsed());

    let seq = seq.clone();
    let start = std::time::Instant::now();
    let builder = LtFmIndexBuilder::new()
        .text_type_is_inferred()
        .set_lookup_table_kmer_size_to_default()
        .set_suffix_array_sampling_ratio_to_default();

    let fm_index = builder.build(seq).unwrap();
    let matches_pos = fm_index.locate(&pattern);
    println!("FM {:?}", start.elapsed());
}

pub fn get_linearized_paths(graph: &HashGraph) -> Vec<Vec<u8>> {
    graph
        .paths
        .iter()
        .map(|(_, path)| {
            path.nodes
                .iter()
                .map(|node| graph.sequence(node.clone()))
                .collect::<Vec<_>>()
                .concat()
        })
        .collect::<Vec<_>>()
}

pub fn build_heuristic(
    linearized_paths: &Vec<Vec<u8>>,
    query: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) {
    let matches = get_matches(linearized_paths, query, chunk_size, indexes);
}

pub fn get_matches(
    linearized_paths: &Vec<Vec<u8>>,
    query_w_prefix: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) -> Vec<Vec<(usize, usize)>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let start = std::time::Instant::now();
    let matches = linearized_paths
        .par_iter()
        .map(|path| {
            seeds
                .iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let myers = Myers::<u64>::new(*seed);
                    let occ: Vec<(usize, u8)> = myers.find_all_end(path, 1).collect();
                    occ.iter()
                        .map(|(pos, _)| (*pos, seed_id))
                        .collect::<Vec<_>>()
                })
                .collect()
        })
        .collect();
    println!("Ukkonen {:?}", start.elapsed());

    let start = std::time::Instant::now();
    let mut old_matches: Vec<Vec<Match>> = indexes
        .par_iter()
        .map(|(path_id, index)| {
            seeds
                .iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let matches_pos = index.locate(seed);
                    matches_pos
                        .iter()
                        .map(|path_pos| Match::init(*path_pos as usize, seed_id, *path_id))
                        .collect::<Vec<_>>()
                })
                .collect()
        })
        .collect();
    println!("FM {:?}", start.elapsed());
    matches
}
