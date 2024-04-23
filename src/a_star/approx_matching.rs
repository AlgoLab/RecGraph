use bio::pattern_matching::{myers::Myers, ukkonen};
use rayon::prelude::*;

use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};

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

pub fn build_heuristic(linearized_paths: &Vec<Vec<u8>>, query: &BString, chunk_size: usize) {
    let matches = get_matches(linearized_paths, query, chunk_size);
    let heus = matches
        .iter()
        .map(|m| get_path_max_chain(m, chunk_size, query.len() / chunk_size, query.len()))
        .collect::<Vec<_>>();
    heus.iter().for_each(|h| {
        println!("{:?}", h);
    });
}

fn get_matches(
    linearized_paths: &Vec<Vec<u8>>,
    query_w_prefix: &BString,
    chunk_size: usize,
) -> Vec<Vec<Match>> {
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
                        .map(|(pos, dist)| Match::new(*pos, *dist, seed_id))
                        .collect::<Vec<_>>()
                })
                .collect()
        })
        .collect();
    println!("Myers {:?}", start.elapsed());

    matches
}

fn get_path_max_chain(
    matches: &Vec<Match>,
    match_len: usize,
    seeds_number: usize,
    query_len: usize,
) -> Vec<usize> {
    let mut chains = vec![Link::new(); matches.len()];
    for i in 0..matches.len() {
        chains[i] = Link::init(1, i, match_len);
        for j in 0..i {
            if matches[j].seed_id < matches[i].seed_id && matches[j].pos < matches[i].pos {
                let new_len = chains[j].len + 1;
                let gap_cost = (matches[i].pos - matches[j].pos)
                    .abs_diff((matches[i].seed_id - matches[j].seed_id) * match_len);
                if gap_cost > match_len {
                    continue;
                }
                let new_score = chains[j].score + match_len - gap_cost - matches[i].dist as usize;

                if new_score > chains[i].score {
                    chains[i] = Link::init(new_len, j, new_score);
                }
            }
        }
    }

    let max_chain_ending_pos = chains
        .iter()
        .enumerate()
        .max_by_key(|x| x.1.score)
        .unwrap()
        .0;
    let mut max_chain = Vec::new();
    let mut current = max_chain_ending_pos;
    while chains[current].pred != current {
        max_chain.push(matches[current].clone());
        current = chains[current].pred;
    }
    max_chain.push(matches[current].clone());
    max_chain.reverse();
    let mut max_chain_seed = vec![None; seeds_number];
    max_chain.iter().for_each(|m| {
        max_chain_seed[m.seed_id] = Some(m);
    });
    let mut sum = 0;
    let mut rec_chain: Vec<_> = max_chain_seed
        .iter()
        .rev()
        .map(|is_match| {
            let update = if is_match.is_none() {
                2
            } else {
                is_match.as_ref().unwrap().dist as usize
            };
            sum += update;
            sum
        })
        .collect();
    rec_chain.reverse();
    let mut heu: Vec<_> = rec_chain
        .iter()
        .flat_map(|&x| std::iter::repeat(x).take(match_len))
        .collect();

    while heu.len() < query_len {
        heu.push(rec_chain[rec_chain.len() - 1]);
    }
    heu.insert(0, rec_chain[0]);
    heu
}

#[derive(Debug, Clone)]
pub struct Match {
    pub pos: usize,
    pub dist: u8,
    pub seed_id: usize,
}

impl Match {
    pub fn new(pos: usize, dist: u8, seed_id: usize) -> Self {
        Self { pos, dist, seed_id }
    }
}

#[derive(Debug, Clone)]
pub struct Link {
    pub len: usize,
    pub pred: usize,
    pub score: usize,
}

impl Link {
    pub fn new() -> Self {
        Link {
            len: 0,
            pred: 0,
            score: 0,
        }
    }
    pub fn init(len: usize, pred: usize, score: usize) -> Self {
        Link { len, pred, score }
    }
}
