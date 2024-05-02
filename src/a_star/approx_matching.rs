use bio::pattern_matching::myers::Myers;
use rayon::prelude::*;

use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};

pub fn get_linearized_paths(graph: &HashGraph) -> Vec<Vec<u8>> {
    let mut lnz_paths = graph
        .paths
        .par_iter()
        .map(|(path_id, path)| {
            let seq = path
                .nodes
                .iter()
                .map(|node| graph.sequence(node.clone()))
                .collect::<Vec<_>>()
                .concat();
            (*path_id, seq)
        })
        .collect::<Vec<_>>();
    lnz_paths.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    lnz_paths.iter().map(|(_, path)| path.to_owned()).collect()
}

pub fn build_heuristic(
    linearized_paths: &Vec<Vec<u8>>,
    query: &BString,
    chunk_size: usize,
    rec_cost: usize,
    max_err: u8,
) -> Vec<Vec<u32>> {
    let matches = get_matches(linearized_paths, query, chunk_size, max_err);
    let chains = matches
        .par_iter()
        .map(|m| get_path_max_chain(m, chunk_size, query.len() / chunk_size, max_err))
        .collect::<Vec<_>>();

    let heus = rec_chain_update(&chains, rec_cost, query.len(), chunk_size);
    heus
}

fn get_matches(
    linearized_paths: &Vec<Vec<u8>>,
    query_w_prefix: &BString,
    chunk_size: usize,
    max_err: u8,
) -> Vec<Vec<Match>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let matches = linearized_paths
        .par_iter()
        .map(|path| {
            seeds
                .par_iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let mut myers = Myers::<u64>::new(*seed);
                    let occ_iter = myers.find_all(path, max_err);
                    occ_iter
                        .map(|(start, _, dist)| Match::new(start, dist, seed_id))
                        .collect::<Vec<_>>()
                })
                .collect()
        })
        .collect();

    matches
}

fn get_path_max_chain(
    matches: &Vec<Match>,
    match_len: usize,
    seeds_number: usize,
    max_err: u8,
) -> Vec<u8> {
    let mut chains = vec![Link::new(); matches.len()];
    for i in 0..matches.len() {
        chains[i] = Link::init(1, i, match_len - matches[i].dist as usize);
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
    let mut max_chain_seed = vec![max_err + 1; seeds_number];
    max_chain.iter().for_each(|m| {
        max_chain_seed[m.seed_id] = m.dist;
    });
    max_chain_seed
}

fn rec_chain_update(
    chains: &Vec<Vec<u8>>,
    rec_cost: usize,
    query_len: usize,
    match_len: usize,
) -> Vec<Vec<u32>> {
    let mut rec_chains = vec![vec![Link::new(); chains[0].len()]; chains.len()];
    let mut best_paths = vec![0; chains[0].len()];
    for j in (0..rec_chains[0].len() - 1).rev() {
        let mut curr_best = None;
        for i in 0..rec_chains.len() {
            let best_path = best_paths[j + 1];
            let score = rec_chains[i][j + 1].score + chains[i][j] as usize;
            let score_rec = rec_chains[best_path][j + 1].score + chains[i][j] as usize + rec_cost;
            if score < score_rec {
                rec_chains[i][j] = Link::init(rec_chains[i][j + 1].len + 1, i, score);
            } else {
                rec_chains[i][j] =
                    Link::init(rec_chains[best_path][j + 1].len + 1, best_path, score_rec);
            }

            if curr_best.is_none()
                || rec_chains[i][j].score <= rec_chains[curr_best.unwrap() as usize][j].score
            {
                curr_best = Some(i);
            }
        }
        best_paths[j] = curr_best.unwrap();
    }

    let mut chains_score: Vec<Vec<_>> = rec_chains
        .iter()
        .map(|chain| {
            chain[1..]
                .iter()
                .flat_map(|link| std::iter::repeat(link.score as u32).take(match_len))
                .collect::<Vec<_>>()
        })
        .collect();

    chains_score.iter_mut().for_each(|chain| {
        while chain.len() < query_len - 1 {
            chain.push(chain[chain.len() - 1]);
        }
        chain.insert(0, chain[0]);
    });

    chains_score
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
