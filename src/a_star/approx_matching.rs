
use bio::pattern_matching::myers::Myers;
use rayon::prelude::*;

use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};

use crate::args_parser::ClArgs;

pub fn get_linearized_paths(graph: &HashGraph) -> Vec<Vec<u8>> {
    let mut lnz_paths = graph
        .paths
        .iter()
        .map(|(path_id, path)| {
            let seq = path
                .nodes
                .iter()
                .map(|node| graph.sequence(node.clone()))
                .collect::<Vec<_>>()
                .concat();
            (*path_id as usize, seq)
        })
        .collect::<Vec<_>>();
    lnz_paths.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    lnz_paths.iter().map(|(_, path)| path.to_owned()).collect()
}

pub fn build_heuristic(
    linearized_paths: &Vec<Vec<u8>>,
    query: &BString,
    chunk_size: usize,
) -> Vec<Vec<usize>> {
    let matches = get_matches(linearized_paths, query, chunk_size);
    let mut heus = matches
        .iter()
        .map(|m| get_path_max_chain(m, chunk_size, query.len() / chunk_size, query.len()))
        .collect::<Vec<_>>();

    let flat_matches: Vec<(&Match, usize)> = matches
        .iter()
        .enumerate()
        .flat_map(|(i, m)| m.iter().map(move |x| (x, i)))
        .collect();

    let rec_cost = ClArgs::parse().base_rec_cost;
    rec_chain_update(
        &flat_matches,
        chunk_size,
        rec_cost as usize,
        query.len() / chunk_size,
        query.len(),
        &mut heus,
    );
    heus
}

fn get_matches(
    linearized_paths: &Vec<Vec<u8>>,
    query_w_prefix: &BString,
    chunk_size: usize,
) -> Vec<Vec<Match>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

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

    while heu.len() < query_len - 1{
        heu.push(heu[heu.len()-1]);
    }
    heu.insert(0, rec_chain[0]);
    heu
}

fn rec_chain_update(
    matches: &Vec<(&Match, usize)>,
    match_len: usize,
    rec_cost: usize,
    seeds_number: usize,
    query_len: usize,
    heuristics: &mut Vec<Vec<usize>>,
) {
    let mut chains = vec![Link::new(); matches.len()];
    for i in 0..matches.len() {
        chains[i] = Link::init(1, i, match_len);
        for j in 0..i {
            if matches[j].0.seed_id < matches[i].0.seed_id && matches[j].0.pos < matches[i].0.pos {
                let new_len = chains[j].len + 1;
                let gap_cost = (matches[i].0.pos - matches[j].0.pos)
                    .abs_diff((matches[i].0.seed_id - matches[j].0.seed_id) * match_len);
                if gap_cost > match_len {
                    continue;
                }
                let new_score = if matches[j].1 == matches[i].1 {
                    chains[j].score + match_len - gap_cost - matches[i].0.dist as usize
                } else {
                    chains[j].score - rec_cost + match_len - gap_cost - matches[i].0.dist as usize
                };

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
    max_chain.iter().for_each(|(m, path)| {
        max_chain_seed[m.seed_id] = Some((m, path));
    });
    let mut sum = 0;
    let mut current_path = max_chain[max_chain.len() - 1].1;
    let mut rec_chain: Vec<_> = max_chain_seed
        .iter()
        .rev()
        .map(|m| {
            let update = if m.is_none() {
                2
            } else if *m.as_ref().unwrap().1 == current_path {
                m.as_ref().unwrap().0.dist as usize
            } else {
                current_path = *m.as_ref().unwrap().1;
                rec_cost + m.as_ref().unwrap().0.dist as usize
            };
            sum += update;
            (sum, current_path)
        })
        .collect();
    rec_chain.reverse();
    let mut heu: Vec<_> = rec_chain
        .iter()
        .flat_map(|&x| std::iter::repeat(x).take(match_len))
        .collect();

    while heu.len() < query_len - 1 {
        heu.push(rec_chain[rec_chain.len() - 1]);
    }
    heu.insert(0, rec_chain[0]);

    heuristics
        .iter_mut()
        .enumerate()
        .for_each(|(path_id, path_heu)| {
            path_heu.iter_mut().enumerate().for_each(|(pos, val)| {
                let (rec_score, path) = heu[pos];
                if path == path_id && rec_score < *val {
                    *val = rec_score;
                }
            });
        });
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
