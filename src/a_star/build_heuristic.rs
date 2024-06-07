use std::{cmp, ops::Neg};

use ahash::AHashMap as HashMap;
use lt_fm_index::LtFmIndex;
use rayon::prelude::*;

use bstr::BString;

pub fn build_heuristic(
    indexes: &Vec<(LtFmIndex, Vec<u32>)>,
    query: &BString,
    chunk_size: usize,
    rec_cost: usize,
) -> (Vec<Vec<u32>>, Vec<HashMap<(usize, usize), (usize, usize)>>) {
    let matches = get_matches(indexes, query, chunk_size);
    let (chains, matches_pos): (Vec<_>, Vec<_>) = matches
        .par_iter()
        .enumerate()
        .map(|(path, m)| {
            get_path_max_chain(
                m,
                &indexes[path].1,
                chunk_size,
                (query.len() - 2) / chunk_size,
            )
        })
        .unzip();
    //println!("______________________");

    let heus = rec_chain_update(&chains, rec_cost, query.len(), chunk_size);

    (heus, matches_pos)
}

fn get_matches(
    indexes: &Vec<(LtFmIndex, Vec<u32>)>,
    query_w_prefix: &BString,
    chunk_size: usize,
) -> Vec<Vec<(usize, usize)>> {
    let query = &BString::from(&query_w_prefix[1..&query_w_prefix.len() - 1]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();
    let matches: Vec<_> = indexes
        .par_iter()
        .map(|(index, _)| {
            let matches: Vec<_> = seeds
                .par_iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let mut occs = index.locate(seed);
                    occs.sort();
                    occs.iter()
                        .map(|start| (*start as usize, seed_id))
                        .collect::<Vec<_>>()
                })
                .collect();
            matches
        })
        .collect();
    matches
}

fn get_path_max_chain(
    matches: &Vec<(usize, usize)>,
    lnz_pos: &Vec<u32>,
    match_len: usize,
    seeds_number: usize,
) -> (Vec<u8>, HashMap<(usize, usize), (usize, usize)>) {
    let mut chains = vec![Link::new(); matches.len()];
    for i in 0..matches.len() {
        let (pos_i, seed_i) = matches[i];
        chains[i] = Link::init(0, i, (seeds_number - 1) as i32 * match_len as i32, 1);
        for j in 0..i {
            let (pos_j, seed_j) = matches[j];
            if seed_j < seed_i && pos_j + match_len - 1 < pos_i {
                //let gap_cost = (pos_i - pos_j).abs_diff((seed_i - seed_j) * match_len);
                /*
                let gap_cost = if seed_i - seed_j == 1 {
                    (pos_i - pos_j).abs_diff((seed_i - seed_j) * match_len)
                } else {
                    let gap_update = (pos_i - pos_j).abs_diff((seed_i - seed_j) * match_len) as i32 - (seed_i as i32 - seed_j as i32 - 1);
                    cmp::max(0, gap_update) as usize
                };


                if gap_cost > match_len {
                    continue;
                }


                let gap_cost = cmp::max(
                    (pos_i - pos_j).abs_diff((seed_i - seed_j) * match_len),
                    seed_i - seed_j - 1,
                );
                */
                let gap_cost = 0;
                let new_score = chains[j].score - match_len as i32 + gap_cost as i32;

                if new_score < chains[i].score {
                    chains[i] = Link::init(gap_cost, j, new_score, chains[j].len + 1);
                }
            }
        }
    }
    if chains.is_empty() {
        return (vec![1; seeds_number], HashMap::new());
    } else {
        let max_chain_ending_pos = chains
            .iter()
            .enumerate()
            .min_by_key(|x| x.1.score)
            .unwrap_or((0, &Link::new()))
            .0;

        let mut max_chain = Vec::new();
        let mut current = max_chain_ending_pos;
        while chains[current].pred != current {
            max_chain.push((chains[current].gap, &matches[current]));
            current = chains[current].pred;
        }

        max_chain.push((chains[current].gap, &matches[current]));
        //max_chain.reverse();
        let mut max_chain_seed = vec![1; seeds_number];

        let mut match_handles = HashMap::new();
        max_chain.iter().for_each(|(score, m)| {
            max_chain_seed[m.1] = *score as u8;
            match_handles.insert(
                (m.0, m.1 * match_len + 1),
                (m.0 + match_len - 1, (m.1 + 1) * match_len),
            );
        });
        (max_chain_seed, merge_matches(&match_handles, lnz_pos))
    }
}

fn merge_matches(
    match_handles: &HashMap<(usize, usize), (usize, usize)>,
    lnz_pos: &Vec<u32>,
) -> HashMap<(usize, usize), (usize, usize)> {
    let mut matches = match_handles.iter().collect::<Vec<_>>();
    let mut merged_matches = HashMap::new();
    if matches.len() > 0 {
        matches.sort_by(|a, b| a.0 .0.cmp(&b.0 .0));
        matches.iter().for_each(|(k, v)| {
            merged_matches.insert((lnz_pos[k.0] as usize, k.1), (lnz_pos[v.0] as usize, v.1));
        });
    }
    merged_matches
}
fn rec_chain_update(
    chains: &Vec<Vec<u8>>,
    rec_cost: usize,
    query_len: usize,
    match_len: usize,
) -> Vec<Vec<u32>> {
    let mut rec_chains = vec![vec![0; chains[0].len()]; chains.len()];
    let mut best_paths = vec![0; chains[0].len()];

    // iterate
    for j in (0..rec_chains[0].len() - 1).rev() {
        let mut curr_best: Option<usize> = None;
        let best_path = best_paths[j + 1];
        for i in 0..rec_chains.len() {
            let score = rec_chains[i][j + 1] + chains[i][j + 1] as usize;
            let score_rec =
                rec_chains[best_path][j + 1] + chains[best_path][j + 1] as usize + rec_cost;

            if score < score_rec {
                rec_chains[i][j] = score;
            } else {
                rec_chains[i][j] = score_rec;
            }

            if curr_best.is_none()
                || rec_chains[i][j] + chains[i][j] as usize
                    <= rec_chains[curr_best.unwrap() as usize][j]
                        + chains[curr_best.unwrap()][j] as usize
            {
                curr_best = Some(i);
            }
        }
        best_paths[j] = curr_best.unwrap();
    }
    let mut chains_score: Vec<Vec<_>> = rec_chains
        .par_iter()
        .map(|chain| {
            chain
                .iter()
                .flat_map(|score| std::iter::repeat(*score as u32).take(match_len))
                .collect::<Vec<_>>()
        })
        .collect();

    chains_score.par_iter_mut().for_each(|chain| {
        while chain.len() < query_len - 1 {
            chain.push(chain[chain.len() - 1]);
        }
        chain.insert(0, chain[0]);
    });

    chains_score
}

#[derive(Debug, Clone)]
pub struct Link {
    pub gap: usize,
    pub pred: usize,
    pub score: i32,
    pub len: usize,
}

impl Link {
    pub fn new() -> Self {
        Link {
            gap: 0,
            pred: 0,
            score: 0,
            len: 0,
        }
    }
    pub fn init(gap: usize, pred: usize, score: i32, len: usize) -> Self {
        Link {
            gap,
            pred,
            score,
            len,
        }
    }
}

impl cmp::PartialEq for Link {
    fn eq(&self, other: &Self) -> bool {
        self == other
    }
}

impl cmp::Eq for Link {}

impl cmp::PartialOrd for Link {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl cmp::Ord for Link {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        if self.len == other.len {
            self.score.neg().cmp(&other.score.neg())
        } else {
            self.len.cmp(&other.len)
        }
    }
}
