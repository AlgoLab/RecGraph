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
                (query.len() - 1) / chunk_size,
            )
        })
        .unzip();

    let heus = rec_chain_update(&chains, rec_cost, query.len(), chunk_size);
    (heus, matches_pos)
}

fn get_matches(
    indexes: &Vec<(LtFmIndex, Vec<u32>)>,
    query_w_prefix: &BString,
    chunk_size: usize,
) -> Vec<Vec<(usize, usize)>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();
    let matches: Vec<_> = indexes
        .par_iter()
        .map(|(index, _)| {
            let matches: Vec<_> = seeds
                .par_iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let occs = index.locate(seed);
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
        chains[i] = Link::init(0, i, match_len);
        for j in 0..i {
            let (pos_j, seed_j) = matches[j];
            if seed_j < seed_i && pos_j + match_len <= pos_i {
                let gap_cost = (pos_i - pos_j).abs_diff((seed_i - seed_j) * match_len);
                if gap_cost > match_len {
                    continue;
                }
                let new_score = chains[j].score + match_len - gap_cost;

                if new_score > chains[i].score {
                    chains[i] = Link::init(gap_cost, j, new_score);
                }
            }
        }
    }

    let max_chain_ending_pos = chains
        .iter()
        .enumerate()
        .max_by_key(|x| x.1.score)
        .unwrap_or((0, &Link::new()))
        .0;
    if max_chain_ending_pos == 0 {
        (vec![1; seeds_number], HashMap::new())
    } else {
        let mut max_chain = Vec::new();
        let mut current = max_chain_ending_pos;
        while chains[current].pred != current {
            max_chain.push((chains[current].val, matches[current].clone()));
            current = chains[current].pred;
        }
        max_chain.push((chains[current].val, matches[current].clone()));
        max_chain.reverse();
        let mut max_chain_seed = vec![1; seeds_number];
        let mut match_handles = HashMap::new();
        max_chain.iter().for_each(|(score, m)| {
            max_chain_seed[m.1] = *score as u8;

            match_handles.insert(
                (lnz_pos[m.0] as usize, m.1 * match_len + 1),
                (lnz_pos[m.0 + match_len - 1] as usize, (m.1 + 1) * match_len),
            );
        });

        (max_chain_seed, match_handles)
    }
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
        let mut curr_best = None;
        for i in 0..rec_chains.len() {
            let best_path = best_paths[j + 1];
            let score = rec_chains[i][j + 1] + chains[i][j + 1] as usize;
            let score_rec = rec_chains[best_path][j + 1] + chains[i][j + 1] as usize + rec_cost;
            if score < score_rec {
                rec_chains[i][j] = score;
            } else {
                rec_chains[i][j] = score_rec;
            }

            if curr_best.is_none() || rec_chains[i][j] <= rec_chains[curr_best.unwrap() as usize][j]
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
    pub val: usize,
    pub pred: usize,
    pub score: usize,
}

impl Link {
    pub fn new() -> Self {
        Link {
            val: 0,
            pred: 0,
            score: 0,
        }
    }
    pub fn init(val: usize, pred: usize, score: usize) -> Self {
        Link { val, pred, score }
    }
}
