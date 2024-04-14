use std::cmp::Ordering;

use bit_vec::BitVec;
use bstr::BString;
use handlegraph::hashgraph::HashGraph;
use longest_increasing_subsequence::lis;
use lt_fm_index::LtFmIndex;
use rayon::prelude::*;
// CHAINING
/// Get the matches chains for each path in the graph, return Vec[Vec[Match]; paths_number]
/// first get the matches for each path, then get the maximum chain of matches for each path
fn get_matches_chains(
    query_w_prefix: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) -> Vec<Vec<usize>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let mut matches: Vec<Vec<Match>> = indexes
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
    let match_chains = matches
        .iter_mut()
        .map(|path_matches| {
            // sort by path position and seed index, faster chain computation
            path_matches.sort_by_key(|m| m.path_pos);
            path_matches.sort_by_key(|m| m.seed_idx);
            get_max_chain(path_matches, seeds.len())
        })
        .collect();
    match_chains
}

/// Get the maximum chain of matches for a path
fn get_max_chain(matches: &Vec<Match>, seeds_number: usize) -> Vec<usize> {
    let matches_in_chain = lis(matches);

    let mut max_chain_seed = BitVec::from_elem(seeds_number, false);
    matches_in_chain.iter().for_each(|m| {
        max_chain_seed.set(matches[*m].seed_idx, true);
    });

    let mut sum = 0;
    let mut max_chain: Vec<usize> = max_chain_seed
        .iter()
        .rev()
        .map(|is_match| {
            sum += is_match as usize;
            sum
        })
        .collect();
    max_chain.reverse();
    max_chain
}


fn get_rec_chain(matches: &Vec<Match>) {
    let mut chains = Vec::new();
    chains.push(matches[0].clone());
    for i in 1..matches.len() {
        let mut max_len = 0;
        for j in 0..i {
            
        }
        
    }
}
/// Get the chaining heuristic for each path in the graph, return Vec[Vec[usize; query.len()]; paths_number]
/// heuristic[i][j] = x, where x is the number of seeds after the j-th that match the i-th path considering the max chain
pub fn get_chaining_sh(
    query: &BString,
    graph: &HashGraph,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) -> Vec<Vec<usize>> {
    let chains = get_matches_chains(query, chunk_size, indexes);
    let seeds_number = (query.len() - 1) / chunk_size;
    let paths_number = graph.paths.len();
    let mut heuristic = vec![vec![0; query.len()]; paths_number];
    heuristic
        .iter_mut()
        .enumerate()
        .for_each(|(path_id, path_heu)| {
            let path_chain = &chains[path_id];
            path_heu[1..(seeds_number-1) * chunk_size]
                .iter_mut()
                .enumerate()
                .for_each(|(pos, val)| {
                    // prefix is now used
                    let idx = pos / chunk_size;
                    let potential = seeds_number - idx - 1;
                    let actual = potential - path_chain[idx+1];
                    *val = actual;
                });
            path_heu[0] = path_heu[1];
        });

    heuristic
}

#[derive(Debug, Clone)]
pub struct Match {
    pub path_pos: usize,
    pub seed_idx: usize,
    pub path_id: usize,
}

impl Match {
    pub fn new() -> Self {
        Match {
            path_pos: 0,
            seed_idx: 0,
            path_id: 0,
        }
    }
    pub fn init(path_pos: usize, seed_idx: usize, path_id: usize) -> Self {
        Match { path_pos, seed_idx, path_id}
    }
}

impl PartialOrd for Match {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.path_pos < other.path_pos && self.seed_idx < other.seed_idx {
            Some(Ordering::Less)
        } else if self.path_pos > other.path_pos && self.seed_idx > other.seed_idx {
            Some(Ordering::Greater)
        } else if self.seed_idx == other.seed_idx {
            Some(Ordering::Equal)
        } else {
            None
        }
    }
}

impl Ord for Match {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
impl PartialEq for Match {
    fn eq(&self, other: &Self) -> bool {
        self.path_pos == other.path_pos && self.seed_idx == other.seed_idx
    }
}

impl Eq for Match {}
