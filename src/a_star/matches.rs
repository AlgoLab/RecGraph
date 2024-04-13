use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{
    handlegraph::HandleGraph,
    hashgraph::{HashGraph, Path},
};
use longest_increasing_subsequence::lis;
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};
use rayon::prelude::*;
use std::cmp::Ordering;

/// Get the FM index for each path in the graph, return Vec[(path_id, fm_index); paths_number]
pub fn get_fm_index(graph: &HashGraph) -> Vec<(usize, LtFmIndex)> {
    let mut fm_indexes: Vec<(usize, LtFmIndex)> = graph
        .paths
        .par_iter()
        .map(|(id, path)| {
            let path_seq = linearize_path(path, graph);
            let builder = LtFmIndexBuilder::new()
                .text_type_is_inferred()
                .set_lookup_table_kmer_size_to_default()
                .set_suffix_array_sampling_ratio_to_default();
            let fm_index = builder.build(path_seq).unwrap();
            (*id as usize, fm_index)
        })
        .collect();
    fm_indexes.sort_by_key(|(id, _)| *id);
    fm_indexes
}

fn linearize_path(path: &Path, graph: &HashGraph) -> Vec<u8> {
    path.nodes
        .iter()
        .map(|node| graph.sequence(*node))
        .collect::<Vec<_>>()
        .concat()
}

/// Get the matches for each path in the graph, return Vec[Vec[bool; seeds_number]; paths_number]
/// matches[i][j] = true if the j-th seed matches the i-th path in some position
fn get_matches(
    query: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) -> Vec<Vec<bool>> {
    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let mut matches = vec![vec![false; seeds.len()]; indexes.len()];

    indexes.iter().for_each(|(path_id, index)| {
        seeds.iter().enumerate().for_each(|(seed_id, seed)| {
            let is_match = index.count(seed) > 0;
            matches[*path_id][seed_id] = is_match;
        });
    });
    matches
}

/// Get the base heuristic for each path in the graph, return Vec[Vec[usize; query.len()]; paths_number]
/// heuristic[i][j] = x, where x is the number of seeds after the j-th that match the i-th path
pub fn get_base_sh(
    query: &BString,
    graph: &HashGraph,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex)>,
) -> Vec<Vec<usize>> {
    let matches = get_matches(query, chunk_size, indexes);
    let seeds_number = query.len() / chunk_size;
    let paths_number = graph.paths.len();
    let mut heuristic = vec![vec![0; query.len()]; paths_number];
    heuristic
        .iter_mut()
        .enumerate()
        .for_each(|(path_id, path_heu)| {
            let path_matches = &matches[path_id];
            path_heu[..seeds_number * chunk_size]
                .iter_mut()
                .enumerate()
                .for_each(|(pos, val)| {
                    let idx = pos / chunk_size;
                    let potential = seeds_number - idx;
                    let actual = potential - path_matches[idx + 1..].iter().filter(|&&x| x).count();
                    *val = actual;
                })
        });

    heuristic
}

#[derive(Debug)]
pub struct GraphMatches {
    pub matches: Vec<PathMatches>,
    pub seeds: Vec<BString>,
    pub paths: Vec<usize>,
}

#[derive(Debug)]
pub struct PathMatches {
    pub path_id: usize,
    pub matches: Vec<bool>,
}

impl PathMatches {
    pub fn new() -> Self {
        PathMatches {
            path_id: 0,
            matches: Vec::new(),
        }
    }
    pub fn init(path_id: usize, matches: Vec<bool>) -> Self {
        PathMatches { path_id, matches }
    }

    pub fn set_position(&mut self, pos: usize) {
        self.matches[pos] = true;
    }

    pub fn get_position(&self, pos: usize) -> bool {
        self.matches[pos]
    }
}
#[derive(Debug, Clone)]
pub struct Match {
    pub path_pos: usize,
    pub seed_idx: usize,
}

impl Match {
    pub fn new() -> Self {
        Match {
            path_pos: 0,
            seed_idx: 0,
        }
    }
    pub fn init(path_pos: usize, seed_idx: usize) -> Self {
        Match { path_pos, seed_idx }
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
