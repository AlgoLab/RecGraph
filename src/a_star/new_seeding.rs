use std::{cmp::Ordering, collections::HashMap, hash::{Hash, Hasher}};

use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{
    handlegraph::HandleGraph,
    hashgraph::{HashGraph, Path},
};
use longest_increasing_subsequence::lis;
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};
use rayon::prelude::*;

/// Get the FM index for each path in the graph, return Vec[(path_id, fm_index); paths_number]
pub fn get_fm_index(graph: &HashGraph) -> Vec<(usize, LtFmIndex, Vec<u64>)> {
    let mut fm_indexes: Vec<(usize, LtFmIndex, Vec<u64>)> = graph
        .paths
        .par_iter()
        .map(|(id, path)| {
            let (path_seq, path_nodes) = linearize_path(path, graph);
            let builder = LtFmIndexBuilder::new()
                .text_type_is_inferred()
                .set_lookup_table_kmer_size_to_default()
                .set_suffix_array_sampling_ratio_to_default();
            let fm_index = builder.build(path_seq).unwrap();
            (*id as usize, fm_index, path_nodes)
        })
        .collect();
    fm_indexes.sort_by_key(|(id, _,_)| *id);
    fm_indexes
}

fn linearize_path(path: &Path, graph: &HashGraph) -> (Vec<u8>, Vec<u64>) {
    let lnz_path = path.nodes
        .iter()
        .map(|node| graph.sequence(*node), )
        .collect::<Vec<_>>()
        .concat();
    let path_nodes = path.nodes.iter().flat_map(|node| vec![node.0 as u64;  graph.sequence(*node).len()]).collect();
    (lnz_path, path_nodes)
}

// CHAINING
/// Get the matches chains for each path in the graph, return Vec[Vec[Match]; paths_number]
/// first get the matches for each path, then get the maximum chain of matches for each path
fn get_matches_chains(
    query_w_prefix: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex, Vec<u64>)>,
    rec_cost: usize,
) -> (Vec<Vec<usize>>, Vec<(usize, usize)>) {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let mut matches: Vec<Vec<Match>> = indexes
        .par_iter()
        .map(|(path_id, index,_)| {
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
    matches.iter().for_each(|x| println!("{:?}", x));
    let prova = get_matches_coord(&matches, indexes);
    prova.iter().for_each(|x| println!("{:?}", x));
    let match_chains: Vec<Vec<usize>> = matches
        .par_iter_mut()
        .map(|path_matches| {
            // sort by path position and seed index, faster chain computation
            path_matches.sort_by_key(|m| m.path_pos);
            path_matches.sort_by_key(|m| m.seed_idx);
            get_max_chain(path_matches, seeds.len())
        })
        .collect();
    
    let mut flat_matches = matches.iter().flatten().collect::<Vec<_>>();
    flat_matches.sort_by_key(|m| m.path_pos);
    flat_matches.sort_by_key(|m| m.seed_idx);
    let rec_chain = get_rec_chain(
        &flat_matches,
        chunk_size,
        rec_cost,
        seeds.len(),
        query.len(),
    );
    (match_chains, rec_chain)
}

fn get_matches_coord(
    matches: &Vec<Vec<Match>>,
    indexes: &Vec<(usize, LtFmIndex, Vec<u64>)>,
) -> HashMap<MatchCoord, BitVec> {
    let mut matches_coord: HashMap<MatchCoord, BitVec> = HashMap::new();
    matches
        .iter()
        .enumerate()
        .for_each(|(path, path_matches)| {
            path_matches
                .iter()
                .for_each(|m| {
                    let node_id = indexes[m.path_id].2[m.path_pos];
                    let coord = MatchCoord::build(node_id, m.path_pos, &indexes[m.path_id].2);
                    if matches_coord.contains_key(&coord) {
                        matches_coord.get_mut(&coord).unwrap().set(path, true);
                    } else {
                        matches_coord.insert(coord.clone(), BitVec::from_elem(indexes.len(), false));
                        matches_coord.get_mut(&coord).unwrap().set(path, true);
                    }
                })
               
        });
        
    matches_coord
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

fn get_rec_chain(
    matches: &Vec<&Match>,
    match_len: usize,
    rec_cost: usize,
    seeds_number: usize,
    query_len: usize,
) -> Vec<(usize, usize)> {
    let mut chains = vec![Link::new(); matches.len()];
    for i in 0..matches.len() {
        chains[i] = Link::init(1, i, match_len);
        for j in 0..i {
            if matches[j].seed_idx < matches[i].seed_idx && matches[j].path_pos < matches[i].path_pos
            {
                let new_len = chains[j].len + 1;
                let gap_cost = (matches[i].path_pos - matches[j].path_pos)
                    .abs_diff((matches[i].seed_idx - matches[j].seed_idx) * match_len);
                if gap_cost > match_len {
                    continue;
                }
                let new_score = if matches[j].path_id == matches[i].path_id {
                    chains[j].score + match_len - gap_cost
                } else {
                    chains[j].score - rec_cost + match_len - gap_cost
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
    max_chain.iter().for_each(|m| {
        max_chain_seed[m.seed_idx] = Some(m);
    });
    let mut sum = 0;
    let mut current_path = max_chain[max_chain.len() - 1].path_id;
    let mut rec_chain: Vec<_> = max_chain_seed
        .iter()
        .rev()
        .map(|is_match| {
            let update = if is_match.is_none() {
                1
            } else if is_match.as_ref().unwrap().path_id == current_path {
                0
            } else {
                current_path = is_match.as_ref().unwrap().path_id;
                rec_cost
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

    while heu.len() < query_len {
        heu.push(rec_chain[rec_chain.len() - 1]);
    }
    heu.insert(0, rec_chain[0]);
    heu
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
/// Get the chaining heuristic for each path in the graph, return Vec[Vec[usize; query.len()]; paths_number]
/// heuristic[i][j] = x, where x is the number of seeds after the j-th that match the i-th path considering the max chain
pub fn get_chaining_sh(
    query: &BString,
    chunk_size: usize,
    indexes: &Vec<(usize, LtFmIndex, Vec<u64>)>,
    rec_cost: usize,
) -> Vec<Vec<usize>> {
    let (chains, rec_chain) = get_matches_chains(query, chunk_size, indexes, rec_cost);
    let seeds_number = (query.len() - 1) / chunk_size;
    let mut heuristic = build_heuristic(&chains, &query, seeds_number, chunk_size);
    update_heuristic_with_rec(&mut heuristic, &rec_chain);
    heuristic
}

fn build_heuristic(
    chains: &[Vec<usize>],
    query: &BString,
    seeds_number: usize,
    chunk_size: usize,
) -> Vec<Vec<usize>> {
    (0..chains.len())
        .into_par_iter()
        .map(|path_id| {
            let path_chain = &chains[path_id];
            let mut path_heu = vec![0; query.len()];

            path_heu[1..(seeds_number - 1) * chunk_size]
                .iter_mut()
                .enumerate()
                .for_each(|(pos, val)| {
                    let idx = pos / chunk_size;
                    let potential = seeds_number - idx - 1;
                    let actual = potential - path_chain[idx + 1];
                    *val = actual;
                });

            path_heu[0] = path_heu[1];
            path_heu
        })
        .collect()
}

fn update_heuristic_with_rec(heuristic: &mut Vec<Vec<usize>>, rec_chain: &Vec<(usize, usize)>) {
    heuristic
        .iter_mut()
        .enumerate()
        .for_each(|(path_id, path_heu)| {
            path_heu.iter_mut().enumerate().for_each(|(pos, val)| {
                let (rec_score, path) = rec_chain[pos];
                if path == path_id && rec_score < *val {
                    *val = rec_score;
                }
            });
        });
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
        Match {
            path_pos,
            seed_idx,
            path_id,
        }
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

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct MatchCoord {
    pub node_id: usize,
    pub offset: usize,
}

impl MatchCoord {
    pub fn empty_new() -> MatchCoord {
        MatchCoord {
            node_id: 0,
            offset: 0,
        }
    }
    pub fn new(position: usize) -> MatchCoord {
        MatchCoord {
            node_id: 0,
            offset: position,
        }
    }
    pub fn build(handle_id: u64, position: usize, handles_pos: &[u64]) -> MatchCoord {
        let mut offset = 0;
        let mut start = position;
        while start > 0 && handles_pos[start - 1] == handle_id {
            start -= 1;
            offset += 1;
        }

        MatchCoord {
            node_id: handle_id as usize,
            offset,
        }
    }

    pub fn equal(&self, other: &MatchCoord) -> bool {
        self.node_id == other.node_id && self.offset == other.offset
    }

    pub fn included(&self, start_coor: &MatchCoord, end_coor: &MatchCoord) -> bool {
        let after_start = match self.node_id.cmp(&start_coor.node_id) {
            Ordering::Greater => true,
            Ordering::Equal => self.offset >= start_coor.offset,
            Ordering::Less => false,
        };

        let before_end = match self.node_id.cmp(&end_coor.node_id) {
            Ordering::Greater => false,
            Ordering::Equal => self.offset <= end_coor.offset,
            Ordering::Less => true,
        };

        after_start && before_end
    }
}