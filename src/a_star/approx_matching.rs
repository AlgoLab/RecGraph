use std::{
    cmp::Ordering,
    collections::HashMap,
    hash::{Hash, Hasher},
    time::Instant,
};

use bio::pattern_matching::myers::Myers;
use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
use rayon::prelude::*;

use crate::args_parser::ClArgs;

pub fn get_linearized_paths_and_handles(graph: &HashGraph) -> Vec<(Vec<u8>, Vec<u8>)> {
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
            let path_nodes = path
                .nodes
                .iter()
                .flat_map(|node| vec![node.0 as u8; graph.sequence(*node).len()])
                .collect::<Vec<_>>();
            (*path_id as usize, seq, path_nodes)
        })
        .collect::<Vec<_>>();
    lnz_paths.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    lnz_paths
        .iter()
        .map(|(_, path, handles)| (path.to_owned(), handles.to_owned()))
        .collect()
}

pub fn build_heuristic(
    linearized_paths: &Vec<(Vec<u8>, Vec<u8>)>,
    query: &BString,
    chunk_size: usize,
) -> Vec<Vec<usize>> {
    let matches = get_matches(linearized_paths, query, chunk_size);

    let start = Instant::now();
    let mut heus = matches
        .iter()
        .map(|m| get_path_max_chain(m, chunk_size, query.len() / chunk_size, query.len()))
        .collect::<Vec<_>>();
    println!("\tnormal: {:?}", start.elapsed());
    let filtered_match = get_matches_coord(&matches, linearized_paths, chunk_size);
    println!("{:?}", filtered_match);
    let start = Instant::now();

    let mut flat_matches: Vec<(&Match, Option<usize>)> = matches
        .iter()
        .enumerate()
        .flat_map(|(i, m)| m.iter().map(move |x| (x, Some(i))))
        .collect();
    flat_matches.sort_by_key(|x| x.0.pos);
    flat_matches.sort_by_key(|x| x.0.seed_id);
    let rec_cost = ClArgs::parse().base_rec_cost;
    rec_chain_update(
        &flat_matches,
        chunk_size,
        rec_cost as usize,
        query.len() / chunk_size,
        query.len(),
        &mut heus,
    );
    println!("\trec: {:?}", start.elapsed());
    heus
}

fn get_matches(
    linearized_paths: &Vec<(Vec<u8>, Vec<u8>)>,
    query_w_prefix: &BString,
    chunk_size: usize,
) -> Vec<Vec<Match>> {
    let query = &BString::from(&query_w_prefix[1..]);

    let seeds: Vec<_> = query.chunks_exact(chunk_size).collect::<Vec<_>>();

    let matches = linearized_paths
        .par_iter()
        .enumerate()
        .map(|(path_id, (path, _))| {
            seeds
                .iter()
                .enumerate()
                .flat_map(|(seed_id, seed)| {
                    let myers = Myers::<u64>::new(*seed);
                    let occ: Vec<(usize, u8)> = myers.find_all_end(path, 1).collect();
                    occ.iter()
                        .map(|(pos, dist)| Match::new(*pos, *dist, seed_id, path_id))
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
    let mut chain_score: Vec<_> = max_chain_seed
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
    chain_score.reverse();
    let mut heu: Vec<_> = chain_score
        .iter()
        .flat_map(|&x| std::iter::repeat(x).take(match_len))
        .collect();

    while heu.len() < query_len - 1 {
        heu.push(heu[heu.len() - 1]);
    }
    heu.insert(0, chain_score[0]);
    heu
}

fn get_matches_coord(
    matches: &Vec<Vec<Match>>,
    indexes: &Vec<(Vec<u8>, Vec<u8>)>,
    match_len: usize,
) -> Vec<(MatchCoord, BitVec)> {
    let mut matches_coord: HashMap<MatchCoord, BitVec> = HashMap::new();
    matches.iter().enumerate().for_each(|(path, path_matches)| {
        path_matches.iter().for_each(|m| {
            let coord = MatchCoord::build(m.pos, &indexes[m.path_id].1, m.seed_id, match_len);
            if matches_coord.contains_key(&coord) {
                matches_coord.get_mut(&coord).unwrap().set(path, true);
            } else {
                matches_coord.insert(coord.clone(), BitVec::from_elem(indexes.len(), false));
                matches_coord.get_mut(&coord).unwrap().set(path, true);
            }
        })
    });

    matches_coord.into_iter().collect()
}

// TODO: consider score for match coord
// Way to get distance rapidly
/*
fn new_rec_chain_update(
    matches: &Vec<(MatchCoord, BitVec)>,
    match_len: usize,
    rec_cost: usize,
    seeds_number: usize,
    query_len: usize,
    heuristics: &mut Vec<Vec<usize>>,
) {
    let mut chains = vec![Link::new(); matches.len()];

    for (i, (m_i, path_i)) in matches.iter().enumerate() {
        chains[i] = Link::init(1, i, match_len);

        for (j, (m_j, path_j)) in matches.iter().enumerate().take(i) {

            if  m_j.seed_id < m_i.seed_id && m_j < m_i{
                let new_len = chains[j].len + 1;
                let gap_cost = (m_i.pos - m_j.pos).abs_diff((m_i.seed_id - m_j.seed_id) * match_len);

                if gap_cost <= match_len {
                    let new_score = if path_i == path_j {
                        chains[j].score + match_len - gap_cost - m_i.dist as usize
                    } else {
                        chains[j].score - rec_cost + match_len - gap_cost - m_i.dist as usize
                    };

                    if new_score > chains[i].score {
                        chains[i] = Link::init(new_len, j, new_score);
                    }
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
        max_chain_seed[m.path_id] = Some((m, path));
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
                if path.unwrap_or(0) == path_id && rec_score < *val {
                    *val = rec_score;
                }
            });
        });
}
*/

fn rec_chain_update(
    matches: &Vec<(&Match, Option<usize>)>,
    match_len: usize,
    rec_cost: usize,
    seeds_number: usize,
    query_len: usize,
    heuristics: &mut Vec<Vec<usize>>,
) {
    let chains = build_chain_links(matches, match_len, rec_cost);

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
                if path.unwrap_or(0) == path_id && rec_score < *val {
                    *val = rec_score;
                }
            });
        });
}

fn build_chain_links(
    matches: &Vec<(&Match, Option<usize>)>,
    match_len: usize,
    rec_cost: usize,
) -> Vec<Link> {
    let mut chains = vec![Link::new(); matches.len()];

    for (i, (m_i, path_i)) in matches.iter().enumerate() {
        chains[i] = Link::init(1, i, match_len);

        for (j, (m_j, path_j)) in matches.iter().enumerate().take(i) {
            if m_j.seed_id < m_i.seed_id && m_j.pos < m_i.pos {
                let new_len = chains[j].len + 1;
                let gap_cost =
                    (m_i.pos - m_j.pos).abs_diff((m_i.seed_id - m_j.seed_id) * match_len);

                if gap_cost <= match_len {
                    let new_score = if path_i == path_j {
                        chains[j].score + match_len - gap_cost - m_i.dist as usize
                    } else {
                        chains[j].score - rec_cost + match_len - gap_cost - m_i.dist as usize
                    };

                    if new_score > chains[i].score {
                        chains[i] = Link::init(new_len, j, new_score);
                    }
                }
            }
        }
    }
    chains
}
#[derive(Debug, Clone)]
pub struct Match {
    pub pos: usize,
    pub dist: u8,
    pub seed_id: usize,
    pub path_id: usize,
}

impl Match {
    pub fn new(pos: usize, dist: u8, seed_id: usize, path_id: usize) -> Self {
        Self {
            pos,
            dist,
            seed_id,
            path_id,
        }
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

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct MatchNode {
    pub node_id: usize,
    pub offset: usize,
}

impl MatchNode {
    pub fn new(node_id: usize, offset: usize) -> MatchNode {
        MatchNode { node_id, offset }
    }
}

impl Ord for MatchNode {
    fn cmp(&self, other: &Self) -> Ordering {
        self.node_id.cmp(&other.node_id)
    }
}

impl PartialOrd for MatchNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct MatchCoord {
    pub node_start: MatchNode,
    pub node_end: MatchNode,
    pub seed_id: usize,
}

impl MatchCoord {
    pub fn empty_new() -> MatchCoord {
        MatchCoord {
            node_start: MatchNode::new(0, 0),
            node_end: MatchNode::new(0, 0),
            seed_id: 0,
        }
    }
    pub fn build(
        position: usize,
        handles_pos: &[u8],
        seed_id: usize,
        match_len: usize,
    ) -> MatchCoord {
        let start_handle_id = handles_pos[position] as usize;
        let mut start_offset = 0;
        let mut start = position;
        while start > 0 && handles_pos[start - 1] == start_handle_id as u8 {
            start -= 1;
            start_offset += 1;
        }

        let end_handle_id = handles_pos[position + match_len - 1] as usize;
        let mut end = position + match_len - 1;
        let mut end_offset = 0;

        while end > 0 && handles_pos[end - 1] == end_handle_id as u8 {
            end -= 1;
            end_offset += 1;
        }

        MatchCoord {
            node_start: MatchNode::new(start_handle_id, start_offset),
            node_end: MatchNode::new(end_handle_id, end_offset),
            seed_id,
        }
    }

    pub fn equal(&self, other: &MatchCoord) -> bool {
        self.node_start == other.node_start && self.node_end == other.node_end
    }
}

impl Ord for MatchCoord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.node_start.cmp(&other.node_start)
    }
}

impl PartialOrd for MatchCoord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_coord_equal() {
        let coord1: MatchCoord = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(1, 0),
            seed_id: 0,
        };
        let coord2 = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(1, 0),
            seed_id: 0,
        };
        let coord3 = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(1, 0),
            seed_id: 1,
        };
        let coord4 = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(1, 1),
            seed_id: 0,
        };
        let coord5 = MatchCoord {
            node_start: MatchNode::new(1, 1),
            node_end: MatchNode::new(1, 0),
            seed_id: 0,
        };
        let coord6 = MatchCoord {
            node_start: MatchNode::new(1, 1),
            node_end: MatchNode::new(1, 1),
            seed_id: 0,
        };
        assert_eq!(coord1, coord2);
        assert_ne!(coord1, coord3);
        assert_ne!(coord1, coord4);
        assert_ne!(coord1, coord5);
        assert_ne!(coord1, coord6);
        assert_ne!(coord4, coord5);
        assert_ne!(coord4, coord6);
        assert_ne!(coord5, coord6);
    }

    // TODO: test for ordering between match coord
    #[test]
    fn test_coord_ordering() {
        let coord1: MatchCoord = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(2, 0),
            seed_id: 0,
        };
        let coord2 = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(1, 0),
            seed_id: 1,
        };
        let coord3 = MatchCoord {
            node_start: MatchNode::new(3, 0),
            node_end: MatchNode::new(4, 0),
            seed_id: 0,
        };
        let coord4 = MatchCoord {
            node_start: MatchNode::new(1, 0),
            node_end: MatchNode::new(2, 10),
            seed_id: 0,
        };
        let coord5 = MatchCoord {
            node_start: MatchNode::new(3, 5),
            node_end: MatchNode::new(5, 6),
            seed_id: 0,
        };
        let coord6 = MatchCoord {
            node_start: MatchNode::new(2, 10),
            node_end: MatchNode::new(6, 0),
            seed_id: 1,
        };
        let coord7 = MatchCoord {
            node_start: MatchNode::new(1, 1),
            node_end: MatchNode::new(1, 0),
            seed_id: 1,
        };
        let coord8 = MatchCoord {
            node_start: MatchNode::new(1, 1),
            node_end: MatchNode::new(1, 1),
            seed_id: 1,
        };

        assert!(!(coord1 < coord2));
        assert!(!(coord2 < coord1));

        assert!(coord1 < coord3);
        assert!(coord4 < coord5);

        assert!(!(coord4 < coord5));
    }
}
