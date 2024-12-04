 use ahash::AHashMap as HashMap;
use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
use lt_fm_index::{blocks::Block3, LtFmIndex};
use rayon::prelude::*;

use super::successors_hashmap::SuccHash;

#[derive(Debug)]
pub struct PathGraph {
    pub lnz: BString,
    pub nws: BitVec,
    pub succ_hash: SuccHash,
    pub handles_ids: Vec<u32>,
    pub ending_positions: Vec<u32>,
    paths_composition: Vec<Vec<(u32, u32)>>,
    pub original_handles: HashMap<u32, (u64, char)>,

    pub common_nodes: Vec<Vec<bool>>,
}

impl PathGraph {
    pub fn empty_new() -> PathGraph {
        PathGraph {
            lnz: BString::from(""),
            nws: BitVec::new(),
            succ_hash: SuccHash::new(),
            handles_ids: Vec::new(),
            ending_positions: Vec::new(),
            paths_composition: Vec::new(),
            original_handles: HashMap::new(),
            common_nodes: Vec::new(),
        }
    }

    pub fn from_hash_graph(graph: &HashGraph) -> PathGraph {
        let mut lnz = BString::from("$");

        let min_node_id: u64 = graph.min_node_id().into();
        let max_node_id: u64 = graph.max_node_id().into();
        let mut visited_handles =
            BitVec::from_elem((max_node_id - min_node_id + 1) as usize * 2 + 1, false);
        let mut last_idx = max_node_id * 2 + 1;
        let mut handles_id_pos: HashMap<u32, (u32, u32)> = HashMap::new();

        let mut succ_hash = SuccHash::new();

        succ_hash.paths_number = graph.paths.len() as u32;
        succ_hash.set_node_paths(0, BitVec::from_elem(succ_hash.paths_number as usize, true));
        let mut nws = BitVec::new();
        nws.push(true);
        let mut handles_ids = vec![0];
        let mut last_path_pos = vec![0; graph.paths.len()];
        let mut paths_composition = vec![Vec::new(); graph.paths.len()];
        let mut path_iterator = graph.paths.iter().collect::<Vec<_>>();
        path_iterator.sort_by(|a, b| a.0.cmp(b.0));
        let mut dup_handles: HashMap<(u64, i32), u32> = HashMap::new();
        let mut original_handles: HashMap<u32, (u64, char)> = HashMap::new();

        let mut handles_in_path: Vec<BitVec> =
            vec![
                BitVec::from_elem((max_node_id - min_node_id + 1) as usize * 2 + 1, false);
                graph.paths.len()
            ];
        //let mut handles_pos_in_path: Vec<HashMap<_,_>> = vec![HashMap::new(); graph.paths.len()];
        path_iterator.iter().for_each(|(id, path)| {
            let mut current_len = 0;
            //handles_pos_in_path[**id as usize].insert(0, 0);
            let mut path_handles: HashMap<u64, i32> = HashMap::new();
            let mut prev_handle_end = 0;
            path.nodes.iter().for_each(|node| {
                let handle_id = if let Some(iter) = path_handles.get(&node.0) {
                    if let Some(handle_id) = dup_handles.get(&(node.0, *iter)) {
                        *handle_id
                    } else {
                        last_idx += 1;
                        dup_handles.insert((node.0, *iter), last_idx as u32);
                        visited_handles.push(false);
                        handles_in_path.iter_mut().for_each(|bitvec| {
                            bitvec.push(false);
                        });
                        last_idx as u32
                    }
                } else {
                    node.0 as u32
                };
                original_handles.insert(
                    handle_id,
                    (
                        node.id().into(),
                        match node.is_reverse() {
                            true => '-',
                            false => '+',
                        },
                    ),
                );
                if !visited_handles[cast_handle_id(handle_id, min_node_id)] {
                    let handle_start = lnz.len() as u32;
                    lnz.append(&mut graph.sequence(*node));
                    handles_ids.append(&mut vec![handle_id; graph.sequence(*node).len()]);
                    let handle_end = lnz.len() as u32 - 1;
                    visited_handles.set(cast_handle_id(handle_id, min_node_id), true);
                    handles_id_pos.insert(handle_id, (handle_start, handle_end));
                    let mut nws_slice = BitVec::from_elem(graph.sequence(*node).len(), false);
                    nws_slice.set(nws_slice.len() - 1, true);
                    nws.append(&mut nws_slice);
                }
                if let Some(iter) = path_handles.get_mut(&node.0) {
                    *iter += 1;
                } else {
                    path_handles.insert(node.0, 1);
                }
                let (handle_start, handle_end) = handles_id_pos.get(&(handle_id)).unwrap();
                paths_composition[**id as usize].push((*handle_start, *handle_end));
                succ_hash.set_node_path(handle_id, **id as u32);
                succ_hash.set_node_successor(prev_handle_end, **id as u32, *handle_start);
                prev_handle_end = handles_id_pos.get(&(handle_id)).unwrap().1;
                handles_in_path[**id as usize].set(cast_handle_id(handle_id, min_node_id), true);
                //handles_pos_in_path[**id as usize].insert(handle_id, current_len);
                current_len += graph.sequence(*node).len();
            });
            last_path_pos[**id as usize] = prev_handle_end;
        });

        lnz.push(b'$');
        let max_handle_id = last_idx as u32 + 1;
        handles_ids.push(max_handle_id);
        nws.push(false);
        succ_hash.set_node_paths(max_handle_id, BitVec::from_elem(graph.paths.len(), true));
        for (path, final_pos) in last_path_pos.iter().enumerate() {
            succ_hash.set_node_successor(*final_pos, path as u32, lnz.len() as u32 - 1);
        }
        let common_nodes = find_rightest_common_node(&handles_in_path);
        //let common_nodes = find_last_common_node(&handles_in_path, &handles_pos_in_path, min_node_id);

        PathGraph {
            lnz,
            nws,
            succ_hash,
            handles_ids,
            ending_positions: last_path_pos,
            paths_composition,
            original_handles,
            common_nodes,
        }
    }

    pub fn get_node_path(&self, node: u32) -> &BitVec {
        self.succ_hash
            .get_paths_node(self.handles_ids[node as usize])
    }

    pub fn extract_path(&self, path_id: usize) -> (BString, Vec<u32>) {
        let tmp: (Vec<_>, Vec<_>) = self.paths_composition[path_id]
            .iter()
            .map(|(start, end)| {
                (
                    &self.lnz[*start as usize..=*end as usize],
                    (*start..=*end).collect::<Vec<_>>(),
                )
            })
            .unzip();
        (BString::from(tmp.0.concat()), tmp.1.concat())
    }

    pub fn get_indexes(&self) -> Vec<(LtFmIndex<u32, Block3<u128>>, Vec<u32>)> {
        let characters_by_index: &[&[u8]] = &[
            b"Aa", b"Cc", b"Gg", b"Tt"
        ];
        let mut indexes = (0..self.succ_hash.paths_number as usize)
            .into_par_iter()
            .enumerate()
            .map(|(idx, path_id)| {
                let (path, positions) = self.extract_path(path_id);
                let index = LtFmIndex::<u32, Block3<u128>>::build(
                    path.to_vec(),
                    characters_by_index,
                    2,
                    4
                ).unwrap();
                
                (idx, index, positions)
            })
            .collect::<Vec<_>>();
        indexes.sort_by(|a, b| a.0.cmp(&b.0));
        indexes
            .into_iter()
            .map(|(_, index, positions)| (index, positions))
            .collect()
    }

    pub fn get_graph_size(&self) -> usize {
        (0..self.succ_hash.paths_number as usize)
            .map(|path_id| {
                let (path, _) = self.extract_path(path_id);
                path.len()
            })
            .sum()
    }
}

fn cast_handle_id(handle_id: u32, min: u64) -> usize {
    (handle_id as u64 - min * 2) as usize
}

/* 
fn reverse_cast_handle_id(handle_id: usize, min: u64) -> u32 {
    (handle_id as u64 + min * 2) as u32
}
*/
pub fn remove_duplicate_paths(graph: &mut HashGraph) -> HashMap<u8, u8> {
    let mut paths_to_remove = Vec::new();
    graph
        .paths
        .iter()
        .enumerate()
        .for_each(|(idx, (path_id, path))| {
            graph.paths.iter().skip(idx).for_each(|(path_id2, path2)| {
                if path_id != path_id2 && path.nodes == path2.nodes {
                    paths_to_remove.push(*path_id2);
                }
            });
        });

    graph
        .paths
        .retain(|path_id, _| !paths_to_remove.contains(path_id));

    // change id back to 0..n
    let mut old_graph = graph.paths.drain().collect::<Vec<_>>();
    old_graph.sort_by(|a, b| a.0.cmp(&b.0));
    let mut original_path_ids = HashMap::new();
    old_graph
        .drain(0..)
        .enumerate()
        .for_each(|(idx, (original_id, path))| {
            graph.paths.insert(idx as i64, path);
            original_path_ids.insert(idx as u8, original_id as u8);
        });
    original_path_ids
}

fn find_rightest_common_node(handles_in_path: &Vec<BitVec>) -> Vec<Vec<bool>> {
    let mut common_nodes: Vec<_> = vec![vec![false; handles_in_path.len()]; handles_in_path.len()];
    handles_in_path.iter().enumerate().for_each(|(i, bitvec)| {
        handles_in_path.iter().enumerate().for_each(|(j, bitvec2)| {
            let mut res = bitvec.clone();
            res.and(bitvec2);
            common_nodes[i][j] = res.any();
        });
    });
    common_nodes
}
/* 
fn find_last_common_node(
    handles_in_path: &Vec<BitVec>,
    handles_pos_in_path: &Vec<HashMap<u32, usize>>,
    min: u64,
) -> Vec<Vec<u32>> {
    let mut common_nodes: Vec<_> = vec![vec![0; handles_in_path.len()]; handles_in_path.len()];
    handles_in_path.iter().enumerate().for_each(|(i, bitvec)| {
        handles_in_path.iter().enumerate().for_each(|(j, bitvec2)| {
            let mut res = bitvec.clone();
            res.and(bitvec2);
            let handle = res
                .iter()
                .enumerate()
                .rev()
                .find(|(_, x)| *x)
                .unwrap_or((0, true))
                .0;
            if handle != 0 {
                let original_handle = reverse_cast_handle_id(handle, min);
                common_nodes[i][j] = *handles_pos_in_path[j].get(&original_handle).unwrap() as u32;
            } else {
                common_nodes[i][j] = 0;
            }
        });
    });
    common_nodes
}
*/