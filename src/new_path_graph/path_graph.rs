use ahash::AHashMap as HashMap;
use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};
use lt_fm_index::LtFmIndex;
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
        }
    }

    pub fn from_hash_graph(graph: &HashGraph) -> PathGraph {
        let mut lnz = BString::from("$");

        let min_node_id: u64 = graph.min_node_id().into();
        let max_node_id: u64 = graph.max_node_id().into();
        let mut visited_handles =
            BitVec::from_elem((max_node_id - min_node_id + 1) as usize * 2 + 1, false);
        let mut handles_id_pos: HashMap<u32, (u32, u32)> = HashMap::new();

        let mut succ_hash = SuccHash::new();

        succ_hash.paths_number = graph.paths.len() as u32;
        succ_hash.set_node_paths(0, BitVec::from_elem(succ_hash.paths_number as usize, true));
        let mut nws = BitVec::new();
        nws.push(true);
        let mut handles_ids = vec![0];
        let mut last_path_pos = vec![0; graph.paths.len()];
        let mut paths_composition = vec![Vec::new(); graph.paths.len()];
        graph.paths.iter().for_each(|(id, path)| {
            let mut prev_handle_end = 0;
            path.nodes.iter().for_each(|node| {
                let handle_id: u32 = node.0 as u32;
                if !visited_handles[cast_handle_id(node.0, min_node_id)] {
                    let handle_start = lnz.len() as u32;
                    lnz.append(&mut graph.sequence(*node));
                    handles_ids.append(&mut vec![handle_id; graph.sequence(*node).len()]);
                    let handle_end = lnz.len() as u32 - 1;
                    visited_handles.set(cast_handle_id(node.0, min_node_id), true);
                    handles_id_pos.insert(handle_id, (handle_start, handle_end));
                    let mut nws_slice =
                        BitVec::from_elem((handle_end + 1 - handle_start) as usize, false);
                    nws_slice.set(nws_slice.len() - 1, true);
                    nws.append(&mut nws_slice);
                }
                let (handle_start, handle_end) = handles_id_pos.get(&(node.0 as u32)).unwrap();
                paths_composition[*id as usize].push((*handle_start, *handle_end));
                succ_hash.set_node_path(handle_id, *id as u32);
                succ_hash.set_node_successor(prev_handle_end, *handle_start);
                prev_handle_end = handles_id_pos.get(&(node.0 as u32)).unwrap().1;
            });
            last_path_pos[*id as usize] = prev_handle_end;
        });

        lnz.push(b'$');
        let max_handle_id = graph.handles_iter().max().unwrap().0 as u32 + 1;
        handles_ids.push(max_handle_id);
        succ_hash.set_node_paths(max_handle_id, BitVec::from_elem(graph.paths.len(), true));
        for final_pos in last_path_pos.iter() {
            succ_hash.set_node_successor(*final_pos, lnz.len() as u32 - 1);
        }
        PathGraph {
            lnz,
            nws,
            succ_hash,
            handles_ids,
            ending_positions: last_path_pos,
            paths_composition,
        }
    }

    pub fn get_node_path(&self, node: u32) -> &BitVec {
        self.succ_hash
            .get_paths_node(self.handles_ids[node as usize])
    }

    fn extract_path(&self, path_id: usize) -> (BString, Vec<u32>) {
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

    pub fn get_indexes(&self) -> Vec<(LtFmIndex, Vec<u32>)> {
        (0..self.succ_hash.paths_number as usize)
            .into_par_iter()
            .map(|path_id| {
                let (path, positions) = self.extract_path(path_id);
                let builder = lt_fm_index::LtFmIndexBuilder::new()
                    .text_type_is_nucleotide_with_noise()
                    .set_suffix_array_sampling_ratio_to_default()
                    .set_lookup_table_kmer_size_to_default();
                (builder.build(path.to_vec()).unwrap(), positions)
            })
            .collect::<Vec<_>>()
    }
}

fn cast_handle_id(handle_id: u64, min: u64) -> usize {
    (handle_id - min) as usize
}
