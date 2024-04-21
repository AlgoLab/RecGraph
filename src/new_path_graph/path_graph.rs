use ahash::AHashMap as HashMap;
use bit_vec::BitVec;
use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::HashGraph};

use super::successors_hashmap::SuccHash;

#[derive(Debug)]
pub struct PathGraph {
    pub lnz: BString,
    pub nws: BitVec,
    pub succ_hash: SuccHash,
}

impl PathGraph {
    pub fn empty_new() -> PathGraph {
        PathGraph {
            lnz: BString::from(""),
            nws: BitVec::new(),
            succ_hash: SuccHash::new(),
        }
    }

    pub fn from_hash_graph(graph: &HashGraph) -> PathGraph {
        let mut lnz = BString::from("$");

        let mut visited_handles = BitVec::from_elem(graph.node_count() * 2 + 1, false);
        let mut handles_id_pos = HashMap::new();
        let mut succ_hash = SuccHash::new();

        succ_hash.paths_number = graph.paths.len();
        succ_hash.set_node_paths(0, BitVec::from_elem(succ_hash.paths_number, true));
        let mut nws = BitVec::new();
        nws.push(true);

        let mut last_path_pos = Vec::new();
        graph.paths.iter().for_each(|(id, path)| {
            let mut prev_handle_end = 0;
            path.nodes.iter().for_each(|node| {
                if !visited_handles[node.0 as usize] {
                    let handle_start = lnz.len();
                    lnz.append(&mut graph.sequence(*node));
                    let handle_end = lnz.len() - 1;
                    visited_handles.set(node.0 as usize, true);
                    handles_id_pos.insert(node.0 as usize, (handle_start, handle_end));
                    let mut nws_slice = BitVec::from_elem(handle_end + 1 - handle_start, false);
                    nws_slice.set(nws_slice.len() - 1, true);
                    nws.append(&mut nws_slice);
                }
                let (handle_start, handle_end) = handles_id_pos.get(&(node.0 as usize)).unwrap();
                succ_hash.set_node_path(*handle_start, *id as usize);
                succ_hash.set_node_path(*handle_end, *id as usize);
                succ_hash.set_node_successor(prev_handle_end, *handle_start);
                prev_handle_end = handles_id_pos.get(&(node.0 as usize)).unwrap().1;
            });
            last_path_pos.push(prev_handle_end);
        });

        lnz.push(b'$');

        succ_hash.set_node_paths(lnz.len() - 1, BitVec::from_elem(graph.paths.len(), true));
        for final_pos in last_path_pos {
            succ_hash.set_node_successor(final_pos, lnz.len() - 1);
        }
        PathGraph {
            lnz,
            nws,
            succ_hash,
        }
    }
}
