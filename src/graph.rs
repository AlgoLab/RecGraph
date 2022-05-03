use bit_vec::BitVec;
use gfa::{gfa::*, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};
use std::collections::HashMap;

pub fn read_graph(file_path: &str, amb_mode: bool) -> LnzGraph {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    create_graph_struct(&graph, amb_mode)
}

pub struct LnzGraph {
    pub lnz: Vec<char>,
    pub nwp: BitVec,
    pub pred_hash: HashMap<usize, Vec<usize>>,
}
fn create_graph_struct(graph: &HashGraph, amb_mode: bool) -> LnzGraph {
    let mut sorted_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_handles.sort();
    if amb_mode {
        sorted_handles.reverse();
        sorted_handles = sorted_handles
            .iter()
            .map(|h| h.flip())
            .collect::<Vec<Handle>>();
    }
    let mut last_index = 1;
    let mut visited_node: HashMap<NodeId, i32> = HashMap::new();
    let mut last_nodes: HashMap<NodeId, i32> = HashMap::new();
    let mut handles_id_position = HashMap::new();
    let mut linearization: Vec<char> = vec!['$'];
    // concateno tutte le sequenze
    for handle in &sorted_handles {
        let start_pos = last_index;
        for c in graph.sequence(*handle) {
            linearization.push(c as char);
            last_index += 1;
        }
        let last_pos = last_index - 1; // last position included in position of current handle in lnz
        handles_id_position.insert(handle.id(), (start_pos, last_pos));
        visited_node.insert(handle.id(), last_index - 1);
        last_nodes.insert(handle.id(), last_index - 1);
    }
    let mut nodes_with_predecessor = BitVec::from_elem(linearization.len() + 1, false);
    let mut predecessor_hash: HashMap<usize, Vec<usize>> = HashMap::new();

    // per ogni sequenza guardo i predecessori e aggiorno valore corrispondente
    for handle in &sorted_handles {
        if visited_node.get(&handle.id()).is_some() {
            if graph
                .handle_edges_iter(*handle, Direction::Left)
                .into_iter()
                .count()
                == 0
            {
                let h_last_idx = get_idx(&visited_node, handle.id());
                let handle_start_pos = h_last_idx as usize - graph.sequence(*handle).len() + 1;
                nodes_with_predecessor.set(handle_start_pos, true);
                update_hash(&mut predecessor_hash, handle_start_pos, 0);
            }
            for predecessor in graph.handle_edges_iter(*handle, Direction::Left) {
                let pred_last_idx = get_idx(&visited_node, predecessor.id());
                let h_last_idx = get_idx(&visited_node, handle.id());
                let handle_start_pos = h_last_idx as usize - graph.sequence(*handle).len() + 1;
                last_nodes.remove(&predecessor.id());
                nodes_with_predecessor.set(handle_start_pos, true);
                update_hash(
                    &mut predecessor_hash,
                    handle_start_pos,
                    pred_last_idx as usize,
                );
            }
        }
    }

    set_last_node(
        &mut linearization,
        &mut nodes_with_predecessor,
        &mut predecessor_hash,
        &last_nodes,
    );

    LnzGraph {
        lnz: linearization,
        nwp: nodes_with_predecessor,
        pred_hash: predecessor_hash,
    }
}

fn update_hash(hashmap: &mut HashMap<usize, Vec<usize>>, k: usize, val: usize) {
    if let Some(arr) = hashmap.get_mut(&k) {
        arr.push(val);
    } else {
        hashmap.insert(k, vec![val]);
    }
}

fn set_last_node(
    linearization: &mut Vec<char>,
    nodes_with_predecessor: &mut BitVec,
    predecessor_hash: &mut HashMap<usize, Vec<usize>>,
    last_nodes: &HashMap<NodeId, i32>,
) {
    linearization.insert(linearization.len(), 'F');
    nodes_with_predecessor.set(linearization.len() - 1, true);
    for (_, idx) in last_nodes.iter() {
        update_hash(predecessor_hash, linearization.len() - 1, *idx as usize);
    }
}

fn get_idx(visited_node: &HashMap<NodeId, i32>, pred_id: NodeId) -> i32 {
    *visited_node.get(&pred_id).unwrap()
}
pub fn get_sorted_handles(file_path: &str, amb_mode: bool) -> Vec<Handle> {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let mut sorted_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_handles.sort();
    if amb_mode {
        sorted_handles.reverse();
        sorted_handles = sorted_handles
            .iter()
            .map(|h| h.flip())
            .collect::<Vec<Handle>>();
    }
    sorted_handles
}

pub fn create_nodes_paths(file_path: &str) -> Vec<Vec<usize>> {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let paths = &graph.paths;
    let mut every_path = vec![];
    for path in paths.keys() {
        every_path.push(*path as usize)
    }

    let mut sorted_handles: Vec<Handle> = graph.handles_iter().collect();
    sorted_handles.sort();
    let mut paths_node: Vec<Vec<usize>> = Vec::new();
    let mut current_position = 0;
    paths_node.insert(0, every_path);
    current_position += 1;
    for handle in sorted_handles.iter() {
        let handle_length = &graph.sequence(*handle).len();
        for i in current_position..current_position + handle_length {
            paths_node.insert(i, vec![]);
            for (path_id, path) in paths.iter() {
                if path.nodes.contains(handle) {
                    paths_node[i].push(*path_id as usize);
                }
            }
        }
        current_position += handle_length;
    }
    let mut every_path = vec![];
    for path in paths.keys() {
        every_path.push(*path as usize)
    }
    paths_node.insert(current_position, every_path);
    for p in paths_node.iter_mut() {
        p.sort_unstable()
    }
    paths_node
}
#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use handlegraph::{handle::Edge, hashgraph::HashGraph, mutablehandlegraph::MutableHandleGraph};

    #[test]
    fn graph_struct_correctly_created() {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.append_handle("A".as_bytes());
        let h2 = graph.append_handle("T".as_bytes());
        let h3 = graph.append_handle("C".as_bytes());
        let h4 = graph.append_handle("G".as_bytes());

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));

        let graph_struct = super::create_graph_struct(&graph, false);
        assert!(graph_struct.nwp[1]);
        assert!(graph_struct.nwp[5]);
        assert_eq!(graph_struct.pred_hash.get(&1).unwrap()[0], 0);
        assert_eq!(graph_struct.pred_hash.get(&5).unwrap()[0], 4);
        assert_eq!(graph_struct.lnz, ['$', 'A', 'T', 'C', 'G', 'F']);
    }
    #[test]
    fn rev_graph_struct_correctly_created() {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.append_handle("A".as_bytes());
        let h2 = graph.append_handle("T".as_bytes());
        let h3 = graph.append_handle("C".as_bytes());
        let h4 = graph.append_handle("G".as_bytes());

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h2, h3));
        graph.create_edge(&Edge(h3, h4));

        let graph_struct = super::create_graph_struct(&graph, true);
        assert!(graph_struct.nwp[1]);
        assert!(graph_struct.nwp[5]);
        assert_eq!(graph_struct.pred_hash.get(&1).unwrap()[0], 0);
        assert_eq!(graph_struct.pred_hash.get(&5).unwrap()[0], 4);
        assert_eq!(graph_struct.lnz, ['$', 'C', 'G', 'A', 'T', 'F']);
    }
    #[test]
    fn handle_id_from_lnz_pos_and_sorted_handles() {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.append_handle("A".as_bytes());
        let h2 = graph.append_handle("TA".as_bytes());
        let h3 = graph.append_handle("CGG".as_bytes());
        let h4 = graph.append_handle("G".as_bytes());
        let h5 = graph.append_handle("TCCCC".as_bytes());

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h3, h5));

        let gs = super::create_graph_struct(&graph, false);
        let mut curr_hand = -1;
        let mut handle_pos = HashMap::new();
        for i in 1..gs.lnz.len() - 1 {
            if gs.nwp[i] {
                curr_hand += 1;
            }
            handle_pos.insert(i, curr_hand);
        }
        assert_eq!(handle_pos.get(&1).unwrap(), &0);
        assert_eq!(handle_pos.get(&2).unwrap(), &1);
        assert_eq!(handle_pos.get(&4).unwrap(), &2);
        assert_eq!(handle_pos.get(&6).unwrap(), &2);
        assert_eq!(handle_pos.get(&7).unwrap(), &3);
        assert_eq!(handle_pos.get(&12).unwrap(), &4);
    }
}
