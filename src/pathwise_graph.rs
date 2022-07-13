use bit_vec::BitVec;
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};
use std::collections::HashMap;
use gfa::{gfa::*, parser::GFAParser};
 
pub struct PathGraph {
    pub lnz: Vec<char>,
    pub nwp: BitVec,
    pub pred_hash: HashMap<usize, Vec<usize>>,
    pub paths_nodes: Vec<BitVec>,
    pub paths_number: usize,
}

impl PathGraph {
    pub fn new() -> PathGraph {
        PathGraph {
            lnz: vec![],
            nwp: BitVec::new(),
            pred_hash: HashMap::new(),
            paths_nodes: vec![],
            paths_number: 0,
        }
    }

    pub fn build(
        lnz: Vec<char>,
        nwp: BitVec,
        pred_hash: HashMap<usize, Vec<usize>>,
        paths_nodes: Vec<BitVec>,
        paths_number: usize,
    ) -> PathGraph {
        PathGraph {
            lnz,
            nwp,
            pred_hash,
            paths_nodes,
            paths_number,
        }
    }

    pub fn to_string(self) {
        println!("Linearization:");
        println!("{:?}", self.lnz);
        println!();

        println!("Nodes with preds:");
        println!("{:?}", self.nwp);
        println!();

        println!("Preds hash:");
        println!("{:?}", self.pred_hash);
        println!();

        println!("Paths of nodes:");
        println!("{:?}", self.paths_nodes);
        println!();

        println!("Number of paths: {}", self.paths_number);
    } 
}
pub fn read_graph_w_path(file_path: &str) -> PathGraph {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    create_path_graph(&graph)
}

pub fn create_path_graph(graph: &HashGraph) -> PathGraph {
    let mut sorted_handles = graph.handles_iter().collect::<Vec<Handle>>();
    sorted_handles.sort();

    //create graph linearization
    let mut last_index = 1;
    let mut visited_node: HashMap<NodeId, i32> = HashMap::new();
    let mut last_nodes: HashMap<NodeId, i32> = HashMap::new();
    let mut linearization: Vec<char> = vec!['$'];
    let mut handles_id_position = HashMap::new();

    for handle in sorted_handles.iter() {
        let start_position = last_index;
        for ch in graph.sequence(*handle) {
            linearization.push(ch as char);
            last_index += 1;
        }
        let end_position = last_index - 1;
        visited_node.insert(handle.id(), end_position);
        last_nodes.insert(handle.id(), end_position);
        handles_id_position.insert(handle.id(), (start_position, end_position));
    }
    linearization.push('F');

    //create nwp and pred_hash
    let mut nodes_with_pred = BitVec::from_elem(linearization.len(), false);
    let mut pred_hash: HashMap<usize, Vec<usize>> = HashMap::new();

    for handle in sorted_handles.iter() {
        if graph
            .handle_edges_iter(*handle, Direction::Left)
            .into_iter()
            .count()
            == 0
        {
            let h_last_idx = get_idx(&visited_node, handle.id());
            let handle_start_pos = h_last_idx as usize - graph.sequence(*handle).len() + 1;
            nodes_with_pred.set(handle_start_pos, true);
            update_hash(&mut pred_hash, handle_start_pos, 0);
        }
        for predecessor in graph.handle_edges_iter(*handle, Direction::Left) {
            let pred_last_idx = get_idx(&visited_node, predecessor.id());
            let h_last_idx = get_idx(&visited_node, handle.id());
            let handle_start_pos = h_last_idx as usize - graph.sequence(*handle).len() + 1;
            last_nodes.remove(&predecessor.id());
            nodes_with_pred.set(handle_start_pos, true);
            update_hash(&mut pred_hash, handle_start_pos, pred_last_idx as usize);
        }
    }

    // create nodes_paths
    let paths = &graph.paths;
    let paths_number = paths.keys().len();
    let mut paths_nodes = vec![BitVec::from_elem(paths_number, false); linearization.len()];

    paths_nodes[0] = BitVec::from_elem(paths_number, true);
    for (path_id, path) in paths.iter() {
        for handle in path.nodes.iter() {
            let (handle_start, handle_end) = handles_id_position.get(&handle.id()).unwrap();
            for idx in *handle_start..*handle_end {
                paths_nodes[idx as usize].set(*path_id as usize, true);
            }
        }
    }
    paths_nodes[linearization.len() - 1] = BitVec::from_elem(paths_number, true);

    PathGraph::build(
        linearization,
        nodes_with_pred,
        pred_hash,
        paths_nodes,
        paths_number,
    )
}

fn update_hash(hashmap: &mut HashMap<usize, Vec<usize>>, k: usize, val: usize) {
    if let Some(arr) = hashmap.get_mut(&k) {
        arr.push(val);
    } else {
        hashmap.insert(k, vec![val]);
    }
}

fn get_idx(visited_node: &HashMap<NodeId, i32>, pred_id: NodeId) -> i32 {
    *visited_node.get(&pred_id).unwrap()
}
