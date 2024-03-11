use bit_vec::BitVec;

use bstr::BString;
use gfa::{gfa::*, parser::GFAParser};
use handlegraph::{
    handle::{Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};

use ahash::AHashMap as HashMap;
pub struct PathGraph {
    pub lnz: BString,
    pub nwp: BitVec,
    pub pred_hash: PredHash,
    pub paths_nodes: Vec<BitVec>,
    pub alphas: Vec<usize>,
    pub paths_number: usize,
    pub nodes_id_pos: Vec<u64>,
}

impl Default for PathGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl PathGraph {
    pub fn new() -> PathGraph {
        PathGraph {
            lnz: BString::new(vec![]),
            nwp: BitVec::new(),
            pred_hash: PredHash::new(),
            paths_nodes: vec![],
            alphas: vec![],
            paths_number: 0,
            nodes_id_pos: vec![],
        }
    }

    pub fn build(
        lnz: BString,
        nwp: BitVec,
        pred_hash: PredHash,
        paths_nodes: Vec<BitVec>,
        alphas: Vec<usize>,
        paths_number: usize,
        nodes_id_pos: Vec<u64>,
    ) -> PathGraph {
        PathGraph {
            lnz,
            nwp,
            pred_hash,
            paths_nodes,
            alphas,
            paths_number,
            nodes_id_pos,
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

#[derive(Debug)]
pub struct PredHash {
    predecessor: HashMap<usize, HashMap<usize, BitVec>>,
}

impl Default for PredHash {
    fn default() -> Self {
        Self::new()
    }
}

impl PredHash {
    pub fn new() -> PredHash {
        PredHash {
            predecessor: HashMap::new(),
        }
    }

    pub fn get_preds_and_paths(&self, curr_node: usize) -> Vec<(usize, BitVec)> {
        let preds = self.predecessor.get(&curr_node).unwrap();
        let mut preds_paths = Vec::new();
        for (pred_pos, pred_paths) in preds.iter() {
            preds_paths.push((*pred_pos, pred_paths.clone()));
        }
        preds_paths
    }

    pub fn set_preds_and_paths(
        &mut self,
        curr_node: usize,
        pred_pos: usize,
        path_id: usize,
        paths_number: usize,
    ) {
        if self.predecessor.get(&curr_node).is_none() {
            self.predecessor.insert(curr_node, HashMap::new());
        }

        if self
            .predecessor
            .get(&curr_node)
            .unwrap()
            .get(&pred_pos)
            .is_none()
        {
            self.predecessor
                .get_mut(&curr_node)
                .unwrap()
                .insert(pred_pos, BitVec::from_elem(paths_number, false));
        }
        self.predecessor
            .get_mut(&curr_node)
            .unwrap()
            .get_mut(&pred_pos)
            .unwrap()
            .set(path_id, true);
    }
}

pub fn read_graph_w_path(file_path: &str, is_reversed: bool) -> PathGraph {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    create_path_graph(&graph, is_reversed)
}

pub fn create_path_graph(graph: &HashGraph, is_reversed: bool) -> PathGraph {
    let mut sorted_handles = graph.handles_iter().collect::<Vec<Handle>>();
    sorted_handles.sort();

    if is_reversed {
        sorted_handles.reverse();
        sorted_handles = sorted_handles
            .iter()
            .map(|h| h.flip())
            .collect::<Vec<Handle>>();
    }
    //create graph linearization
    let mut last_index = 1;
    let mut visited_node: HashMap<NodeId, i32> = HashMap::new();
    let mut linearization: Vec<u8> = vec![b'$'];
    let mut handles_id_position = HashMap::new();
    let mut nodes_id_pos = Vec::new();
    nodes_id_pos.push(0);
    for handle in sorted_handles.iter() {
        let start_position = last_index;
        for ch in graph.sequence(*handle) {
            linearization.push(ch);
            nodes_id_pos.push(handle.id().into());
            last_index += 1;
        }
        let end_position = last_index - 1;
        visited_node.insert(handle.id(), end_position);
        handles_id_position.insert(handle.id(), (start_position, end_position));
    }
    linearization.push(b'F');
    nodes_id_pos.push(0);

    //create nwp, pred_hash,nodes paths and
    let mut nodes_with_pred = BitVec::from_elem(linearization.len(), false);
    let mut pred_hash_struct = PredHash::new();

    let paths_set = &graph.paths;
    let mut paths = Vec::new();
    for (_id, path) in paths_set.iter() {
        paths.push(path)
    }
    for (id, path) in paths_set.iter() {
        paths[*id as usize] = path
    }

    //let paths = &graph.paths;
    let paths_number = paths_set.keys().len();
    let mut alphas = vec![paths_number + 1; linearization.len()];
    let mut paths_nodes = vec![BitVec::from_elem(paths_number, false); linearization.len()];

    paths_nodes[0] = BitVec::from_elem(paths_number, true);
    alphas[0] = 0;
    alphas[linearization.len() - 1] = 0;
    for (path_id, path) in paths.iter().enumerate() {
        let path_nodes = if is_reversed {
            path.nodes.iter().rev().collect::<Vec<&Handle>>()
        } else {
            path.nodes.iter().collect::<Vec<&Handle>>()
        };

        for (pos, handle) in path_nodes.iter().enumerate() {
            let (handle_start, handle_end) = handles_id_position.get(&handle.id()).unwrap();
            let handle_start = *handle_start as usize;
            let handle_end = *handle_end as usize;

            for idx in handle_start..=handle_end {
                paths_nodes[idx].set(path_id, true);
                if alphas[idx] == paths_number + 1 {
                    alphas[idx] = path_id;
                }
            }

            if !nodes_with_pred[handle_start] {
                nodes_with_pred.set(handle_start, true);
            }

            if pos == 0 {
                pred_hash_struct.set_preds_and_paths(handle_start, 0, path_id, paths_number)
            } else {
                //ricava handle id pos prima, ricava suo handle end e aggiorna hash
                let pred = path_nodes[pos - 1];
                let pred_end = handles_id_position.get(&pred.id()).unwrap().1;
                pred_hash_struct.set_preds_and_paths(
                    handle_start,
                    pred_end as usize,
                    path_id,
                    paths_number,
                );

                // se ultimo nodo path aggiorna anche F
                if pos == path_nodes.iter().len() - 1 {
                    pred_hash_struct.set_preds_and_paths(
                        linearization.len() - 1,
                        handle_end,
                        path_id,
                        paths_number,
                    );
                }
            }
        }
    }
    nodes_with_pred.set(linearization.len() - 1, true);
    paths_nodes[linearization.len() - 1] = BitVec::from_elem(paths_number, true);

    PathGraph::build(
        BString::new(linearization),
        nodes_with_pred,
        pred_hash_struct,
        paths_nodes,
        alphas,
        paths_number,
        nodes_id_pos,
    )
}

pub fn create_reverse_path_graph(forward_graph: &PathGraph) -> PathGraph {
    // create reverse predecessor
    let mut nodes_with_pred_rev = BitVec::from_elem(forward_graph.lnz.len(), false);
    let mut pred_hash_struct_rev = PredHash::new();

    for (node, predecessors) in forward_graph.pred_hash.predecessor.iter() {
        for (pred, paths) in predecessors.iter() {
            if !nodes_with_pred_rev[*pred] {
                nodes_with_pred_rev.set(*pred, true);
            }
            for (path_id, path) in paths.iter().enumerate() {
                if path {
                    pred_hash_struct_rev.set_preds_and_paths(
                        *pred,
                        *node,
                        path_id,
                        forward_graph.paths_number,
                    );
                }
            }
        }
    }

    PathGraph::build(
        forward_graph.lnz.clone(),
        nodes_with_pred_rev,
        pred_hash_struct_rev,
        forward_graph.paths_nodes.clone(),
        forward_graph.alphas.clone(),
        forward_graph.paths_number,
        forward_graph.nodes_id_pos.clone(),
    )
}

pub fn nodes_displacement_matrix(graph: &PathGraph, rev_graph: &PathGraph) -> Vec<Vec<i32>> {
    let paths = &graph.paths_nodes;

    let dfe = get_distance_from_end(graph);

    let dfs = get_distance_from_start(rev_graph);

    let mut ndm = vec![vec![0; paths.len()]; paths.len()];
    for i in 0..paths.len() {
        for j in 0..paths.len() {
            if i == j {
                ndm[i][j] = 0;
            } else {
                let distance = (dfs[i] - dfs[j]).abs() + (dfe[i] - dfe[j]).abs();
                let displacement = distance as i32;
                ndm[i][j] = displacement;
            }
        }
    }

    ndm
}
fn get_distance_from_start(graph: &PathGraph) -> Vec<isize> {
    let nwp = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let lnz_len = graph.lnz.len();
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[0] = 0;
    for (p, _) in pred_hash.get_preds_and_paths(0) {
        r_values[p] = 1;
    }

    for i in 1..lnz_len - 1 {
        if r_values[i] == -1 || r_values[i] > r_values[i - 1] + 1 {
            r_values[i] = r_values[i - 1] + 1;
        }
        if nwp[i] {
            for (p, _) in pred_hash.get_preds_and_paths(i) {
                if r_values[p] == -1 || r_values[p] > r_values[i] + 1 {
                    r_values[p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values
}
fn get_distance_from_end(graph: &PathGraph) -> Vec<isize> {
    let nwp = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let lnz_len = graph.lnz.len();
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[lnz_len - 1] = 0;

    for (p, _) in pred_hash.get_preds_and_paths(lnz_len - 1) {
        r_values[p] = 1;
    }

    for i in (1..lnz_len - 1).rev() {
        if r_values[i] == -1 || r_values[i] > r_values[i + 1] + 1 {
            r_values[i] = r_values[i + 1] + 1;
        }
        if nwp[i] {
            for (p, _) in pred_hash.get_preds_and_paths(i) {
                if r_values[p] == -1 || r_values[p] > r_values[i] + 1 {
                    r_values[p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values
}
