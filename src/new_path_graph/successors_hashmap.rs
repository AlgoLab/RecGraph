use ahash::AHashMap as HashMap;
use bit_vec::BitVec;

#[derive(Debug)]
pub struct SuccHash {
    succecessor: HashMap<usize, Vec<usize>>, // key: node, values: node_idx of succ in paths_in_nodes
    paths_in_nodes: HashMap<usize, BitVec>,
    pub paths_number: usize,
}

impl Default for SuccHash {
    fn default() -> Self {
        Self::new()
    }
}

impl SuccHash {
    pub fn new() -> SuccHash {
        SuccHash {
            succecessor: HashMap::new(),
            paths_in_nodes: HashMap::new(),
            paths_number: 0,
        }
    }

    pub fn get_node_succs_and_paths(&self, node: usize) -> Vec<(usize, &BitVec)> {
        let succs = self.succecessor.get(&node).unwrap();
        let mut result = Vec::new();
        for succ in succs {
            let paths = self.paths_in_nodes.get(succ).unwrap();
            result.push((*succ, paths));
        }
        result
    }

    pub fn get_succ_in_path(&self, node: usize, path: usize) -> Vec<usize> {
        let succs = self.succecessor.get(&node).unwrap();
        let mut res = Vec::new();
        for succ in succs {
            let paths = self.paths_in_nodes.get(succ).unwrap();
            if paths[path] {
                res.push(*succ);
            }
        }
        res
    }

    pub fn get_paths_node(&self, node: usize) -> &BitVec {
        self.paths_in_nodes.get(&node).unwrap()
    }
    pub fn set_node_successor(&mut self, curr_node: usize, succ_pos: usize) {
        if self.succecessor.get(&curr_node).is_none() {
            self.succecessor.insert(curr_node, Vec::new());
        }

        self.succecessor.get_mut(&curr_node).unwrap().push(succ_pos);
    }

    pub fn set_node_paths(&mut self, curr_node: usize, paths: BitVec) {
        self.paths_in_nodes.insert(curr_node, paths);
    }

    pub fn set_node_path(&mut self, curr_node: usize, path: usize) {
        if self.paths_in_nodes.get(&curr_node).is_none() {
            self.paths_in_nodes
                .insert(curr_node, BitVec::from_elem(self.paths_number, false));
        }
        self.paths_in_nodes
            .get_mut(&curr_node)
            .unwrap()
            .set(path, true);
    }
}
