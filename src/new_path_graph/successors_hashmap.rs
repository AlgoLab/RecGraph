use ahash::AHashMap as HashMap;
use bit_vec::BitVec;

#[derive(Debug)]
pub struct SuccHash {
    successor: HashMap<u32, Vec<u32>>, // key: node, values: node_idx of succ in paths_in_nodes
    paths_in_nodes: HashMap<u32, BitVec>,
    pub paths_number: u32,
}

impl Default for SuccHash {
    fn default() -> Self {
        Self::new()
    }
}

impl SuccHash {
    pub fn new() -> SuccHash {
        SuccHash {
            successor: HashMap::new(),
            paths_in_nodes: HashMap::new(),
            paths_number: 0,
        }
    }

    pub fn get_node_succs(&self, node: u32) -> &Vec<u32> {
        match self.successor.get(&node) {
            Some(succs) => succs,
            None => {
                println!("Node {node}");
                panic!("Node not found")
            }
        }
    }

    pub fn get_paths_node(&self, node: u32) -> &BitVec {
        self.paths_in_nodes.get(&node).unwrap()
    }
    pub fn set_node_successor(&mut self, curr_node: u32, succ_pos: u32) {
        if self.successor.get(&curr_node).is_none() {
            self.successor.insert(curr_node, Vec::new());
        }
        if !self.successor.get(&curr_node).unwrap().contains(&succ_pos) {
            self.successor.get_mut(&curr_node).unwrap().push(succ_pos);
        }
    }

    pub fn set_node_paths(&mut self, curr_node: u32, paths: BitVec) {
        self.paths_in_nodes.insert(curr_node, paths);
    }

    pub fn set_node_path(&mut self, curr_node: u32, path: u32) {
        if self.paths_in_nodes.get(&curr_node).is_none() {
            self.paths_in_nodes.insert(
                curr_node,
                BitVec::from_elem(self.paths_number as usize, false),
            );
        }
        self.paths_in_nodes
            .get_mut(&curr_node)
            .unwrap()
            .set(path as usize, true);
    }

    pub fn get_all_paths(&self) -> &HashMap<u32, bit_vec::BitVec> {
        &self.paths_in_nodes
    }
}
