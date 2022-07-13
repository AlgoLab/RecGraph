use bit_vec::BitVec;
use gfa::{gfa::*, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};
use std::collections::HashMap;

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

    //create nwp, pred_hash and nodes paths
    let mut nodes_with_pred = BitVec::from_elem(linearization.len(), false);
    let mut pred_hash: HashMap<usize, Vec<usize>> = HashMap::new();

    let paths = &graph.paths;
    let paths_number = paths.keys().len();
    let mut paths_nodes = vec![BitVec::from_elem(paths_number, false); linearization.len()];

    paths_nodes[0] = BitVec::from_elem(paths_number, true);
    for (path_id, path) in paths.iter() {
        for (pos, handle) in path.nodes.iter().enumerate() {
            let (handle_start, handle_end) = handles_id_position.get(&handle.id()).unwrap();
            for idx in *handle_start..=*handle_end {
                paths_nodes[idx as usize].set(*path_id as usize, true);
            }
            if !nodes_with_pred[*handle_start as usize] {
                nodes_with_pred.set(*handle_start as usize, true);
            }
            if pos == 0 {
                update_hash(&mut pred_hash, *handle_start as usize, 0);
            } else if pos == path.nodes.iter().len() - 1 {
                update_hash(
                    &mut pred_hash,
                    linearization.len() - 1,
                    *handle_end as usize,
                );
            } else {
                //ricava handle id pos prima, ricava suo handle end e aggiorna hash
                let pred = path.nodes[pos - 1];
                let pred_end = handles_id_position.get(&pred.id()).unwrap().1;
                update_hash(&mut pred_hash, *handle_start as usize, pred_end as usize)
            }
        }
    }
    nodes_with_pred.set(linearization.len() - 1, true);
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
        if !arr.contains(&val) {
            arr.push(val);
        }
    } else {
        hashmap.insert(k, vec![val]);
    }
}

#[cfg(test)]
mod tests {
    use handlegraph::{
        handle::Edge, hashgraph::HashGraph, hashgraph::Path,
        mutablehandlegraph::MutableHandleGraph, pathgraph::PathHandleGraph,
    };
    use std::collections::HashMap;

    #[test]
    fn pathwise_graph_correctly_created() {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.append_handle("A".as_bytes());
        let h2 = graph.append_handle("T".as_bytes());
        let h3 = graph.append_handle("C".as_bytes());
        let h4 = graph.append_handle("G".as_bytes());

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));

        let p1 = graph.create_path_handle(&['1' as u8], false);
        let p2 = graph.create_path_handle(&['2' as u8], false);

        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);

        let graph_struct = super::create_path_graph(&graph);

        assert_eq!(graph_struct.paths_number, 2);
        assert_eq!(graph_struct.pred_hash.get(&1).unwrap()[0], 0);
        assert_eq!(graph_struct.pred_hash.get(&5).unwrap()[0], 4);
        assert_eq!(graph_struct.lnz, ['$', 'A', 'T', 'C', 'G', 'F']);
        assert_eq!(graph_struct.nwp[2], true);

        let paths_h2 = &graph_struct.paths_nodes[2];
        assert_eq!(paths_h2[0], true);
        assert_eq!(paths_h2[1], false);

        let paths_start = &graph_struct.paths_nodes[0];
        let paths_end = &graph_struct.paths_nodes[5];

        assert!(paths_start[0]);
        assert!(paths_start[1]);
        assert!(paths_end[0]);
        assert!(paths_end[1]);
    }
    #[test]
    pub fn multiple_starts_and_ends_pathwise() {
        let mut graph: HashGraph = HashGraph::new();
        let h1 = graph.append_handle("A".as_bytes());
        let h1_bis = graph.append_handle("B".as_bytes());

        let h2 = graph.append_handle("T".as_bytes());
        let h3 = graph.append_handle("C".as_bytes());
        let h4 = graph.append_handle("G".as_bytes());
        let h4_bis = graph.append_handle("H".as_bytes());

        graph.create_edge(&Edge(h1, h2));
        graph.create_edge(&Edge(h1, h3));
        graph.create_edge(&Edge(h2, h4));
        graph.create_edge(&Edge(h3, h4));
        graph.create_edge(&Edge(h1_bis, h4_bis));

        let p1 = graph.create_path_handle(&['1' as u8], false);
        let p2 = graph.create_path_handle(&['2' as u8], false);
        let p3 = graph.create_path_handle(&['3' as u8], false);

        graph.append_step(&p1, h1);
        graph.append_step(&p1, h2);
        graph.append_step(&p1, h4);
        graph.append_step(&p2, h1);
        graph.append_step(&p2, h3);
        graph.append_step(&p2, h4);
        graph.append_step(&p3, h1_bis);
        graph.append_step(&p3, h4_bis);

        let graph_struct = super::create_path_graph(&graph);

        assert_eq!(graph_struct.paths_number, 3);
        assert_eq!(graph_struct.pred_hash.get(&1).unwrap()[0], 0);
        assert!(graph_struct.pred_hash.get(&7).unwrap().contains(&6));
        assert!(graph_struct.pred_hash.get(&7).unwrap().contains(&5));

        let paths_h2 = &graph_struct.paths_nodes[3];
        assert_eq!(paths_h2[0], true);
        assert_eq!(paths_h2[1], false);

        let paths_start = &graph_struct.paths_nodes[0];
        let paths_end = &graph_struct.paths_nodes[7];

        assert!(paths_start[0]);
        assert!(paths_start[1]);
        assert!(paths_end[0]);
        assert!(paths_end[1]);
    }
}
