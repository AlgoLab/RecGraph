use gfa::{
    parser::GFAParser,
    gfa::*,
};
use handlegraph::{hashgraph::HashGraph, handlegraph::HandleGraph, handle::{Handle, Direction, NodeId}};

pub fn read_graph(file_path: &str) {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let mut handles: Vec<Handle> = graph.handles_iter().collect();
    handles.sort();
    
    let mut linearization: Vec<(char, Vec<i32>)> = Vec::new();
    let mut last_index = 0;
    let mut visited_node: Vec<(NodeId, i32)> = Vec::new();
    
    for handle in &handles {
            for c in graph.sequence(*handle) {
                linearization.push((c as char, vec![]));
                last_index += 1;
            }
            
            visited_node.push((handle.id(), last_index-1));

            // devo avere il primo carattere che si riferisce al predecessore
            
    }

    for handle in handles {
        if visited_node.iter().any(|&(i,_)| i==handle.id()) {
            for predecessor in graph.handle_edges_iter(handle, Direction::Left) {
                // ho ultimo indice dell'handle, devo recuperare il primo, a quel punto aggiungo alla lista di i32 la posizione dell'ultimo carattere di 
                // predecessor (che ho visitato prima perchè grafo è ordinato)
                let pred_last_idx = get_idx(&visited_node, predecessor.id());
                let h_last_idx = get_idx(&visited_node, handle.id());
                linearization[h_last_idx as usize-graph.sequence(handle).len()+1].1.push(pred_last_idx);
            }
        }
    }
    println!("{:?}", linearization);
}

fn get_idx (visited_node: &Vec<(NodeId, i32)>, pred_id: NodeId) -> i32{
    for (node_id, idx) in visited_node.iter() {
        if *node_id == pred_id {
            return *idx
        }
    }
    panic!("unable to find predecessor")
}