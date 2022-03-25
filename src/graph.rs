use std::collections::HashMap;

use gfa::{gfa::*, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};

pub fn get_linearization(file_path: &str) -> Vec<(char, Vec<usize>)> {
    let (graph, sorted_handles) = read_graph(file_path);

    let mut last_index = 1;
    let mut visited_node: HashMap<NodeId, i32> = HashMap::new();
    let mut linearization: Vec<(char, Vec<usize>)> = vec![('$', vec![])];
    let mut last_nodes: HashMap<NodeId, i32> = HashMap::new();

    // concateno tutte le sequenze
    for handle in &sorted_handles {
        for c in graph.sequence(*handle) {
            linearization.push((c as char, vec![]));
            last_index += 1;
        }
        visited_node.insert(handle.id(), last_index - 1);
        last_nodes.insert(handle.id(), last_index - 1);
    }
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
                linearization[h_last_idx as usize - graph.sequence(*handle).len() + 1]
                    .1
                    .push(0);
            }
            for predecessor in graph.handle_edges_iter(*handle, Direction::Left) {
                let pred_last_idx = get_idx(&visited_node, predecessor.id());
                let h_last_idx = get_idx(&visited_node, handle.id());
                linearization[h_last_idx as usize - graph.sequence(*handle).len() + 1]
                    .1
                    .push(pred_last_idx as usize);
                last_nodes.remove(&predecessor.id());
            }
        }
    }
    linearization.insert(linearization.len(), ('F', vec![]));
    set_last_nodes(&last_nodes, &mut linearization);
    linearization
}

fn get_idx(visited_node: &HashMap<NodeId, i32>, pred_id: NodeId) -> i32 {
    *visited_node.get(&pred_id).unwrap()
}
fn set_last_nodes(last_nodes: &HashMap<NodeId, i32>, lnz: &mut [(char, Vec<usize>)]) {
    for (_, idx) in last_nodes.iter() {
        lnz[lnz.len() - 1].1.push(*idx as usize)
    }
}
fn read_graph(file_path: &str) -> (HashGraph, Vec<Handle>) {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let mut handles: Vec<Handle> = graph.handles_iter().collect();
    handles.sort();
    (graph, handles)
}
