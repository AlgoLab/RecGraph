use gfa::{gfa::*, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Handle, NodeId},
    handlegraph::HandleGraph,
    hashgraph::HashGraph,
};

pub fn get_linearization(file_path: &str) -> Vec<(char, Vec<usize>)> {
    let (graph, sorted_handles) = read_graph(file_path);

    let mut last_index = 1;
    let mut visited_node: Vec<(NodeId, i32)> = Vec::new();
    let mut linearization: Vec<(char, Vec<usize>)> = Vec::new();
    linearization.push(('$', vec![]));

    // concateno tutte le sequenze
    for handle in &sorted_handles {
        for c in graph.sequence(*handle) {
            linearization.push((c as char, vec![]));
            last_index += 1;
        }

        visited_node.push((handle.id(), last_index - 1));
    }
    // per ogni sequenza guardo i predecessori e aggiorno valore corrispondente
    for handle in &sorted_handles {
        if visited_node.iter().any(|&(i, _)| i == handle.id()) {
            if graph
                .handle_edges_iter(*handle, Direction::Left)
                .into_iter()
                .count()
                == 0
            {
                let h_last_idx = get_idx(&visited_node, handle.id());
                linearization[h_last_idx as usize - graph.sequence(*handle).len() + 1]
                    .1
                    .push(0 as usize);
            }
            for predecessor in graph.handle_edges_iter(*handle, Direction::Left) {
                let pred_last_idx = get_idx(&visited_node, predecessor.id());
                let h_last_idx = get_idx(&visited_node, handle.id());
                linearization[h_last_idx as usize - graph.sequence(*handle).len() + 1]
                    .1
                    .push(pred_last_idx as usize);
            }
        }
    }
    linearization
}

fn get_idx(visited_node: &Vec<(NodeId, i32)>, pred_id: NodeId) -> i32 {
    for (node_id, idx) in visited_node.iter() {
        if *node_id == pred_id {
            return *idx;
        }
    }
    panic!("unable to find predecessor")
}

fn read_graph(file_path: &str) -> (HashGraph, Vec<Handle>) {
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();

    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let mut handles: Vec<Handle> = graph.handles_iter().collect();
    handles.sort();
    (graph, handles)
}
