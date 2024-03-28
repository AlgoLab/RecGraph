
use bstr::BString;
use petgraph::graph::NodeIndex;
use petgraph::Graph;
use petgraph::algo::astar;

pub fn demo(graph:&BString, seq: &BString) {
    let mut g = Graph::new();
    
    g.add_node((0, 0));
    for char in 1..seq.len() {
        g.add_node((0, char));
        g.add_edge( NodeIndex::new(char-1), NodeIndex::new(char), 1);
    }
    for node in 1..graph.len() {
        for char in 0..seq.len() {
            g.add_node((node, char));
            if char == 0 {
                g.add_edge( NodeIndex::new((node-1)*seq.len()+(char)), NodeIndex::new(node*seq.len()+char), 1);
            } else {
                let weight = if graph[char] == seq[node] {
                    0
                } else {
                    1
                };
                g.add_edge( NodeIndex::new((node-1)*seq.len()+(char-1)), NodeIndex::new(node*seq.len()+char), weight);
                g.add_edge( NodeIndex::new((node-1)*seq.len()+(char)), NodeIndex::new(node*seq.len()+char), 1);
                g.add_edge( NodeIndex::new((node)*seq.len()+(char-1)), NodeIndex::new(node*seq.len()+char), 1);
            }
            

        }
    }
    
    let start = g.node_indices().nth(0).unwrap();
    let end = g.node_indices().last().unwrap();

    let path = astar(&g, start, |finish| finish == end, |e| *e.weight(), |_| 0);
    println!("{:?}", path)
}

