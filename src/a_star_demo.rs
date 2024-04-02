
use bstr::BString;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use lt_fm_index::LtFmIndexBuilder;
use petgraph::graph::NodeIndex;
use petgraph::Graph;
use petgraph::algo::astar;

use crate::args_parser::ClArgs;

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

pub fn fm_index_demo(query: &BString) {
    let file_path = ClArgs::parse().graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let mut paths_seq = Vec::new();
    graph.paths.iter().for_each(|(id,path)|{
        let mut path_seq = String::new();   
        path.nodes.iter().for_each(|node| {
            let node_seq = graph.sequence(*node);
            path_seq.push_str(&String::from_utf8(node_seq).unwrap());
        });
        paths_seq.push((id, path_seq));
    });

    
    for (id, seq) in paths_seq {
        let builder = LtFmIndexBuilder::new()
            .text_type_is_inferred()
            .set_lookup_table_kmer_size_to_default()
            .set_suffix_array_sampling_ratio_to_default();
        let fm_index = builder.build(seq.as_bytes().to_vec()).unwrap();
        let seeds: Vec<&[u8]> = query.chunks(10).collect();
        let mut path_matches = Vec::new();
        for seed in seeds {
            let matches = fm_index.locate(seed);
            path_matches.push(matches);
        }
        println!("Path: {}", id);
        println!("\tMatches: {:?}", path_matches);
        
    }
}