use bstr::BString;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;

use crate::args_parser::ClArgs;
use crate::pathwise_graph::create_path_graph;

use super::a_star_visit::{self, AStarNode};
use super::matches::get_base_sh;


pub fn a_star_demo(query: &BString) -> AStarNode{
    let file_path = ClArgs::parse().graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = create_path_graph(&graph, false);

    // compute basic heuristic
    let crumbs = get_base_sh(query, &graph);

    
    // navigate graph
    let (end_pos, mut alignment_graph) = a_star_visit::exec(query, crumbs, &path_graph);
    
    let mut align = end_pos;    
    while  align.parent != (0,0,0) {
        println!("{:?}  {} {}", align, path_graph.lnz[align.node], query[align.pos]);
        align = alignment_graph.remove(&(align.parent.0,align.parent.1,align.parent.2,)).unwrap();
         
      }
    align
}




