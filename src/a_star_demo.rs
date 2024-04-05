use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::f32::INFINITY;
use std::path;

use bstr::BString;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::handlegraph::HandleGraph;
use handlegraph::hashgraph::HashGraph;
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};
use needletail::parser;
use petgraph::algo::astar;
use petgraph::graph::NodeIndex;
use petgraph::Graph;

use crate::args_parser::ClArgs;
use crate::pathwise_graph::create_path_graph;

pub fn demo(graph: &BString, seq: &BString) {
    let mut g = Graph::new();

    g.add_node((0, 0));
    for char in 1..seq.len() {
        g.add_node((0, char));
        g.add_edge(NodeIndex::new(char - 1), NodeIndex::new(char), 1);
    }
    for node in 1..graph.len() {
        for char in 0..seq.len() {
            g.add_node((node, char));
            if char == 0 {
                g.add_edge(
                    NodeIndex::new((node - 1) * seq.len() + (char)),
                    NodeIndex::new(node * seq.len() + char),
                    1,
                );
            } else {
                let weight = if graph[char] == seq[node] { 0 } else { 1 };
                g.add_edge(
                    NodeIndex::new((node - 1) * seq.len() + (char - 1)),
                    NodeIndex::new(node * seq.len() + char),
                    weight,
                );
                g.add_edge(
                    NodeIndex::new((node - 1) * seq.len() + (char)),
                    NodeIndex::new(node * seq.len() + char),
                    1,
                );
                g.add_edge(
                    NodeIndex::new((node) * seq.len() + (char - 1)),
                    NodeIndex::new(node * seq.len() + char),
                    1,
                );
            }
        }
    }

    let start = g.node_indices().nth(0).unwrap();
    let end = g.node_indices().last().unwrap();

    let path = astar(&g, start, |finish| finish == end, |e| *e.weight(), |_| 0);
    
    
}

pub fn a_star_demo(query: &BString) -> AStarNode{
    let file_path = ClArgs::parse().graph_path;
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(file_path).unwrap();
    let graph: HashGraph = HashGraph::from_gfa(&gfa);
    let path_graph = create_path_graph(&graph, false);

    // compute basic heuristic
    let crumbs = get_base_sh(query, &graph);

    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut open_set = BinaryHeap::new();
    
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path);
        open_set.push(node.clone());
        alignment_graph.insert((node.node, node.pos, node.path),node); 
    }


    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        
        let current_node = open_set.pop().unwrap();
        if current_node.pos == query.len() -1 {
            println!("FOUND");
            end_pos = Some(current_node);
            break;
        } 
        // add neigh of current node (EDIT ops + rec) if not outside graph
        if current_node.node+1 < path_graph.lnz.len() && current_node.pos+1 < query.len() {

        
        let (m_x, ins, del) = if !path_graph.nwp_rev[current_node.node]{
            let match_mis = if path_graph.lnz[current_node.node+1] == query[current_node.pos+1] {
                0
            } else {
                1
            };
            let h = crumbs[current_node.path][current_node.pos+1];

            (AStarNode::init(
                current_node.node + 1,
                current_node.pos + 1,
                current_node.path,
                current_node.g + match_mis,
                h,
                &current_node
            ),
            
            AStarNode::init(
                current_node.node + 1,
                current_node.pos,
                current_node.path,
                current_node.g + 1,
                current_node.h,
                &current_node
            ),

            AStarNode::init(
                current_node.node,
                current_node.pos + 1,
                current_node.path,
                current_node.g + 1,
                h,
                &current_node
            )
        )
        } else {
            let (mut m_x, mut ins, mut del) = (AStarNode::new(), AStarNode::new(), AStarNode::new());
            path_graph.pred_hash_rev.get_preds_and_paths(current_node.node).iter().for_each(|(succ, paths)|{
                if paths[current_node.path] {
                    let match_mis = if path_graph.lnz[*succ] == query[current_node.pos+1] {
                        0
                    } else {
                        1
                    };
                    let h = crumbs[current_node.path][current_node.pos+1];
        
                    m_x = AStarNode::init(
                        *succ,
                        current_node.pos + 1,
                        current_node.path,
                        current_node.g + match_mis,
                        h,
                        &current_node
                    );
                    
                    ins = AStarNode::init(
                        *succ,
                        current_node.pos,
                        current_node.path,
                        current_node.g + 1,
                        current_node.h,
                        &current_node
                    );
        
                   del = AStarNode::init(
                        current_node.node,
                        current_node.pos + 1,
                        current_node.path,
                        current_node.g + 1,
                        h,
                        &current_node
                    );
                                      
                }
            }); 
            (m_x, ins, del)  
        };

        let old_m_x = alignment_graph.remove(&(m_x.node,m_x.pos,m_x.path));
        let old_ins = alignment_graph.remove(&(ins.node,ins.pos,ins.path));
        let old_del = alignment_graph.remove(&(del.node,del.pos,del.path));

        if old_m_x.is_none() {
            open_set.push(m_x.clone());
            alignment_graph.insert((m_x.node,m_x.pos,m_x.path), m_x.clone());
        } else if m_x.g < old_m_x.as_ref().unwrap().g{
            open_set.push(m_x.clone());
            alignment_graph.insert((m_x.node,m_x.pos,m_x.path), m_x.clone());
        } else {
            alignment_graph.insert((m_x.node,m_x.pos,m_x.path), old_m_x.unwrap());
        }

        if old_ins.is_none() {
            open_set.push(ins.clone());
            alignment_graph.insert((ins.node,ins.pos,ins.path), ins.clone());
        } else if ins.g < old_ins.as_ref().unwrap().g{
            open_set.push(ins.clone());
            alignment_graph.insert((ins.node,ins.pos,ins.path), ins.clone());
        } else {
            alignment_graph.insert((ins.node,ins.pos,ins.path), old_ins.unwrap());
        }

        if old_del.is_none() {
            open_set.push(del.clone());
            alignment_graph.insert((del.node,del.pos,del.path), del.clone());
        } else if del.g < old_del.as_ref().unwrap().g{
            open_set.push(del.clone());
            alignment_graph.insert((del.node,del.pos,del.path), del.clone());
        } else {
            alignment_graph.insert((del.node,del.pos,del.path), old_del.unwrap());
        }

    }
    }
    
     
    let mut align = end_pos.unwrap();
    while  align.parent != (0,0,0) {
        align = alignment_graph.remove(&(align.parent.0,align.parent.1,align.parent.2,)).unwrap();
        println!("{:?}  {} {}", align, path_graph.lnz[align.node], query[align.pos]);  
      }
    
    align
}

fn get_fm_index(graph: &HashGraph) -> Vec<(usize, LtFmIndex)> {
    let mut paths_index = Vec::new();

    graph.paths.iter().for_each(|(id, path)| {
        let mut path_seq = String::new();
        path.nodes.iter().for_each(|node| {
            let node_seq = graph.sequence(*node);
            path_seq.push_str(&String::from_utf8(node_seq).unwrap());
        });
        let builder = LtFmIndexBuilder::new()
            .text_type_is_inferred()
            .set_lookup_table_kmer_size_to_default()
            .set_suffix_array_sampling_ratio_to_default();
        let fm_index = builder.build(path_seq.as_bytes().to_vec()).unwrap();
        paths_index.push((*id as usize, fm_index));
    });
    paths_index
}

fn get_matches(graph: &HashGraph, query: &BString) -> Vec<Vec<bool>> {
    let indexes = get_fm_index(graph);
    let chunk_size = 50;
    let seeds: Vec<_> = query.chunks(chunk_size).collect::<Vec<_>>();
    
    let mut matches = vec![vec![false; seeds.len()]; indexes.len()];

    indexes.iter().for_each(|(path_id, index)| {
        seeds.iter().enumerate().for_each(|(seed_id, seed)| {
            let is_match = index.count(seed) > 0;
            matches[*path_id][seed_id] = is_match;
        });
    });
    matches
}
fn get_base_sh(query: &BString, graph: &HashGraph) ->  Vec<Vec<usize>>{
    let chunk_size = 50;
    let matches = get_matches(graph, query);
    let seeds_number = query.len()/chunk_size;
    let paths_number = graph.paths.len();
    let mut heuristic = vec![vec![0; query.len()];paths_number];
    heuristic.iter_mut().enumerate().for_each(|(path_id, path_heu)| {
        let path_matches = &matches[path_id];
        path_heu.iter_mut().enumerate().for_each(|(pos, val)| {
            let idx = pos/chunk_size;
            let potential = seeds_number - idx;
            let actual = potential -  path_matches[idx+1..].iter().filter(|&&x| x).count();
            *val = actual;

        })
    });

    heuristic
}

#[derive(Debug, Clone)]
pub struct AStarNode{
    pub node: usize, 
    pub pos: usize,
    pub path: usize, 
    pub g: usize,
    pub h: usize,
    pub parent: (usize, usize, usize)
}
impl AStarNode {
    pub fn new() -> Self{
        AStarNode{
            node:0,
            pos: 0,
            path: 0,
            g: 0,
            h: 0,
            parent: (0,0,0)
        }
    }

    pub fn new_path(path: usize) -> Self {
        let mut node = AStarNode::new();
        node.path = path;
        node
    }

    pub fn init(node: usize, pos: usize, path: usize, g:usize, h: usize, parent: &AStarNode) -> Self {
        AStarNode{
            node,
            pos, 
            path, 
            g, 
            h,
            parent: (parent.node, parent.pos, parent.path)
        }
    }

    
}

impl Ord for AStarNode {
    fn cmp(&self, other: &Self) -> Ordering {
        // Ordinamento basato sul campo curr_best_score
        (self.g+self.h)
            .partial_cmp(&(other.g + other.h))
            .unwrap_or(Ordering::Equal).reverse()
    }
}

impl PartialOrd for AStarNode  {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for AStarNode  {
    fn eq(&self, other: &Self) -> bool {
        self.node == other.node && self.pos == other.pos && self.path == other.path
}
}

impl Eq for AStarNode  {
}