use crate::args_parser::ClArgs;
use crate::new_path_graph::path_graph::PathGraph;
use ahash::AHashMap as HashMap;
use bit_vec::BitVec;
use bstr::BString;
use pheap::PairingHeap as FibHeap;
use std::hash::Hash;
use std::{fmt::Debug, hash::Hasher};

pub fn exec(
    query: &BString,
    crumbs: &Vec<Vec<usize>>,
    path_graph: &PathGraph,
) -> (Coord, HashMap<Coord, AStarNode>) {
    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut best_score_per_position = HashMap::new();
    let mut open_set = FibHeap::new();
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path, &crumbs);
        let node_coord = Coord::init(0, 0, path);
        open_set.insert(node_coord.clone(), node.g + node.h);
        best_score_per_position.insert((node_coord.node, node_coord.pos), (node.g, path));
        alignment_graph.insert(node_coord, node);
    }

    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        let (current_node_coord, _) = open_set.delete_min().unwrap();
        let current_node = alignment_graph.get(&current_node_coord).unwrap().clone();
        if current_node_coord.pos == query.len() - 1 {
            end_pos = Some(current_node_coord);
            break;
        }
        // add neigh of current node (EDIT ops + rec) if not outside graph
        if current_node_coord.node + 1 < path_graph.lnz.len()
            && current_node_coord.pos + 1 < query.len()
        {
            if !path_graph.nws[current_node_coord.node] {
                let match_mis = if path_graph.lnz[current_node_coord.node + 1]
                    == query[current_node_coord.pos + 1]
                {
                    0
                } else {
                    1
                };

                push_neigh(
                    match_mis,
                    &current_node,
                    &current_node_coord,
                    &mut open_set,
                    &mut alignment_graph,
                    &crumbs,
                    current_node_coord.node + 1,
                );
            } else {
                path_graph
                    .succ_hash
                    .get_node_succs_and_paths(current_node_coord.node)
                    .iter()
                    .for_each(|(succ, paths)| {
                        if paths[current_node_coord.path] {
                            let match_mis =
                                if path_graph.lnz[*succ] == query[current_node_coord.pos + 1] {
                                    0
                                } else {
                                    1
                                };

                            push_neigh(
                                match_mis,
                                &current_node,
                                &current_node_coord,
                                &mut open_set,
                                &mut alignment_graph,
                                &crumbs,
                                *succ,
                            );
                        }
                    });
                let paths = path_graph.succ_hash.get_paths_node(current_node_coord.node);
                add_multi_recs(
                    &current_node,
                    &current_node_coord,
                    &mut best_score_per_position,
                    &crumbs,
                    paths,
                    &mut open_set,
                    &mut alignment_graph,
                );
            }
        }
    }
    if end_pos.is_none() {
        panic!("No path found");
    }
    (end_pos.unwrap().clone(), alignment_graph)
}

fn get_neighbours(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    crumbs: &Vec<Vec<usize>>,
    match_mis: usize,
) -> (AStarNode, AStarNode, AStarNode) {
    let h = crumbs[current_node_coord.path][current_node_coord.pos + 1];

    let m_x = AStarNode::init(current_node.g + match_mis, h, &current_node_coord);

    let ins = AStarNode::init(current_node.g + 1, current_node.h, &current_node_coord);

    let del = AStarNode::init(current_node.g + 1, h, &current_node_coord);

    (m_x, ins, del)
}

fn update_open_set(
    open_set: &mut FibHeap<Coord, usize>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    new_node: &AStarNode,
    new_node_coord: Coord,
) {
    if let Some(old_node) = alignment_graph.remove(&new_node_coord) {
        if new_node.g < old_node.g {
            open_set.insert(new_node_coord.clone(), new_node.g + new_node.h);
            alignment_graph.insert(new_node_coord.clone(), new_node.clone());
        } else {
            alignment_graph.insert(new_node_coord.clone(), old_node);
        }
    } else {
        open_set.insert(new_node_coord.clone(), new_node.g + new_node.h);
        alignment_graph.insert(new_node_coord.clone(), new_node.clone());
    }
}

// TODO: add ony if path is present in node
fn add_multi_recs(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    best_score_per_position: &mut HashMap<(usize, usize), (usize, usize)>,
    crumbs: &Vec<Vec<usize>>,
    paths: &BitVec,
    open_set: &mut FibHeap<Coord, usize>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
) {
    let rec_cost = ClArgs::parse().base_rec_cost as usize;
    if let Some((score, path)) =
        best_score_per_position.get(&(current_node_coord.node, current_node_coord.pos))
    {
        if score + rec_cost < current_node.g {
            let rec_node = AStarNode::init(
                *score + rec_cost,
                crumbs[*path][current_node_coord.pos],
                &current_node_coord,
            );
            let rec_node_coord =
                Coord::init(current_node_coord.node, current_node_coord.pos, *path);
            update_open_set(open_set, alignment_graph, &rec_node, rec_node_coord);
        } else {
            if score > &current_node.g {
                best_score_per_position.insert(
                    (current_node_coord.node, current_node_coord.pos),
                    (current_node.g, current_node_coord.path),
                );
            }
        }
    } else {
        best_score_per_position.insert(
            (current_node_coord.node, current_node_coord.pos),
            (current_node.g, current_node_coord.path),
        );
        for (path, is_in) in paths.iter().enumerate() {
            if is_in && crumbs[path][current_node_coord.pos] <= current_node.h {
                let rec_node = AStarNode::init(
                    current_node.g + rec_cost,
                    crumbs[path][current_node_coord.pos],
                    &current_node_coord,
                );
                let rec_node_coord =
                    Coord::init(current_node_coord.node, current_node_coord.pos, path);
                update_open_set(open_set, alignment_graph, &rec_node, rec_node_coord);
            }
        }
    }
}

fn push_neigh(
    match_mis: usize,
    current_node: &AStarNode,
    current_node_coord: &Coord,
    open_set: &mut FibHeap<Coord, usize>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    crumbs: &Vec<Vec<usize>>,
    succ: usize,
) {
    let (m_x, ins, del) = get_neighbours(&current_node, &current_node_coord, &crumbs, match_mis);
    // change for semiglobal alignment
    if match_mis == 0 {
        update_open_set(
            open_set,
            alignment_graph,
            &m_x,
            Coord::init(succ, current_node_coord.pos + 1, current_node_coord.path),
        );
    } else {
        update_open_set(
            open_set,
            alignment_graph,
            &m_x,
            Coord::init(succ, current_node_coord.pos + 1, current_node_coord.path),
        );
        update_open_set(
            open_set,
            alignment_graph,
            &ins,
            Coord::init(succ, current_node_coord.pos, current_node_coord.path),
        );
        update_open_set(
            open_set,
            alignment_graph,
            &del,
            Coord::init(
                current_node_coord.node,
                current_node_coord.pos + 1,
                current_node_coord.path,
            ),
        );
    }
}
#[derive(Debug, Clone)]
pub struct Coord {
    pub node: usize,
    pub pos: usize,
    pub path: usize,
}

impl Coord {
    pub fn new() -> Self {
        Coord {
            node: 0,
            pos: 0,
            path: 0,
        }
    }

    pub fn init(node: usize, pos: usize, path: usize) -> Self {
        Coord { node, pos, path }
    }
}

impl Hash for Coord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node.hash(state);
        self.pos.hash(state);
        self.path.hash(state);
    }
}

impl PartialEq for Coord {
    fn eq(&self, other: &Self) -> bool {
        self.node == other.node && self.pos == other.pos && self.path == other.path
    }
}

impl Eq for Coord {}

impl Copy for Coord {}

#[derive(Debug, Clone)]
pub struct AStarNode {
    pub g: usize,
    pub h: usize,
    pub parent: Coord,
}
impl AStarNode {
    pub fn new() -> Self {
        AStarNode {
            g: 0,
            h: 0,
            parent: Coord::new(),
        }
    }

    pub fn new_path(path: usize, heu: &Vec<Vec<usize>>) -> Self {
        let mut node = AStarNode::new();
        node.h = heu[path][0];
        node
    }

    pub fn init(g: usize, h: usize, parent: &Coord) -> Self {
        AStarNode {
            g,
            h,
            parent: parent.clone(),
        }
    }
}
