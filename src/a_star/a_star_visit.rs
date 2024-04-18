use crate::args_parser::ClArgs;
use crate::new_path_graph::path_graph::PathGraph;
use ahash::AHashMap as HashMap;
use bstr::BString;
use pheap::PairingHeap as FibHeap;
use std::hash::Hash;
use std::{cmp::Ordering, fmt::Debug, hash::Hasher};

pub fn exec(
    query: &BString,
    crumbs: Vec<Vec<usize>>,
    path_graph: &PathGraph,
) -> (AStarNode, HashMap<Coord, AStarNode>) {
    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut best_score_per_position = HashMap::new();
    let mut open_set = FibHeap::new();
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path, &crumbs);
        open_set.insert(node.clone(), node.g + node.h);
        best_score_per_position.insert((node.coord.node, node.coord.pos), (node.g, path));
        alignment_graph.insert(node.coord, node);
    }

    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        let (current_node, _) = open_set.delete_min().unwrap();
        if current_node.coord.pos == query.len() - 1 {
            end_pos = Some(current_node);
            break;
        }
        // add neigh of current node (EDIT ops + rec) if not outside graph
        if current_node.coord.node + 1 < path_graph.lnz.len()
            && current_node.coord.pos + 1 < query.len()
        {
            if !path_graph.nws[current_node.coord.node] {
                let match_mis = if path_graph.lnz[current_node.coord.node + 1]
                    == query[current_node.coord.pos + 1]
                {
                    0
                } else {
                    1
                };
                let neigh = get_neighbours(
                    &current_node,
                    &crumbs,
                    match_mis,
                    current_node.coord.node + 1,
                );

                update_open_set(&mut open_set, &mut alignment_graph, &neigh.0);
                update_open_set(&mut open_set, &mut alignment_graph, &neigh.1);
                update_open_set(&mut open_set, &mut alignment_graph, &neigh.2);
            } else {
                path_graph
                    .succ_hash
                    .get_node_succs_and_paths(current_node.coord.node)
                    .iter()
                    .for_each(|(succ, paths)| {
                        if paths[current_node.coord.path] {
                            let match_mis =
                                if path_graph.lnz[*succ] == query[current_node.coord.pos + 1] {
                                    0
                                } else {
                                    1
                                };
                            let (m_x, ins, del) =
                                get_neighbours(&current_node, &crumbs, match_mis, *succ);

                            update_open_set(&mut open_set, &mut alignment_graph, &m_x);
                            update_open_set(&mut open_set, &mut alignment_graph, &ins);
                            update_open_set(&mut open_set, &mut alignment_graph, &del);
                        }
                    });
            };

            add_multi_recs(
                &current_node,
                &mut best_score_per_position,
                &crumbs,
                path_graph.succ_hash.paths_number,
                &mut open_set,
                &mut alignment_graph,
            )
        }
    }
    (end_pos.unwrap(), alignment_graph)
}

fn get_neighbours(
    current_node: &AStarNode,
    crumbs: &Vec<Vec<usize>>,
    match_mis: usize,
    succ: usize,
) -> (AStarNode, AStarNode, AStarNode) {
    let h = crumbs[current_node.coord.path][current_node.coord.pos + 1];

    let m_x = AStarNode::init(
        Coord::init(succ, current_node.coord.pos + 1, current_node.coord.path),
        current_node.g + match_mis,
        h,
        &current_node.coord,
    );

    let ins = AStarNode::init(
        Coord::init(succ, current_node.coord.pos, current_node.coord.path),
        current_node.g + 1,
        current_node.h,
        &current_node.coord,
    );

    let del = AStarNode::init(
        Coord::init(
            current_node.coord.node,
            current_node.coord.pos + 1,
            current_node.coord.path,
        ),
        current_node.g + 1,
        h,
        &current_node.coord,
    );

    (m_x, ins, del)
}

fn update_open_set(
    open_set: &mut FibHeap<AStarNode, usize>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    new_node: &AStarNode,
) {
    if let Some(old_node) = alignment_graph.remove(&new_node.coord) {
        if new_node.g < old_node.g {
            open_set.insert(new_node.clone(), new_node.g + new_node.h);
            alignment_graph.insert(new_node.coord, new_node.clone());
        } else {
            alignment_graph.insert(new_node.coord, old_node);
        }
    } else {
        open_set.insert(new_node.clone(), new_node.g + new_node.h);
        alignment_graph.insert(new_node.coord, new_node.clone());
    }
}

// TODO: add ony if path is present in node
fn add_multi_recs(
    current_node: &AStarNode,
    best_score_per_position: &mut HashMap<(usize, usize), (usize, usize)>,
    crumbs: &Vec<Vec<usize>>,
    paths_number: usize,
    open_set: &mut FibHeap<AStarNode, usize>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
) {
    let rec_cost = ClArgs::parse().base_rec_cost as usize;
    if let Some((score, path)) =
        best_score_per_position.get(&(current_node.coord.node, current_node.coord.pos))
    {
        if score + rec_cost < current_node.g {
            let rec_node = AStarNode::init(
                Coord::init(current_node.coord.node, current_node.coord.pos, *path),
                *score + rec_cost, // change +1 to rec
                crumbs[*path][current_node.coord.pos],
                &current_node.coord,
            );
            update_open_set(open_set, alignment_graph, &rec_node);
        } else {
            if score > &current_node.g {
                best_score_per_position.insert(
                    (current_node.coord.node, current_node.coord.pos),
                    (current_node.g, current_node.coord.path),
                );
            }
        }
    } else {
        best_score_per_position.insert(
            (current_node.coord.node, current_node.coord.pos),
            (current_node.g, current_node.coord.path),
        );
        for path in 0..paths_number {
            let rec_node = AStarNode::init(
                Coord::init(current_node.coord.node, current_node.coord.pos, path),
                current_node.g + rec_cost,
                crumbs[path][current_node.coord.pos],
                &current_node.coord,
            );
            update_open_set(open_set, alignment_graph, &rec_node);
        }
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
    pub coord: Coord,
    pub g: usize,
    pub h: usize,
    pub parent: Coord,
}
impl AStarNode {
    pub fn new() -> Self {
        AStarNode {
            coord: Coord::new(),
            g: 0,
            h: 0,
            parent: Coord::new(),
        }
    }

    pub fn new_path(path: usize, heu: &Vec<Vec<usize>>) -> Self {
        let mut node = AStarNode::new();
        node.coord.path = path;
        node.h = heu[path][0];
        node
    }

    pub fn init(coord: Coord, g: usize, h: usize, parent: &Coord) -> Self {
        AStarNode {
            coord,
            g,
            h,
            parent: parent.clone(),
        }
    }
}

impl Ord for AStarNode {
    fn cmp(&self, other: &Self) -> Ordering {
        // Ordinamento basato sul campo curr_best_score
        (self.g + self.h)
            .partial_cmp(&(other.g + other.h))
            .unwrap_or(Ordering::Equal)
            .reverse()
    }
}

impl PartialOrd for AStarNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for AStarNode {
    fn eq(&self, other: &Self) -> bool {
        self.coord.node == other.coord.node
            && self.coord.pos == other.coord.pos
            && self.coord.path == other.coord.path
    }
}

impl Eq for AStarNode {}
