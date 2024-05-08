use crate::new_path_graph::path_graph::PathGraph;
use ahash::AHashMap as HashMap;
use bstr::BString;
use pheap::PairingHeap as FibHeap;
use std::hash::Hash;
use std::{fmt::Debug, hash::Hasher};

pub fn exec(
    query: &BString,
    heuristic: &(Vec<Vec<u32>>, Vec<HashMap<(usize, usize), (usize, usize)>>),
    path_graph: &PathGraph,
    is_local: bool,
    rec_cost: u32,
) -> (Coord, HashMap<Coord, AStarNode>) {
    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut open_set = FibHeap::new();
    let (crumbs, match_handles) = heuristic;
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path, &crumbs);
        let node_coord = Coord::init(0, 0, path as u8);
        open_set.insert(node_coord.clone(), node.g + node.h);
        alignment_graph.insert(node_coord, node);
    }

    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        let (current_node_coord, _) = open_set.delete_min().unwrap();
        let current_node = alignment_graph.get(&current_node_coord).unwrap().clone();
        if current_node_coord.pos == query.len() as u32 - 1
            && (current_node_coord.node
                == path_graph.ending_positions[current_node_coord.path as usize] as u32
                || is_local)
        {
            end_pos = Some(current_node_coord);
            break;
        }

        if current_node_coord.node + 1 < path_graph.lnz.len() as u32
            && current_node_coord.pos + 1 < query.len() as u32
        {
            if false {
                let mut skip_ahead = current_node.clone();
                skip_ahead.parent = current_node_coord.clone();

                let (skip_ahead_node, skip_ahead_pos) = match_handles
                    [current_node_coord.path as usize]
                    .get(&(
                        current_node_coord.node as usize,
                        current_node_coord.pos as usize,
                    ))
                    .unwrap();
                let skip_ahead_coord = Coord::init(
                    *skip_ahead_node as u32,
                    *skip_ahead_pos as u32,
                    current_node_coord.path,
                );

                skip_ahead.h =
                    crumbs[skip_ahead_coord.path as usize][skip_ahead_coord.pos as usize];
                update_open_set(
                    &mut open_set,
                    &mut alignment_graph,
                    &skip_ahead,
                    skip_ahead_coord,
                )
            } else {
                if !path_graph.nws[current_node_coord.node as usize] {
                    let match_mis = if path_graph.lnz[current_node_coord.node as usize + 1]
                        == query[current_node_coord.pos as usize + 1]
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
                        is_local,
                    );
                } else {
                    path_graph
                        .succ_hash
                        .get_node_succs(current_node_coord.node)
                        .iter()
                        .for_each(|succ| {
                            let paths = path_graph.get_node_path(current_node_coord.node);
                            if paths[current_node_coord.path as usize] {
                                let match_mis = if path_graph.lnz[*succ as usize]
                                    == query[current_node_coord.pos as usize + 1]
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
                                    *succ,
                                    is_local,
                                );
                            }
                        });
                }
                new_multi_rec(
                    &current_node,
                    &current_node_coord,
                    &crumbs,
                    &mut open_set,
                    &mut alignment_graph,
                    &path_graph,
                    rec_cost,
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
    crumbs: &Vec<Vec<u32>>,
    match_mis: u32,
    rec_path: Option<(u32, u32)>, // path, rec_cost
) -> (AStarNode, AStarNode, AStarNode) {
    if let Some((new_path, rec_score)) = rec_path {
        let m_x = AStarNode::init(
            current_node.g + rec_score + match_mis,
            crumbs[new_path as usize][current_node_coord.pos as usize + 1],
            &current_node_coord,
        );

        let ins = AStarNode::init(
            current_node.g + rec_score + 1,
            crumbs[new_path as usize][current_node_coord.pos as usize],
            &current_node_coord,
        );

        let del = AStarNode::init(
            current_node.g + rec_score + 1,
            crumbs[new_path as usize][current_node_coord.pos as usize + 1],
            &current_node_coord,
        );
        (m_x, ins, del)
    } else {
        let h = crumbs[current_node_coord.path as usize][current_node_coord.pos as usize + 1];

        let m_x = AStarNode::init(current_node.g + match_mis, h, &current_node_coord);

        let ins = AStarNode::init(current_node.g + 1, current_node.h, &current_node_coord);

        let del = AStarNode::init(current_node.g + 1, h, &current_node_coord);

        (m_x, ins, del)
    }
}

fn update_open_set(
    open_set: &mut FibHeap<Coord, u32>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    new_node: &AStarNode,
    new_node_coord: Coord,
) {
    let old_node = alignment_graph.get(&new_node_coord);
    if old_node.is_none() || old_node.unwrap().g > new_node.g {
        open_set.insert(new_node_coord, new_node.g + new_node.h);
        alignment_graph.insert(new_node_coord, new_node.clone());
    }
}

fn new_multi_rec(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    crumbs: &Vec<Vec<u32>>,
    open_set: &mut FibHeap<Coord, u32>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    path_graph: &PathGraph,
    rec_cost: u32,
) {
    let paths = path_graph.get_node_path(current_node_coord.node);
    paths.iter().enumerate().for_each(|(path, is_in)| {
        let rec_h = crumbs[path][current_node_coord.pos as usize];
        if is_in && path != current_node_coord.path as usize && rec_h <= current_node.h {
            let rec_node = AStarNode::init(current_node.g + rec_cost, rec_h, &current_node_coord);
            update_open_set(
                open_set,
                alignment_graph,
                &rec_node,
                Coord::init(current_node_coord.node, current_node_coord.pos, path as u8),
            );
        }
    })
}
/*
fn add_multi_recs(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    crumbs: &Vec<Vec<u32>>,
    open_set: &mut FibHeap<Coord, u32>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    path_graph: &PathGraph,
    query: &BString,
    rec_cost: u32,
) {
    let succs = path_graph
        .succ_hash
        .get_node_succs(current_node_coord.node);
    succs.iter().for_each(|succ| {
        let paths = path_graph.get_node_path(current_node_coord.node);
        paths.iter().enumerate().for_each(|(path, is_in)| {
            if is_in
                && path as u32 != current_node_coord.path
                && crumbs[path as usize][current_node_coord.pos as usize]
                    <= crumbs[current_node_coord.path as usize][current_node_coord.pos as usize]
            {
                let match_mis = if path_graph.lnz[*succ as usize]
                    == query[current_node_coord.pos as usize + 1]
                {
                    0
                } else {
                    1
                };
                let (m_x, ins, del) = get_neighbours(
                    &current_node,
                    &current_node_coord,
                    &crumbs,
                    match_mis,
                    Some((path as u32, rec_cost)),
                );

                update_open_set(
                    open_set,
                    alignment_graph,
                    &m_x,
                    Coord::init(*succ, current_node_coord.pos + 1, path as u32),
                );

                update_open_set(
                    open_set,
                    alignment_graph,
                    &ins,
                    Coord::init(*succ, current_node_coord.pos, path as u32),
                );

                update_open_set(
                    open_set,
                    alignment_graph,
                    &del,
                    Coord::init(
                        current_node_coord.node,
                        current_node_coord.pos + 1,
                        path as u32,
                    ),
                );
            }
        });
    });
}
 */
fn push_neigh(
    match_mis: u32,
    current_node: &AStarNode,
    current_node_coord: &Coord,
    open_set: &mut FibHeap<Coord, u32>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    crumbs: &Vec<Vec<u32>>,
    succ: u32,
    is_local: bool,
) {
    let (m_x, mut ins, del) =
        get_neighbours(&current_node, &current_node_coord, &crumbs, match_mis, None);

    if current_node_coord.pos == 0 && is_local {
        ins.g = 0;
    }

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
#[derive(Debug, Clone)]
pub struct Coord {
    pub node: u32,
    pub pos: u32,
    pub path: u8,
}

impl Coord {
    pub fn new() -> Self {
        Coord {
            node: 0,
            pos: 0,
            path: 0,
        }
    }

    pub fn init(node: u32, pos: u32, path: u8) -> Self {
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
    pub g: u32,
    pub h: u32,
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

    pub fn new_path(path: usize, heu: &Vec<Vec<u32>>) -> Self {
        let mut node = AStarNode::new();
        node.h = heu[path][0];
        node
    }

    pub fn init(g: u32, h: u32, parent: &Coord) -> Self {
        AStarNode {
            g,
            h,
            parent: parent.clone(),
        }
    }
}
