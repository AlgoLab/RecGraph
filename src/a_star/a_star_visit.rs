use crate::new_path_graph::path_graph::PathGraph;
use ahash::AHashMap as HashMap;
use bstr::BString;
use std::cmp::Reverse;
//use pheap::PairingHeap as FibHeap;
use priority_queue::PriorityQueue as FibHeap;
use std::fmt::Debug;
use std::hash::Hash;

pub fn exec(
    query: &BString,
    heuristic: &mut (Vec<Vec<u32>>, Vec<HashMap<(usize, usize), (usize, usize)>>),
    path_graph: &PathGraph,
    is_local: bool,
    rec_cost: u32,
) -> (Coord, HashMap<Coord, AStarNode> /* , Vec<Coord>*/) {
    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut open_set: FibHeap<Coord, Reverse<u32>> = FibHeap::new();
    let crumbs: &mut Vec<Vec<u32>> = &mut heuristic.0;
    let match_handles = &mut heuristic.1;
    //let mut explored_pos = Vec::new();
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path, &crumbs);
        let node_coord = Coord::init(0, 0, path as u8);
        open_set.push(node_coord.clone(), Reverse(node.g + node.h));
        alignment_graph.insert(node_coord, node);
    }

    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        let (current_node_coord, _) = open_set.pop().unwrap();
        let current_node = alignment_graph.get(&current_node_coord).unwrap().clone();
        //explored_pos.push(current_node_coord.clone());
        if current_node_coord.pos == query.len() as u32 - 2
            && (current_node_coord.node
                == path_graph.ending_positions[current_node_coord.path as usize] as u32
                || is_local)
        {
            end_pos = Some(current_node_coord);
            //let perc_explored = (alignment_graph.len() as f32 / (query.len() as f32 *path_graph.get_graph_size() as f32)) * 100.0;
            //println!("Explored {:.2}%", perc_explored);
            break;
        }

        if current_node_coord.node + 1 < path_graph.lnz.len() as u32
            && current_node_coord.pos + 1 < query.len() as u32
        {
            if let Some((skip_ahead_node, skip_ahead_pos)) =
                match_handles[current_node_coord.path as usize].get(&(
                    current_node_coord.node as usize,
                    current_node_coord.pos as usize,
                ))
            {
                let mut skip_ahead = current_node.clone();
                skip_ahead.parent = current_node_coord.clone();
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
                );
                /*
                crumbs.iter_mut().zip(match_handles.iter_mut()).enumerate().for_each(|(path,(crumb, match_handle))| {
                    if path == current_node_coord.path as usize {
                        update_path_heuristic(
                            crumb,
                            match_handle,
                            (
                                current_node_coord.node as usize,
                                current_node_coord.pos as usize,
                            ),
                        );
                    }

                });
                */

                update_path_heuristic(
                    &mut crumbs[skip_ahead_coord.path as usize],
                    &mut match_handles[skip_ahead_coord.path as usize],
                    (
                        current_node_coord.node as usize,
                        current_node_coord.pos as usize,
                    ),
                );
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
                    let succ = path_graph
                        .succ_hash
                        .get_node_succs(current_node_coord.node, current_node_coord.path as u32);
                    let match_mis = if path_graph.lnz[succ as usize]
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
                        succ,
                        is_local,
                    );
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
    (
        end_pos.unwrap().clone(),
        alignment_graph, /*explored_pos*/
    )
}

fn get_neighbours(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    crumbs: &Vec<Vec<u32>>,
    match_mis: u32,
) -> (AStarNode, AStarNode, AStarNode) {
    let h = crumbs[current_node_coord.path as usize][current_node_coord.pos as usize + 1];

    let m_x = AStarNode::init(current_node.g + match_mis, h, &current_node_coord);

    let ins = AStarNode::init(
        current_node.g + 1,
        crumbs[current_node_coord.path as usize][current_node_coord.pos as usize],
        &current_node_coord,
    );

    let del = AStarNode::init(current_node.g + 1, h, &current_node_coord);

    (m_x, ins, del)
}

fn update_open_set(
    open_set: &mut FibHeap<Coord, Reverse<u32>>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    new_node: &AStarNode,
    new_node_coord: Coord,
) {
    if let Some(old_node) = alignment_graph.get_mut(&new_node_coord) {
        if old_node.g > new_node.g {
            open_set.push(new_node_coord, Reverse(new_node.g + new_node.h));
            *old_node = new_node.clone();
        }
    } else {
        open_set.push(new_node_coord, Reverse(new_node.g + new_node.h));
        alignment_graph.insert(new_node_coord, new_node.clone());
    }
}

fn new_multi_rec(
    current_node: &AStarNode,
    current_node_coord: &Coord,
    crumbs: &Vec<Vec<u32>>,
    open_set: &mut FibHeap<Coord, Reverse<u32>>,
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

fn push_neigh(
    match_mis: u32,
    current_node: &AStarNode,
    current_node_coord: &Coord,
    open_set: &mut FibHeap<Coord, Reverse<u32>>,
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    crumbs: &Vec<Vec<u32>>,
    succ: u32,
    is_local: bool,
) {
    let (m_x, mut ins, del) =
        get_neighbours(&current_node, &current_node_coord, &crumbs, match_mis);

    if current_node_coord.pos == 0 && is_local {
        ins.g = 0;
    }

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

    update_open_set(
        open_set,
        alignment_graph,
        &m_x,
        Coord::init(succ, current_node_coord.pos + 1, current_node_coord.path),
    );
}

fn update_path_heuristic(
    crumbs: &mut Vec<u32>,
    match_handles: &mut HashMap<(usize, usize), (usize, usize)>,
    last_match: (usize, usize),
) {
    if match_handles.contains_key(&last_match) {
        crumbs[..last_match.1].iter_mut().for_each(|score| {
            *score += 1;
        });

        match_handles.remove(&last_match);
    }
}

#[derive(Debug, Clone, Hash)]
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
