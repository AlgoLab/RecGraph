use crate::pathwise_graph::PathGraph;
use ahash::AHashMap as HashMap;
use bstr::BString;
use pheap::PairingHeap as FibHeap;
use std::{cmp::Ordering, f32::consts::E, fmt::Debug};
pub fn exec(
    query: &BString,
    crumbs: Vec<Vec<usize>>,
    path_graph: &PathGraph,
) -> (AStarNode, HashMap<(usize, usize, usize), AStarNode>) {
    // init A* data structure, each path possible starting point
    let mut alignment_graph = HashMap::new();
    let mut open_set = FibHeap::new();
    for path in 0..crumbs.len() {
        let node = AStarNode::new_path(path);
        open_set.insert(node.clone(), node.g + node.h);
        alignment_graph.insert((node.node, node.pos, node.path), node);
    }

    // use PathGraph to navigate graph
    let mut end_pos = None;
    while !open_set.is_empty() {
        let (current_node, _) = open_set.delete_min().unwrap();
        if current_node.pos == query.len() - 1 {
            end_pos = Some(current_node);
            break;
        }
        // add neigh of current node (EDIT ops + rec) if not outside graph
        if current_node.node + 1 < path_graph.lnz.len() && current_node.pos + 1 < query.len() {
            if !path_graph.nwp_rev[current_node.node] {
                let match_mis =
                    if path_graph.lnz[current_node.node + 1] == query[current_node.pos + 1] {
                        0
                    } else {
                        1
                    };
                let neigh =
                    get_neighbours(&current_node, &crumbs, match_mis, current_node.node + 1);
                update_open_set(&mut open_set, &mut alignment_graph, &neigh.0);
                update_open_set(&mut open_set, &mut alignment_graph, &neigh.1);
                update_open_set(&mut open_set, &mut alignment_graph, &neigh.2);
            } else {
                path_graph
                    .pred_hash_rev
                    .get_preds_and_paths(current_node.node)
                    .iter()
                    .for_each(|(succ, paths)| {
                        if paths[current_node.path] {
                            let match_mis = if path_graph.lnz[*succ] == query[current_node.pos + 1]
                            {
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
    let (m_x, ins, del) = {
        let h = crumbs[current_node.path][current_node.pos + 1];

        let m_x = AStarNode::init(
            succ,
            current_node.pos + 1,
            current_node.path,
            current_node.g + match_mis,
            h,
            &current_node,
        );

        let ins = AStarNode::init(
            succ,
            current_node.pos,
            current_node.path,
            current_node.g + 1,
            current_node.h,
            &current_node,
        );

        let del = AStarNode::init(
            current_node.node,
            current_node.pos + 1,
            current_node.path,
            current_node.g + 1,
            h,
            &current_node,
        );

        (m_x, ins, del)
    };
    (m_x, ins, del)
}

fn update_open_set(
    open_set: &mut FibHeap<AStarNode, usize>,
    alignment_graph: &mut HashMap<(usize, usize, usize), AStarNode>,
    new_node: &AStarNode,
) {
    if let Some(old_node) = alignment_graph.remove(&(new_node.node, new_node.pos, new_node.path)) {
        if new_node.g < old_node.g {
            open_set.insert(new_node.clone(), new_node.g + new_node.h);
            alignment_graph.insert(
                (new_node.node, new_node.pos, new_node.path),
                new_node.clone(),
            );
        } else {
            alignment_graph.insert((new_node.node, new_node.pos, new_node.path), old_node);
        }
    } else {
        open_set.insert(new_node.clone(), new_node.g + new_node.h);
        alignment_graph.insert(
            (new_node.node, new_node.pos, new_node.path),
            new_node.clone(),
        );
    }
}

#[derive(Debug, Clone)]
pub struct AStarNode {
    pub node: usize,
    pub pos: usize,
    pub path: usize,
    pub g: usize,
    pub h: usize,
    pub parent: (usize, usize, usize),
}
impl AStarNode {
    pub fn new() -> Self {
        AStarNode {
            node: 0,
            pos: 0,
            path: 0,
            g: 0,
            h: 0,
            parent: (0, 0, 0),
        }
    }

    pub fn new_path(path: usize) -> Self {
        let mut node = AStarNode::new();
        node.path = path;
        node
    }

    pub fn init(
        node: usize,
        pos: usize,
        path: usize,
        g: usize,
        h: usize,
        parent: &AStarNode,
    ) -> Self {
        AStarNode {
            node,
            pos,
            path,
            g,
            h,
            parent: (parent.node, parent.pos, parent.path),
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
        self.node == other.node && self.pos == other.pos && self.path == other.path
    }
}

impl Eq for AStarNode {}
