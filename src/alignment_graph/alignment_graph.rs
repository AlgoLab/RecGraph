use std::{collections::HashSet, rc::Rc};

use ahash::AHashMap as HashMap;
use pheap::PairingHeap;

#[derive(Debug)]
pub struct AlignmentGraph {
    pub coord_set: HashSet<Coord>,
    pub open_set: PairingHeap<Rc<Coord>, usize>,
    pub explored: HashMap<Rc<Coord>, Rc<Coord>>,
}

#[derive(Debug, Clone, Hash, Eq, PartialEq)]
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

impl AlignmentGraph {
    pub fn new() -> Self {
        AlignmentGraph {
            coord_set: HashSet::new(),
            open_set: PairingHeap::new(),
            explored: HashMap::new(),
        }
    }
}
