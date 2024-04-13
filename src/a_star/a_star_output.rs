use ahash::AHashMap as HashMap;
use bstr::BString;

use crate::{build_cigar, pathwise_graph::PathGraph};

use super::a_star_visit::{AStarNode, Coord};

pub fn build_gaf(
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    end_pos: &AStarNode,
    path_graph: &PathGraph,
    query: &BString,
) {
    let mut align = end_pos.clone();
    let ed = align.g;
    let mut cigar = Vec::new();
    while (align.coord.node, align.coord.pos) != (0, 0) {
        if align.parent.node != align.coord.node {
            if align.parent.pos < align.coord.pos {
                if path_graph.lnz[align.coord.node] == query[align.coord.pos] {
                    cigar.push('D');
                } else {
                    cigar.push('d');
                }
            } else {
                cigar.push('U');
            }
        } else {
            cigar.push('L');
        }
        align = alignment_graph.remove(&align.parent).unwrap();
    }
    cigar.reverse();
    println!(
        "{:?}\tbest path: {}\tED {}",
        build_cigar::build_cigar(&cigar),
        align.coord.path,
        ed
    )
}
