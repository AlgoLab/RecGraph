use ahash::AHashMap as HashMap;
use bstr::BString;

use crate::{build_cigar, pathwise_graph::PathGraph};

use super::a_star_visit::AStarNode;

pub fn build_gaf(
    alignment_graph: &mut HashMap<(usize, usize, usize), AStarNode>,
    end_pos: &AStarNode,
    path_graph: &PathGraph,
    query: &BString,
) {
    let mut align = end_pos.clone();
    let ed = align.g;
    let mut cigar = Vec::new();
    while (align.node, align.pos) != (0, 0) {
        println!("{:?}", align);
        if align.parent.0 != align.node {
            if align.parent.1 < align.pos {
                if path_graph.lnz[align.node] == query[align.pos] {
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
        align = alignment_graph
            .remove(&(align.parent.0, align.parent.1, align.parent.2))
            .unwrap();
    }
    cigar.reverse();
    println!(
        "{:?}\tbest path: {}\tED {}",
        build_cigar::build_cigar(&cigar),
        align.path,
        ed
    )
}
