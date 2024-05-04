use ahash::AHashMap as HashMap;
use bstr::BString;

use crate::{build_cigar, new_path_graph::path_graph::PathGraph};

use super::a_star_visit::{AStarNode, Coord};

pub fn build_gaf(
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    end_pos: &Coord,
    path_graph: &PathGraph,
    query: &BString,
    is_local: bool,
) -> String {
    let mut align_coord = end_pos.clone();
    let mut align = alignment_graph.remove(&align_coord).unwrap();
    let mut last_g = -1;
    let ed = align.g;
    let mut cigar = Vec::new();
    let mut recs = Vec::new();
    while align_coord.pos != 0 {
        if align.parent.path != align_coord.path {
            recs.push(format!(
                "REC paths {} - {}\t pos {:?}",
                align.parent.path,
                align_coord.path,
                (align_coord.node, align_coord.pos)
            ));
        }
        if align.parent.node != align_coord.node {
            if align.parent.pos != align_coord.pos {
                if path_graph.lnz[align_coord.node as usize] == query[align_coord.pos as usize] {
                    cigar.push('D');
                } else {
                    cigar.push('d');
                }
            } else {
                cigar.push('U');
            }
        } else {
            if last_g == -1 || last_g != align.g as i32 {
                cigar.push('L');
            } else {
                cigar.push('D');
            }
        }
        align_coord = align.parent.clone();
        last_g = align.g as i32;
        align = alignment_graph.remove(&align_coord).unwrap();
    }

    if !is_local {
        while align_coord.node != 0 {
            cigar.push('U');
            align_coord = align.parent.clone();
            align = alignment_graph.remove(&align_coord).unwrap();
        }
    }
    cigar.reverse();

    let recs_out_string = recs.join("\t");
    let output = format!(
        "{}\tbest path: {}\tED {}\t{}",
        build_cigar::build_cigar(&cigar),
        align_coord.path,
        ed,
        recs_out_string
    );
    output
}
