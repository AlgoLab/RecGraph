use ahash::AHashMap as HashMap;
use bstr::BString;
use lt_fm_index::LtFmIndex;

use crate::{build_cigar, new_path_graph::path_graph::PathGraph};

use super::a_star_visit::{AStarNode, Coord};
use std::io::Write;

pub fn build_gaf(
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    end_pos: &Coord,
    path_graph: &PathGraph,
    query: &BString,
    is_local: bool,
    matches_in_path: &Vec<HashMap<(u32, u32), u32>>,
    match_len: usize,
    indexes: &Vec<(LtFmIndex, Vec<u32>)>,
) -> String {
    let mut align_coord = end_pos.clone();
    let mut align = alignment_graph.remove(&align_coord).unwrap();
    let ed = align.g;
    let mut cigar = Vec::new();
    let mut recs: Vec<String> = Vec::new();
    let mut path_align = Vec::new();
    while align_coord.pos != 0 {
        if align.parent.path != align_coord.path {
            recs.push(format!(
                "REC paths {} - {}\t pos {:?}",
                align.parent.path,
                align_coord.path,
                (align_coord.node, align_coord.pos)
            ));
        } else if align_coord.pos - 1 > align.parent.pos {
            let mut idx = align_coord.pos - align.parent.pos;
            while idx > 0 {
                cigar.push('D');
                idx -= 1;
            }
            let match_end = matches_in_path[align_coord.path as usize]
                .get(&(align_coord.node, align_coord.pos))
                .unwrap();
            let mut idx = 0;
            while idx < match_len {
                let lnz_pos =
                    indexes[align_coord.path as usize].1[*match_end as usize - idx] as usize;
                path_align.push(path_graph.handles_ids[lnz_pos]);
                idx += 1;
            }
        } else if align.parent.node != align_coord.node {
            if align.parent.pos != align_coord.pos {
                if path_graph.lnz[align_coord.node as usize] == query[align_coord.pos as usize] {
                    cigar.push('D');
                } else {
                    cigar.push('d');
                }
            } else {
                cigar.push('U');
            }
            path_align.push(path_graph.handles_ids[align_coord.node as usize]);
        } else {
            cigar.push('L');
        }
        align_coord = align.parent.clone();
        align = alignment_graph.remove(&align_coord).unwrap();
    }

    if !is_local {
        while align_coord.node != 0 {
            cigar.push('U');
            path_align.push(path_graph.handles_ids[align_coord.node as usize]);
            align_coord = align.parent.clone();
            align = alignment_graph.remove(&align_coord).unwrap();
        }
    }
    cigar.reverse();
    path_align.reverse();
    path_align.dedup();
    let alignment = path_align
        .iter()
        .map(|x| path_graph.original_handles.get(x).unwrap().to_string())
        .collect::<Vec<_>>()
        .join(">");
    let recs_out_string = recs.join("\t");
    let output = format!(
        "{}\t{}\tbest path: {}\tED {}\t{}",
        alignment,
        build_cigar::build_cigar(&cigar),
        align_coord.path,
        ed,
        recs_out_string,
    );

    output
}

pub fn save_coords(outfile: &str, coords: &Vec<Vec<Coord>>) {
    let mut file = std::fs::File::create(outfile).unwrap();
    let out = coords
        .iter()
        .map(|coord| {
            coord
                .iter()
                .map(|c| format!("{}\t{}\t{}", c.node, c.pos, c.path))
                .collect::<Vec<String>>()
                .join("\n")
        })
        .collect::<Vec<String>>();

    file.write_all(out.join("\n\n").as_bytes()).unwrap();
}

pub fn get_matches_end_in_path(
    matches: &HashMap<(u32, u32), (u32, u32, u32)>,
) -> HashMap<(u32, u32), u32> {
    let mut matches_in_path = HashMap::new();
    if matches.len() > 0 {
        matches.iter().for_each(|(_, v)| {
            matches_in_path.insert((v.0, v.1), v.2);
        });
    }
    matches_in_path
}
