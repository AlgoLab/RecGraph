use ahash::AHashMap as HashMap;
use bstr::BString;
use lt_fm_index::{blocks::Block3, LtFmIndex};

use crate::{build_cigar::build_cigar, new_path_graph::path_graph::PathGraph};

use super::a_star_visit::{AStarNode, Coord};
use std::{io::Write, time::Duration};

pub fn build_gaf(
    alignment_graph: &mut HashMap<Coord, AStarNode>,
    end_pos: &Coord,
    path_graph: &PathGraph,
    query: &BString,
    is_local: bool,
    matches_in_path: &Vec<HashMap<(u32, u32), u32>>,
    match_len: usize,
    indexes: &Vec<(LtFmIndex<u32, Block3<u128>>, Vec<u32>)>,
    name: &BString,
    align_time: Duration,
    original_path_ids: &HashMap<u8, u8>,
    explored_cells: f32,
) -> Gaf {
    let mut align_coord = end_pos.clone();
    let mut align = alignment_graph.remove(&align_coord).unwrap();
    let ed = align.g;
    let mut cigar = Vec::new();
    let mut paths: Vec<_> = vec![(
        *original_path_ids.get(&align_coord.path).unwrap(),
        path_graph.handles_ids[align_coord.node as usize],
        align_coord.pos,
    )];

    let end_pos = align_coord.node;
    let mut path_align = Vec::new();
    let mut residue_matches = 0;
    let mut alignment_len = 0;
    let mut path_len = 0;
    while align_coord.pos != 0 {
        let mut rec = false;
        if align.parent.path != align_coord.path {
            paths.push((
                *original_path_ids.get(&align.parent.path).unwrap(),
                path_graph.handles_ids[align.parent.node as usize],
                align.parent.pos,
            ));
            rec = true;
        } else if align_coord.pos - 1 > align.parent.pos {
            let mut idx = align_coord.pos - align.parent.pos;
            while idx > 0 {
                cigar.push('D');
                idx -= 1;
                residue_matches += 1;
                alignment_len += 1;
                path_len += 1;
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
            alignment_len += 1;

            if align.parent.pos != align_coord.pos {
                if path_graph.lnz[align_coord.node as usize] == query[align_coord.pos as usize] {
                    cigar.push('D');
                    residue_matches += 1;
                } else {
                    cigar.push('d');
                }
            } else {
                cigar.push('U');
            }
            path_align.push(path_graph.handles_ids[align_coord.node as usize]);
            path_len += 1;
        } else {
            alignment_len += 1;
            cigar.push('L');
        }
        let rec_number = if rec {
            align_coord.rec - 1
        } else {
            align_coord.rec
        };
        align_coord = Coord {
            node: align.parent.node,
            pos: align.parent.pos,
            path: align.parent.path,
            rec: rec_number,
        };
        align = alignment_graph.remove(&align_coord).unwrap();
    }

    if !is_local {
        while align_coord.node != 0 {
            cigar.push('U');
            path_align.push(path_graph.handles_ids[align_coord.node as usize]);
            path_len += 1;
            align_coord = Coord {
                node: align.parent.node,
                pos: align.parent.pos,
                path: align.parent.path,
                rec: align_coord.rec,
            };
            align = alignment_graph.remove(&align_coord).unwrap();
            alignment_len += 1;
        }
    }
    cigar.reverse();
    path_align.reverse();
    path_align.dedup();
    paths.reverse();
    let cigar_str = build_cigar(&cigar);
    let comments = format!(
        "{}\t{}\t{}\t{}\t{:.4}%",
        ed,
        cigar_str,
        build_path_composition(&paths, path_graph),
        align_time.as_micros(),
        explored_cells
    );
    let alignment: Vec<(u64, char)> = path_align
        .iter()
        .map(|x| *path_graph.original_handles.get(x).unwrap())
        .collect::<Vec<_>>();
    let start_pos = align_coord.node;
    let path_start = get_node_offset(start_pos, path_graph) as usize;
    Gaf::new(
        name.to_string(),
        query.len() - 2,
        0,
        query.len() - 2,
        '+',
        alignment,
        path_len + path_start + get_node_distance_from_end(end_pos, path_graph) as usize,
        path_start,
        path_len + path_start + 1,
        residue_matches,
        alignment_len,
        255,
        comments,
    )
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

fn build_path_composition(paths: &Vec<(u8, u32, u32)>, path_graph: &PathGraph) -> String {
    paths
        .iter()
        .map(|x| format!("{}:{}:{}", x.0, get_node_handle(x.1, path_graph), x.2))
        .collect::<Vec<String>>()
        .join(",")
}

fn get_node_handle(node: u32, path_graph: &PathGraph) -> u32 {
    path_graph.original_handles.get(&node).unwrap().0 as u32
}

fn get_node_offset(node: u32, path_graph: &PathGraph) -> u32 {
    let mut offset = 0;
    let mut start = node;
    while node > 0
        && path_graph.handles_ids[start as usize] == path_graph.handles_ids[node as usize]
    {
        offset += 1;
        start -= 1;
    }
    offset
}

fn get_node_distance_from_end(node: u32, path_graph: &PathGraph) -> u32 {
    let mut offset = 0;
    let mut start = node;
    while node > 0
        && path_graph.handles_ids[start as usize] == path_graph.handles_ids[node as usize]
    {
        offset += 1;
        start += 1;
    }
    offset
}
pub struct Gaf {
    pub query_name: String,
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub path_matching: Vec<(u64, char)>,
    pub path_len: usize,
    pub path_start: usize,
    pub path_end: usize,
    pub residue_matches: usize,
    pub alignment_len: usize,
    pub mapq: u8,
    pub comments: String,
}

impl Gaf {
    pub fn new(
        query_name: String,
        query_len: usize,
        query_start: usize,
        query_end: usize,
        strand: char,
        path_matching: Vec<(u64, char)>,
        path_len: usize,
        path_start: usize,
        path_end: usize,
        residue_matches: usize,
        alignment_len: usize,
        mapq: u8,
        comments: String,
    ) -> Self {
        Gaf {
            query_name,
            query_len,
            query_start,
            query_end,
            strand,
            path_matching,
            path_len,
            path_start,
            path_end,
            residue_matches,
            alignment_len,
            mapq,
            comments,
        }
    }
    pub fn to_string(&self) -> String {
        let path = self
            .path_matching
            .iter()
            .map(|x| {
                format!(
                    "{}{}",
                    match x.1 {
                        '+' => ">",
                        '-' => "<",
                        _ => panic!("Invalid orientation"),
                    },
                    x.0
                )
            })
            .collect::<Vec<String>>()
            .join("");

        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n@CO\t{}",
            self.query_name,
            self.query_len,
            self.query_start,
            self.query_end,
            self.strand,
            path,
            self.path_len,
            self.path_start,
            self.path_end,
            self.residue_matches,
            self.alignment_len,
            self.mapq,
            self.comments
        )
    }
}
