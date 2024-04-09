use bstr::BString;
use handlegraph::{handlegraph::HandleGraph, hashgraph::{HashGraph, Path}};
use lt_fm_index::{LtFmIndex, LtFmIndexBuilder};
use rayon::prelude::*;

const CHUNK_SIZE: usize = 50; // change this later

/// Get the FM index for each path in the graph, return Vec[(path_id, fm_index); paths_number]
fn get_fm_index(graph: &HashGraph) -> Vec<(usize, LtFmIndex)> {
    let mut fm_indexes: Vec<(usize, LtFmIndex)>= graph.paths.par_iter().map(|(id, path)| {
        let path_seq = linearize_path(path, graph);
        let builder = LtFmIndexBuilder::new()
            .text_type_is_inferred()
            .set_lookup_table_kmer_size_to_default()
            .set_suffix_array_sampling_ratio_to_default();
        let fm_index = builder.build(path_seq).unwrap();
        (*id as usize, fm_index)
    }).collect();
    fm_indexes.sort_by_key(|(id, _)| *id);  
    fm_indexes
}

fn linearize_path(path: &Path, graph: &HashGraph)  -> Vec<u8>{
    path.nodes.iter().map(|node| {
        graph.sequence(*node)
    }).collect::<Vec<_>>().concat()
}

/// Get the matches for each path in the graph, return Vec[Vec[bool; seeds_number]; paths_number]
/// matches[i][j] = true if the j-th seed matches the i-th path in some position
fn get_matches(graph: &HashGraph, query: &BString) -> Vec<Vec<bool>> {
    let indexes = get_fm_index(graph);
    
    let seeds: Vec<_> = query.chunks(CHUNK_SIZE).collect::<Vec<_>>();
    
    let mut matches = vec![vec![false; seeds.len()]; indexes.len()];

    indexes.iter().for_each(|(path_id, index)| {
        seeds.iter().enumerate().for_each(|(seed_id, seed)| {
            let is_match = index.count(seed) > 0;
            matches[*path_id][seed_id] = is_match;
        });
    });
    matches
}

/// Get the base heuristic for each path in the graph, return Vec[Vec[usize; query.len()]; paths_number]
/// heuristic[i][j] = x, where x is the number of seeds after the j-th that match the i-th path 
pub fn get_base_sh(query: &BString, graph: &HashGraph) ->  Vec<Vec<usize>>{
    let matches = get_matches(graph, query);
    let seeds_number = query.len()/CHUNK_SIZE;
    let paths_number = graph.paths.len();
    let mut heuristic = vec![vec![0; query.len()];paths_number];
    heuristic.iter_mut().enumerate().for_each(|(path_id, path_heu)| {
        let path_matches = &matches[path_id];
        path_heu.iter_mut().enumerate().for_each(|(pos, val)| {
            let idx = pos/CHUNK_SIZE;
            let potential = seeds_number - idx;
            let actual = potential -  path_matches[idx+1..].iter().filter(|&&x| x).count();
            *val = actual;

        })
    });

    heuristic
}

#[derive(Debug)]
pub struct GraphMatches {
    pub matches: Vec<PathMatches>,
    pub seeds: Vec<BString>,
    pub paths: Vec<usize>,
}


#[derive(Debug)]
pub struct PathMatches {
    pub path_id: usize,
    pub matches: Vec<bool>,
}

impl PathMatches {
    pub fn new() -> Self {
        PathMatches {
            path_id: 0,
            matches: Vec::new(),
        }
    }
    pub fn init(path_id: usize, matches: Vec<bool>) -> Self {
        PathMatches {
            path_id,
            matches,
        }
    }

    pub fn set_position(&mut self, pos: usize) {
        self.matches[pos] = true;
    }

    pub fn get_position(&self, pos: usize) -> bool {
        self.matches[pos]
    }
}