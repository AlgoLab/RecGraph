//! RecGraph is an exact sequence to variation graph aligner that allows controlled recombinations.
//! More information at [rsPOA](https://github.com/AlgoLab/recgraph)

/// Library Interface
pub mod api;

/// Command Line Interface
pub mod args_parser;

/// .gaf file creation
pub mod gaf_output;

pub mod build_cigar;
/// Dynamic programming matrix
pub mod dp_matrix;

pub mod pathwise_alignment_recombination;
/// Pathwise graph creation
pub mod pathwise_graph;
pub mod recombination_output;
/// Score matrix for each alignment type
pub mod score_matrix;
/// Read preparation for POA
pub mod sequences;
/// Various and miscellaneous
pub mod utils;

pub mod node_displacement;

pub mod a_star;

pub mod new_path_graph;
