/// Application program interface of rsPOA
pub mod api;
/// Command Line Interface
pub mod args_parser;
/// Interface for path managment
pub mod bitfield_path;
/// .gaf file creation
pub mod gaf_output;
/// LnzGraph creation
pub mod graph;
/// Score matrix for each alignment type
pub mod score_matrix;
/// DEMO
pub mod pathwise_alignment;
/// Read preparation for POA
pub mod sequences;
/// Various and miscellaneous
pub mod utils;
pub mod gap_global_abpoa;
pub mod gap_local_poa;
pub mod global_abpoa;
pub mod local_poa;
