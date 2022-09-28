//! rsPOA is an implementation of the POA algorithm, used to align a sequence to a pangenome graph.   
//! More information at [rsPOA](https://github.com/AlgoLab/rspoa)

/// Application program interface of rsPOA
pub mod api;
/// Command Line Interface
pub mod args_parser;
/// Interface for path managment
pub mod bitfield_path;
/// .gaf file creation
pub mod gaf_output;
/// adaptive banded POA with gap opening and penalty
pub mod gap_global_abpoa;
/// adaptive banded local POA with gap opening and penalty
pub mod gap_local_poa;
/// adaptive banded POA, with avx2 instructions
pub mod global_abpoa;
/// LnzGraph creation
pub mod graph;
/// adaptive banded local POA, with avx2 instructions
pub mod local_poa;
/// DEMO
pub mod pathwise_alignment;
pub mod pathwise_alignment_gap;
pub mod pathwise_alignment_gap_semi;
pub mod pathwise_alignment_output;
pub mod pathwise_alignment_recombination;
pub mod pathwise_alignment_semiglobal;
/// Pathwise graph creation
pub mod pathwise_graph;
pub mod recombination_output;
/// Score matrix for each alignment type
pub mod score_matrix;
/// Read preparation for POA
pub mod sequences;
/// Various and miscellaneous
pub mod utils;
