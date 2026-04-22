use clap::{builder::OsStr, Parser};
use serde::Serialize;

#[derive(Parser, Debug)]
#[clap(author = "Davide Monti <d.monti11@campus.unimib.it>", version, about = "RecGraph", long_about = None)]
struct Args {
    // Input graph
    #[clap(
        help_heading = "I/O",
        short = 'g',
        long = "graph",
        help = "Input graph (in .gfa format) (- for stdin)",
        //required = true,
    )]
    graph_path: String,

    // Input sequence
    #[clap(
        help_heading = "I/O",
        short = 'q',
        long = "query",
        help = "Query reads (in .fasta/.fastq format) (- for stdin)",
        //required = true,
    )]
    sequence_path: String,

    // Output path
    #[clap(
        help_heading = "I/O",
        short = 'o',
        long = "out_file",
        default_value = "standard output",
        help = "Output alignment file"
    )]
    out_file: String,

    // Alignment mode
    #[clap(
        help_heading = "Alignment",
        short = 'm',
        long = "mapping-mode",
        help = "If set, use local alignment mode. If not set, use global alignment mode."
    )]
    alignment_mode: bool,

    #[clap(
        help_heading = "Alignment",
        short = 'a',
        long = "amb-strand",
        help = "If set, try aligning both the query and its reverse and complement"
    )]
    amb_mode: bool,

    // Match score
    #[clap(
        help_heading = "Alignment",
        short = 'M',
        long = "match",
        default_value_t = 1,
        help = "Match score [NOT IMPLEMENTED]"
    )]
    match_score: i32,

    // Mismatch score
    #[clap(
        help_heading = "Alignment",
        short = 'X',
        long = "mismatch",
        default_value_t = 1,
        help = "Mismatch penalty"
    )]
    mismatch_score: i32,

    // Gap open
    #[clap(
        help_heading = "Alignment",
        long = "O1",
        default_value_t = 0,
        help = "Open gap penalty 1 (for dual affine gap)"
    )]
    gap_open_1: i32,

    #[clap(
        help_heading = "Alignment",
        long = "O2",
        default_value_t = 0,
        help = "Open gap penalty 2 (for dual affine gap)"
    )]
    gap_open_2: i32,
    // Gap extension
    #[clap(
        help_heading = "Alignment",
        long = "E1",
        default_value_t = 1,
        help = "Gap extension penalty 1 (for dual affine gap)"
    )]
    gap_ext_1: i32,

    #[clap(
        help_heading = "Alignment",
        long = "E2",
        default_value_t = 1,
        help = "Gap extension penalty 2 (for dual affine gap)"
    )]
    gap_ext_2: i32,

    #[clap(
        help_heading = "Recombination",
        short = 'k',
        long = "recombinations",
        default_value_t = 0,
        help = "Number of recombinations to be performed on the graph. 0 means no recombination."
    )]
    rec_number: i32,

    //Base recombination cost
    #[clap(
        help_heading = "Recombination",
        short = 'r',
        long = "fixed-rec-cost",
        default_value_t = 4,
        help = "Recombination cost"
    )]
    base_rec_cost: i32,

    //Seed length
    #[clap(
        help_heading = "A-Star",
        short = 's',
        long = "seed-len",
        default_value_t = 8,
        help = "Seed length"
    )]
    seed_len: i32,

    #[clap(
        help_heading = "A-Star",
        short = 'e',
        long = "estimate-function",
        default_value = "chaining",
        help = "Choose the estimate function to be used"
    )]
    est_function: EstimateFunction,
}

#[derive(clap::ValueEnum, Clone, Default, Debug, Serialize, PartialEq)]
pub enum EstimateFunction {
    #[default]
    Chaining,
    Seeding,
    Fast,
}

impl Into<OsStr> for EstimateFunction {
    fn into(self) -> OsStr {
        match self {
            EstimateFunction::Chaining => OsStr::from("Chaining"),
            EstimateFunction::Seeding => OsStr::from("Seeding"),
            EstimateFunction::Fast => OsStr::from("Fast"),
        }
    }
}
#[derive(Debug)]
pub struct ScoringParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open_1: i32,
    pub gap_ext_1: i32,
    pub gap_open_2: i32,
    pub gap_ext_2: i32,
}

impl ScoringParams {
    pub fn new(
        match_score: i32,
        mismatch_score: i32,
        gap_open_1: i32,
        gap_ext_1: i32,
        gap_open_2: i32,
        gap_ext_2: i32,
    ) -> Self {
        ScoringParams {
            match_score,
            mismatch_score,
            gap_open_1,
            gap_ext_1,
            gap_open_2,
            gap_ext_2,
        }
    }
    pub fn base() -> Self {
        ScoringParams {
            match_score: 0,
            mismatch_score: 1,
            gap_open_1: 0,
            gap_ext_1: 1,
            gap_open_2: 1,
            gap_ext_2: 0,
        }
    }

    pub fn to_tuple(&self) -> (i32, i32, i32, i32, i32, i32) {
        (
            self.match_score,
            self.mismatch_score,
            self.gap_open_1,
            self.gap_ext_1,
            self.gap_open_2,
            self.gap_ext_2,
        )
    }
}
pub struct ClArgs {
    pub sequence_path: String,
    pub graph_path: String,
    pub alignment_mode: bool,
    pub scoring_params: ScoringParams,
    pub rec_number: i32,
    pub base_rec_cost: i32,
    pub out_file: String,
    pub seed_len: i32,
    pub max_rec: i32,
    pub amb_strand: bool,
    pub est_function: EstimateFunction,
}

impl ClArgs {
    pub fn parse() -> ClArgs {
        let args = Args::parse();
        ClArgs {
            sequence_path: args.sequence_path,
            graph_path: args.graph_path,
            alignment_mode: args.alignment_mode,
            scoring_params: ScoringParams::new(
                args.match_score,
                args.mismatch_score,
                args.gap_open_1,
                args.gap_ext_1,
                args.gap_open_2,
                args.gap_ext_2,
            ),
            rec_number: args.rec_number,
            base_rec_cost: args.base_rec_cost,
            out_file: args.out_file,
            seed_len: args.seed_len,
            max_rec: args.rec_number,
            amb_strand: args.amb_mode,
            est_function: args.est_function,
        }
    }
}
