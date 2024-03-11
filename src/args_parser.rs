use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author = "Davide Monti <d.monti11@campus.unimib.it>", version, about = "RecGraph", long_about = None)]
struct Args {
    // Input graph
    #[clap(
        help_heading = "I/O",
        short = 'g',
        long = "graph",
        help = "Input graph (in .gfa format)",
        required = true
    )]
    graph_path: String,

    // Input sequence
    #[clap(
        help_heading = "I/O",
        short = 'q',
        long = "query",
        help = "Query reads (in .fasta/.fastq format)",
        required = true
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

    // Match score
    #[clap(
        help_heading = "Alignment",
        short = 'M',
        long = "match",
        default_value_t = 2,
        help = "Match score"
    )]
    match_score: i32,

    // Mismatch score
    #[clap(
        help_heading = "Alignment",
        short = 'X',
        long = "mismatch",
        default_value_t = 4,
        help = "Mismatch penalty"
    )]
    mismatch_score: i32,

    // Gap open
    #[clap(
        help_heading = "Alignment",
        short = 'O',
        long = "open-gap",
        default_value_t = 0,
        help = "Open gap penalty [NOT YET IMPLEMENTED, always 0]"
    )]
    gap_open: i32,

    // Gap extension
    #[clap(
        help_heading = "Alignment",
        short = 'E',
        long = "gap-extension",
        default_value_t = 4,
        help = "Gap extension penalty"
    )]
    gap_ext: i32,

    #[clap(
        help_heading = "Recombination",
        short = 'k',
        long = "recombinations",
        default_value_t = 0,
        help = "Number of recombinations to be performed on the graph. 0 means no recombination.[NOT YET IMPLEMENTED, always 1]"
    )]
    rec_number: i32,

    // Recombination constant multiplier
    #[clap(
        help_heading = "Recombination",
        short = 'd',
        long = "displ-multi",
        default_value_t = 0.1,
        help = "Displacement multiplier"
    )]
    multi_rec_cost: f32,

    //Base recombination cost
    #[clap(
        help_heading = "Recombination",
        short = 'r',
        long = "fixed-rec-cost",
        default_value_t = 4,
        help = "Recombination cost, determined with -d as r + d*(displacement_length)"
    )]
    base_rec_cost: i32,

    //Maximum displacement allowed
    // TODO: not yet used
    #[clap(
        help_heading = "Recombination",
        short = 'x',
        long = "max-displacement",
        default_value_t = 1,
        help = "Maximum displacement allowed between the two recombination extremities.[NOT YET IMPLEMENTED]"
    )]
    max_displacement: i32,
}

pub struct ClArgs {
    pub sequence_path: String,
    pub graph_path: String,
    pub alignment_mode: bool,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub rec_number: i32,
    pub multi_rec_cost: f32,
    pub base_rec_cost: i32,
    pub max_displacement: i32,
    pub out_file: String,
}

impl ClArgs {
    pub fn parse() -> ClArgs {
        let args = Args::parse();
        ClArgs {
            sequence_path: args.sequence_path,
            graph_path: args.graph_path,
            alignment_mode: args.alignment_mode,
            match_score: args.match_score,
            mismatch_score: -args.mismatch_score,
            gap_open: -args.gap_open,
            gap_ext: -args.gap_ext,
            rec_number: args.rec_number,
            multi_rec_cost: args.multi_rec_cost,
            base_rec_cost: args.base_rec_cost,
            max_displacement: args.max_displacement,
            out_file: args.out_file,
        }
    }
}
