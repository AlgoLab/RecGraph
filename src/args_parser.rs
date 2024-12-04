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
        //required = true,
        default_value = "example/tests/simple.gfa"
    )]
    graph_path: String,

    // Input sequence
    #[clap(
        help_heading = "I/O",
        short = 'q',
        long = "query",
        help = "Query reads (in .fasta/.fastq format)",
        //required = true,
        default_value = "example/tests/simple.fa"
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
        help = "Mismatch penalty [NOT IMPLEMENTED]"
    )]
    mismatch_score: i32,

    // Gap open
    #[clap(
        help_heading = "Alignment",
        short = 'O',
        long = "open-gap",
        default_value_t = 0,
        help = "Open gap penalty [NOT IMPLEMENTED]"
    )]
    gap_open: i32,

    // Gap extension
    #[clap(
        help_heading = "Alignment",
        short = 'E',
        long = "gap-extension",
        default_value_t = 1,
        help = "Gap extension penalty [NOT IMPLEMENTED]"
    )]
    gap_ext: i32,

    #[clap(
        help_heading = "Maximal number of recombinations allowed",
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
        default_value_t = 3,
        help = "Seed length"
    )]
    seed_len: i32,

    //Max num error per seed
    #[clap(
        help_heading = "A-Star",
        short = 'e',
        long = "err-max",
        default_value_t = 0,
        help = "Set the maximum number of errors allowed per seed between 0,1 or 2 [NOT IMPLEMENTED]"
    )]
    mex_err_seed: u8,
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
    pub base_rec_cost: i32,
    pub out_file: String,
    pub seed_len: i32,
    pub mex_err_seed: u8,
    pub max_rec: i32,
    pub amb_strand: bool
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
            base_rec_cost: args.base_rec_cost,
            out_file: args.out_file,
            seed_len: args.seed_len,
            mex_err_seed: args.mex_err_seed,
            max_rec: args.rec_number,
            amb_strand: args.amb_mode
        }
    }
}
