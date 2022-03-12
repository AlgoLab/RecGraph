use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author = "Davide Monti <d.monti11@campus.unimib.it>", version, about = "POA in rust", long_about = None)]
struct Args {
    // Alignment mode
    #[clap(
        short = 'm',
        long = "aln-mode",
        default_value_t = 0,
        help = "0: global, 1: local, 2: edit distance, 3: affine gap"
    )]
    alignment_mode: i32,

    /// Match score
    #[clap(short = 'M', long = "match", default_value_t = 2, help = "Match score")]
    match_score: i32,

    /// Mismatch score
    #[clap(
        short = 'X',
        long = "mismatch",
        default_value_t = 4,
        help = "Mismatch penalty"
    )]
    mismatch_score: i32,

    // Matrix type
    #[clap(
        short = 't',
        long = "matrix",
        default_value = "none",
        help = "Scoring matrix file, if '-t' is used '-M' and '-X' are not used"
    )]
    matrix: String,

    // Gap open
    #[clap(
        short = 'O',
        long = "gap-open",
        default_value_t = 4,
        help = "Gap opening penalty"
    )]
    gap_open: i32,

    //Gap extension
    #[clap(
        short = 'E',
        long = "gap-ext",
        default_value_t = 2,
        help = "Gap extension penalty"
    )]
    gap_extension: i32,
}

pub fn get_match_mismatch() -> (i32, i32) {
    let args = Args::parse();
    (args.match_score, -args.mismatch_score)
}

pub fn get_matrix_type() -> String {
    let args = Args::parse();
    args.matrix
}

pub fn get_gap_open_gap_ext() -> (i32, i32) {
    let args = Args::parse();
    (-args.gap_open, -args.gap_extension)
}

pub fn get_align_mode() -> i32 {
    let args = Args::parse();
    args.alignment_mode
}
