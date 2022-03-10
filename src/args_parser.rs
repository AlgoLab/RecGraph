use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Match score
    #[clap(short = 'M', long = "match", default_value_t = 1)]
    match_score: i32,

    /// Mismatch score
    #[clap(short = 'X', long = "mismatch", default_value_t = -1)]
    #[structopt(allow_hyphen_values = true)]
    mismatch_score: i32,

    // Matrix type
    #[clap(short = 't', long = "matrix", default_value = "none")]
    matrix: String,

    // Alignment mode
    #[clap(short = 'm', long = "aln-mode", default_value_t = 0)]
    alignment_mode: i32,
}

pub fn get_match_mismatch() -> (i32, i32) {
    let args = Args::parse();
    (args.match_score, args.mismatch_score)
}

pub fn get_matrix_type() -> String {
    let args = Args::parse();
    args.matrix
}

pub fn get_align_mode() -> i32 {
    let args = Args::parse();
    args.alignment_mode
}
