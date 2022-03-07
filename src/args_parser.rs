use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Match score
    #[clap(short = 'M', long = "match")]
    match_score: i32,

    /// Mismatch score
    #[clap(short = 'X', long = "mismatch", default_value_t = 1)]
    #[structopt(allow_hyphen_values=true)]
    mismatch_score: i32,
}

pub fn get_match_mismatch() -> (i32, i32) {
    let args = Args::parse();
    (args.match_score, args.mismatch_score)
}
