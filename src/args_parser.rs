use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author = "Davide Monti <d.monti11@campus.unimib.it>", version, about = "POA in rust", long_about = None)]
struct Args {
    #[clap(
        help_heading = "I/O",
        help = "Sequences to align file path",
        default_value = "sequences.fa"
    )]
    sequence_path: String,
    #[clap(
        help_heading = "I/O",
        help = "Graph file path",
        default_value = "prova.gfa"
    )]
    graph_path: String,

    #[clap(
        help_heading = "I/O",
        short = 'o',
        long = "out_file",
        default_value = "standard output",
        help = "Specifies the output file, if not indicated prints in standard output"
    )]
    out_file: String,

    // Alignment mode
    #[clap(
        help_heading = "Alignment",
        short = 'm',
        long = "aln-mode",
        default_value_t = 0,
        help = "0: global, 1: local, 2: affine gap, 3: local gap, 4: pathwise alignment[DEMO]"
    )]
    alignment_mode: i32,

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

    // Matrix type
    #[clap(
        help_heading = "Alignment",
        short = 't',
        long = "matrix",
        default_value = "none",
        help = "Scoring matrix file, if '-t' is used '-M' and '-X' are not used, you should set gap penalties in this case"
    )]
    matrix: String,

    // Gap open
    #[clap(
        help_heading = "Alignment",
        short = 'O',
        long = "gap-open",
        default_value_t = 4,
        help = "Gap opening penalty"
    )]
    gap_open: i32,

    //Gap extension
    #[clap(
        help_heading = "Alignment",
        short = 'E',
        long = "gap-ext",
        default_value_t = 2,
        help = "Gap extension penalty"
    )]
    gap_extension: i32,

    // Recombination constant multiplier
    #[clap(
        help_heading = "Alignment",
        short = 'r',
        long = "multi-rec-cost",
        default_value_t = 0.1,
        help = "Displacement multiplier"
    )]
    multi_rec_cost: f32,

    //Base recombination cost
    #[clap(
        help_heading = "Alignment",
        short = 'R',
        long = "base-rec-cost",
        default_value_t = 4,
        help = "Recombination cost, determined with -r as R + r(displacement)"
    )]
    base_rec_cost: i32,

    //Recombination band width
    #[clap(
        help_heading = "Alignment",
        short = 'B',
        long = "rec-band-width",
        default_value_t = 1.0,
        help = "Recombinatio band width"
    )]
    rec_band_width: f32,

    //Ambigous strand mode
    #[clap(
        help_heading = "Alignment",
        possible_values = &["true", "false"],
        default_value = "false",
        short = 's',
        long = "amb-strand",
        help = "Ambigous strand mode, try reverse complement if alignment score is too low"
    )]
    amb_strand: String,

    //set banding parameter, with f set the number of extra bases added (b+f*L)
    #[clap(
        help_heading = "Adaptive banded",
        default_value_t = 1,
        short = 'b',
        long = "extra-b",
        help = "First adaptive banding par, set < 0 to disable adaptive banded"
    )]
    extra_b: i32,

    #[clap(
        help_heading = "Adaptive banded",
        default_value_t = 0.01,
        short = 'f',
        long = "extra-f",
        help = "Second adaptive banding par, number of basis added to both side of the band = b+f*L, l = length of the sequence"
    )]
    extra_f: f32,
}
pub fn get_base_multi_recombination_cost() -> (i32, f32) {
    let args = Args::parse();
    (args.base_rec_cost, args.multi_rec_cost)
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

pub fn get_sequence_path() -> String {
    let args = Args::parse();
    args.sequence_path
}

pub fn get_b_f() -> (f32, f32) {
    let args = Args::parse();
    (args.extra_b as f32, args.extra_f)
}

pub fn get_graph_path() -> String {
    let args = Args::parse();
    args.graph_path
}

pub fn get_amb_strand_mode() -> bool {
    let args = Args::parse();
    let amb_strand = args.amb_strand.as_str();
    matches!(amb_strand, "true")
}

pub fn get_out_file() -> String {
    let args = Args::parse();
    args.out_file
}

pub fn get_recombination_band_width() -> f32 {
    let args = Args::parse();
    args.rec_band_width
}
