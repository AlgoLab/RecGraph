# RecAlign
RecAlign is a sequence-to-graph aligner written in Rust. Differently from most aligners, RecAlign is an exact approach that exploits the A* algorithm for computing an **optimal** alignment between a string and a variation graph. Moreover, RecAlign can allowmultiple recombinations in the alignment in a controlled (i.e., non heuristic) way - in other words, it can perform optimal alignment to path not included in the input graphs. This follows directly from the observation that a pangenome graph includes a set of related individuals that are represented as paths of the graph.

## Installation
Install [`rust`](https://doc.rust-lang.org/cargo/getting-started/installation.html), then clone and install RecAlign:
```
git clone https://github.com/AlgoLab/RecGGraph.git
cd RecGraph 
git checkout a_star
cargo build --release
```
#### Static binaries
For user convenience, we provide static binaries for x86_64 linux (see [Release](https://github.com/AlgoLab/RecGraph/releases/tag/v1.0.1)).

## Usage
RecAlign requires as input a variation graph in `.gfa` format and a set of sequences (reads) in `.fasta` format and computes the alignment in `.gaf` format. To run RecAlign, run:
```
cargo run --release -q <reads.fa> -g <graph.gfa> -o <alignments.gaf>
```

#### Example
```
# if you built with cargo, from the root of this repo
cargo run --release -- -q example/reads.fa -g example/graph.gfa -o align.gaf

# if you have the precompiled binary
./recalign-x86_64 -q example/reads.fa -g example/graph.gfa -o align.gaf
```

#### Parameters
RecAlign also allows to set multiple parameters to tweak the dynamic programming alignment procedure. Here the list of parameters (please check also `--help`): 
```
I/O:
  -g, --graph <GRAPH_PATH>     Input graph (in .gfa format) [default: example/tests/simple.gfa]
  -q, --query <SEQUENCE_PATH>  Query reads (in .fasta/.fastq format) [default: example/tests/simple.fa]
  -o, --out_file <OUT_FILE>    Output alignment file [default: "standard output"]

Alignment:
  -m, --mapping-mode               If set, use local alignment mode. If not set, use global alignment mode.
  -a, --amb-strand                 If set, try aligning both the query and its reverse and complement
  -M, --match <MATCH_SCORE>        Match score [NOT IMPLEMENTED] [default: 1]
  -X, --mismatch <MISMATCH_SCORE>  Mismatch penalty [NOT IMPLEMENTED] [default: 1]
  -O, --open-gap <GAP_OPEN>        Open gap penalty [NOT IMPLEMENTED] [default: 0]
  -E, --gap-extension <GAP_EXT>    Gap extension penalty [NOT IMPLEMENTED] [default: 1]
  -t, --threads <THREADS>          Number of threads to use 

Recombination:
  -k, --recombinations <REC_NUMBER>     Number of recombinations to be performed on the graph. 0 means no recombination. [default: 0]
  -r, --fixed-rec-cost <BASE_REC_COST>  Recombination cost [default: 4]

A-Star:
  -s, --seed-len <SEED_LEN>
          Seed length [default: 8]
  -e, --estimate-function <EST_FUNCTION>
          Choose the estimate function to be used [default: chaining] [possible values: chaining, seeding, fast]
```
