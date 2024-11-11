# RecAlign
RecAlign is a sequence-to-graph aligner written in Rust. Differently from most aligners, RecGraph is an exact approach that exploits the A* algorithm for computing an **optimal** alignment between a string and a variation graph. Moreover, RecGraph can allowmultiple recombinations in the alignment in a controlled (i.e., non heuristic) way - in other words, it can perform optimal alignment to path not included in the input graphs. This follows directly from the observation that a pangenome graph includes a set of related individuals that are represented as paths of the graph.

## Installation
Install [`rust`](https://doc.rust-lang.org/cargo/getting-started/installation.html), then clone and install RecGraph:
```
git clone https://github.com/AlgoLab/RecGraph.git
cd RecGraph 
git checkout a_star
cargo build --release
```

## Usage
RecGraph requires as input a variation graph in `.gfa` format and a set of sequences (reads) in `.fasta` format and computes the alignment in `.gaf` format. To run RecGraph, run:
```
cargo run --release -q <reads.fa> -g <graph.gfa> -o <alignments.gaf>
```

#### Example
```
cargo run --release -- -m -q example/reads.fa -g example/graph.gfa -o align.gaf
```

## Parameters
RecAlign also allows to set multiple parameters to tweak the dynamic programming alignment procedure. Here the list of parameters (please check also `--help`): 
```
I/O:
    -g, --graph <GRAPH_PATH>                Path to the input graph file (in .gfa format) 
    -q, --query <SEQUENCE_PATH>             Path to the query reads file (in .fasta/.fastq format) 
    -o, --out_file <OUT_FILE>               Output alignment file [default: "standard output"]
                                
 Alignment:
    -m, --mapping-mode                      If set, perform semiglobal alignment, else perform global alignment mode.
    -s, --seed-len <SEED_LEN>               Seed length [default: 3]
    -M, --match <MATCH_SCORE>               Match score [default: 2]
    -X, --mismatch <MISMATCH_SCORE>         Mismatch penalty [default: 4]
    -E, --gap-extension <GAP_EXT>           Gap penalty [default: 4]

Recombination:
    -k, --recombinations <REC_NUMBER>       Number of recombinations to be performed on the graph. 0 means no recombination. 
                                            [default: 0]
    -r, --fixed-rec-cost <BASE_REC_COST>    Recombination cost [default: 4]
```
