# RecGraph
RecGraph is a sequence-to-graph aligner written in Rust. Differently from most aligners, RecGraph is an exact approach that implements a dynamic programming algorithm for computing an **optimal** alignment between a string and a variation graph. Moreover, RecGraph can allow recombinations in the alignment in a controlled (i.e., non heuristic) way - in other words, it can perform optimal alignment to path not included in the input graphs. This follows directly from the observation that a pangenome graph includes a set of related individuals that are represented as paths of the graph.

## Installation
We support 4 different ways to obtain and use RecGraph:
* [downloading static binaries](#static-binaries)
* [building via cargo](#compilation)
* [installing via conda](#installation-from-conda)
* [getting a docker image](#docker-image)

#### Static binaries
For user convenience, we provide static binaries for x86_64 linux and windows systems (see [Releases](https://github.com/AlgoLab/RecGraph/releases)).

#### Compilation
Install [`rust`](https://doc.rust-lang.org/cargo/getting-started/installation.html), then clone and install RecGraph:
```
git clone https://github.com/AlgoLab/RecGraph.git
cd RecGraph
cargo build --release
```

#### Installation from conda
RecGraph is available on bioconda:
```
conda create -n recgraph -c conda-forge -c bioconda recgraph
```

#### Docker image
We provide a docker image, hosted on [docker.hub](https://hub.docker.com/r/algolab/recgraph):
```
docker pull algolab/recgraph
docker run algolab/recgraph --help
```

## Usage
RecGraph requires as input a variation graph in `.gfa` format and a set of sequences (reads) in `.fasta` format and computes the alignment in `.gaf` format. To run RecGraph, run:
```
cargo run --release <reads.fa> <graph.gfa> > <alignments.gaf>
```

#### Example
```
# if you built with cargo, from the root of this repo
cargo run --release -- -m 1 example/reads.fa example/graph.gfa > align.gaf

# if you have the precompiled binary
./recgraph_linux_x86-64 -m 1 example/reads.fa example/graph.gfa > align.gaf

# if you have the conda version, within the correct environment
recgraph -m 1 example/reads.fa example/graph.gfa > align.gaf

# if you use docker, please bind the volume into the container (-v)
docker run -v $(pwd)/example:/data algolab/recgraph -m1 reads.fa graph.gfa > align.gaf
```

### Alignment modes
RecGraph can be run in several different modes (`-m` flag):
* `-m [0,1,2,3]` performs the classical POA (global, local, affine gap, and local gap)
* `-m [4,5]` performs global/semiglobal alignment in pathwise mode (i.e., following the paths of the graph)
* `-m [8,9]` performs global/semiglobal alignment in recombination mode (i.e., allowing weighted recombinations)

`-m 6` and `-m 7` are experimental and are not fully tested yet. They perform global/semiglobal alignment with affine gap in pathwise mode.

### Other parameters
RecGraph also allows to set multiple parameters to tweak the dynamic programming alignment procedure. Here the list of parameters (please check also `--help`): 
```
    -M, --match <MATCH_SCORE>                Match score [default: 2]
    -X, --mismatch <MISMATCH_SCORE>          Mismatch penalty [default: 4]
    -O, --gap-open <GAP_OPEN>                Gap opening penalty [default: 4]
    -E, --gap-ext <GAP_EXTENSION>            Gap extension penalty [default: 2]
    -R, --base-rec-cost <BASE_REC_COST>      Recombination cost,
                                             determined with -r as R + r*(displacement_length) [default: 4]
    -r, --multi-rec-cost <MULTI_REC_COST>    Displacement multiplier [default: 0.1]
    -B, --rec-band-width <REC_BAND_WIDTH>    Recombination band width [default: 1]
    -b, --extra-b <EXTRA_B>                  First adaptive banding par,
                                             set < 0 to disable adaptive banded [default: 1]
    -f, --extra-f <EXTRA_F>                  Second adaptive banding par, number of basis added to both side of
                                             the band = b+f*L, l = length of the sequence [default: 0.01]
    -t, --matrix <MATRIX>                    Scoring matrix file, if '-t' is used, '-M' and '-X' are not used
                                             and you should set gap penalties in this case [default: none]
```

### Library
RecGraph can also be used as a library for your project. To do so, add these lines to your `Cargo.toml`:
```
[dependencies]
RecGraph = { git = "https://github.com/AlgoLab/RecGraph" }
```
You can use the functions defined in the [`api.rs`](https://github.com/AlgoLab/RecGraph/blob/1b513973c1145015ed626abc975e276970d2a60e/src/api.rs) file (e.g., by adding `use RecGraph::api::*` to your file). All the functions require just a read (as a string) and the graph (as an HashGraph). Other parameters are optional.
