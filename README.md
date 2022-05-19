# rspoa
POA in Rust
## Usage
After cloning the repository from GitHub, move to its main directory and compile the package with
```
cargo build --release
```
To use rspoa you need a .gfa file and a .fasta file containing the reads you want to align to the graph. In order to run rspoa:
```
cargo run --release <reads.fa> <graph.gfa>
```
You can set different parameters and choose between different alignment types:
```
OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

I/O:
    <SEQUENCE_PATH>          Sequences to align file path 
    <GRAPH_PATH>             Graph file path 
     -o, --out_file <OUT_FILE>    Specifies the output file, if not indicated prints in standard
                                 output [default: "standard output"]

Alignment:
    -m, --aln-mode <ALIGNMENT_MODE>    0: global, 1: local, 2: affine gap, 3: local gap, 
                                       4: pathwise alignment[DEMO], 
    -M, --match <MATCH_SCORE>          Match score [default: 2]
    -X, --mismatch <MISMATCH_SCORE>    Mismatch penalty [default: 4]
    -O, --gap-open <GAP_OPEN>          Gap opening penalty [default: 4]
    -E, --gap-ext <GAP_EXTENSION>      Gap extension penalty [default: 2]
    -s, --amb-strand <AMB_STRAND>      Ambigous strand mode, try reverse complement if
                                       alignment score is too low [default: false] [possible values:
                                       true, false]
    -t, --matrix <MATRIX>              Scoring matrix file, if '-t' is used '-M' and '-X' are not
                                       used. You should set appropriate gap penalties if aln-mode is 
                                       2 or 3[default: none]

Adaptive banded:
    -b, --extra-b <EXTRA_B>    First adaptive banding parameter, set < 0 to disable adaptive banded
                               [default: 1]
    -f, --extra-f <EXTRA_F>    Second adaptive banding parameter, number of basis added to both side of
                               the band = b+f*L, l = length of the sequence [default: 0.01]
```

The output is in .gaf format, in reference to the graph in input.

```
>read1	50	0	50	+	>1>3>5>6>8>9>11>12>13>15>16>18>19	50	0	11	49	*	*	8M,1M,1M,3M,1M,4M,5M,2M,8M,1M,4M,1M,11M
```
### Library
rsPOA can also be used as a library for your project. 
In order to do so add to your Cargo.toml file:
```
[dependencies]
rspoa= { git = "https://github.com/AlgoLab/rspoa" }
```
In this case you can use the functions inside the api.rs file, they need only the read as a string and the graph as an HashGraph, other parameters can be set or left by default.
You can use them by adding ```use rspoa::api::*``` in your file.
